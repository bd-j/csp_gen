subroutine csp_gen(mass_ssp,lbol_ssp,spec_ssp,pset,tage,csp_spec)
  !
  ! Return the spectrum of a composite stellar population
  !
  ! sfh=1: tau model
  ! sfh=2: tabulated SFH (from file)
  ! sfh=3: tabulated SFH (stored in sfhtab arr)
  ! sfh=4: delayed tau model
  ! sfh=5: custom SFH (see Simha et al. 2013)

  use sps_vars, only: ntfull
  implicit none
  
  ! Build a structure containing useful units, numbers, and switches for the
  ! weight calculations
  call convert_sfhparams(pset,tage,sfhpars)
     
  ! ----- Get SFH weights -----
  
  ! Tau and delayed-tau
  if ((pset%sfh.EQ.1).OR.(pset%sfh.EQ.4)) THEN
     sfhpars%type = pset%sfh
     ! only calculate SFH weights for SSPs up to tage (plus the next one)
     imax = clip(locate(time_full, log10(sfhpars%tage)) + 1, 1, ntfull)
     total_weights = sfh_weight(sfhpars, 1, imax)
     total_weights = total_weights / sum(total_weights)
  endif

  ! Simha
  if (pset%sfh.eq.5) then
     ! only calculate SFH weights for SSPs up to tage (plus the next one)
     imax = clip(locate(time_full, log10(sfhpars%tage)) + 1, 1, ntfull)
     ! delayed-tau model portion
     sfhpars%type = 4
     w1 = sfh_weight(sfhpars, 1, imax)
     ! linear portion.  Need to set use_simha_limits flag to get corect limits
     sfhpars%type = 5
     sfhpars%use_simha_limits = 1
     w2 = sfh_weight(sfhpars, 1, imax)
     sfhpars%use_simha_limits = 0
     ! normalize and sum.  need to be careful of divide by zero here, if all
     ! linear or all delay-tau
     total_weights = w1 / sum(w1) + simha_norm(pset) * (w2 / sum(w2))
  endif

  ! Add constant and burst weights
  if ((pset%const.gt.0).or.(pset%fburst.gt.tiny_number)) then
     ! Constant
     sfhpars%type = 0
     w1 = sfh_weight(sfhpars, 1, ntfull)
     ! burst.  These weights come pre-normalized to 1 Msun
     w2 = burst_weight(pset%tburst)
     ! sum with proper relative normalization
     total_weights = (1 - pset%const - pset%fburst) * total_weights + &
          pset%const * (weights1 / sum(weights1)) + &
          pset%fburst * weights2 )
  endif

  ! Tabular
  if (pset%sfh.EQ.2.OR.pset%sfh.EQ.3) then
     total_weights = 0.
     call setup_tabular()
     ! linearly interpolate in the bins
     sfhpars%type = 5
     ! loop over each bin
     do i=1,ntabsfh-1
        ! mass formed in this bin assuming linear
        mass = (sfhtab(i,2) + sfhtab(i+1, 2)) / (2 * (sfhtab(i+1,1) - sfhtab(i, 1)))
        ! min and max ssps to consider
        imin = clip(locate(time_full, log10(sfhtab(i, 1)) - 1, 1, ntfull)
        imax = clip(locate(time_full, log10(sfhtab(i+1, 1)) + 1, 1, ntfull)
        ! set integration limits
        sfhpars%tq = sfhtab(i, 1)
        sfhpars%tage = sfhtab(i+1, 1)
        ! get the weights for this bin in the tabulated sfh and add to the
        ! total weight, after normalizing
        w = sfh_weight(sfh_pars, imin, imax)
        total_weight = total_weight + w / sum(w) * mass
     enddo
  endif

  ! Now weight each SSP by `total_weight` and sum
  csp_spec =
  csp_mass =
  csp_lbol = 
  
end subroutine csp_gen

subroutine convert_sfhparams(pset, tage, sfh)
  ! Convert the pset values to yrs and pre-calculate some useful lookback
  ! times, and store in a structure.  Note that this subroutine uses but does
  ! not alter the pset values.
  !
  ! Inputs
  ! ----------
  !
  ! pset:
  !    The parameter set.
  !
  ! tage:
  !    The age (in forward time, not lookback time) at which you are
  !    calculating the composite spectrum, in Gyr.
  !
  ! Outputs
  ! ----------
  ! sfh:
  !    An SFHPARAMS structure, containing SFH parameters in yrs, and some
  !    useful lookback times.
  !       - `tq` is the truncation time, in lookback time
  !       - `t0` is the zero crossing time for a linear SFH, in lookback time.
  !
  use sps_vars, only: pset
  implicit none
  
  type(PARAMS), intent(in) :: pset
  real(SP), intent(in) :: tage
  
  type(SFHPARAMS), intent(inout) :: sfh

  ! Convert units from Gyr to yr
  sfh%tage = tage * 1e9
  sfh%tau = pset%tau * 1e9
  sfh%sf_trunc = pset%sf_trunc * 1e9
  sfh%sf_slope = pset%sf_slope / 1e9

  ! convert sf_trunc to to t_lookback
  if ((sfh%sf_trunc.le.0).or.(sfh%sf_trunc.gt.sfh%tage)) then
     sfh%tq = 0
  else
     sfh%tq = sfh%tage - sfh%sf_trunc
  endif
  ! get zero crossing time (in lookback time), avoiding divison by zero
  if (abs(sfh%sf_slope).gt.tiny_number) then
     sfh%t0 = sfh%tq - 1. / sfh%sf_slope
  else
     sfh%t0 = 0
  endif
  ! If the zero crossing time is outside the linear regime, set it to zero.
  if ((sfh%t0.gt.sfh%tq).or.(sfh%t0.le.0)) then
     sfh%t0 = 0
  endif

end subroutine convert_sfhparams


function simha_norm(pset)
  ! Get the relative normalization of the linear and delyed tau portions of the
  ! simha SFH
  !
  
  implicit none
  type(PARAMS), intent(in) :: pset

  real(SP) :: Tmax, Tz
  real(SP) :: mass_tau, mass_linear
  real(SP) :: m, sfr_q

  real(SP) intent(out) :: simha_norm

  m = pset%sf_slope
  
  ! mass in the delayed tau portion
  if ((pset%sf_trunc.le.0).or.(pset%sf_trunc.gt.pset%tage)) then
     Tmax = pset%tage
  else
     Tmax = pset%sf_trunc
  endif
  mass_tau = pset%tau * gammainc(2, Tmax/pset%tau)

  ! sfr at sf_trunc
  sfr_q = (Tmax/pset%tau) * exp(-Tmax/pset%tau)
  
  ! Mass in the linear portion.
  ! This is integral of (1 - m * (T - Tmax)) from Tmax to Tzero
  if (m.gt.0) then
     Tz = Tmax + 1.0/m
  else
     Tz = pset%tage
  endif
  if ((Tz.lt.Tmax).or.(Tz.gt.pset%tage)) then
     Tz = pset%tage
  endif
  mass_linear = (Tz - Tmax) - m/2.*(Tz**2 + Tmax**2) + m*Tz*Tmax
  
  ! now compute normalization relative to 1 Msun in the delayed-tau portion
  simha_norm = sfr_q * (mass_linear / mass_tau)
  
end function simha_norm
