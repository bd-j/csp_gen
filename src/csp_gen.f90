subroutine csp_gen(mass_ssp, lbol_ssp, spec_ssp, pset, tage,&
                   mass_csp, lbol_csp, spec_csp)
  !
  ! Return the spectrum (and mass and lbol) of a composite stellar population.
  !
  ! Inputs
  ! ---------
  !
  ! mass_ssp, lbol_ssp, spec_ssp:
  !   The (surviving) stellar masses, bolometric luminosities, and spectra of
  !   the SSPs.
  !
  ! pset:
  !   A `PARAMS` structure containing the SFH parameters
  !
  ! tage:
  !   The age (in Gyr) at which the spectrum is desired.  Note that this can be
  !   different than pset%tage if the latter is 0.
  !
  ! Outputs
  ! ---------
  !
  ! mass_csp, lbol_csp, spec_csp:
  !   The (surviving) stellar masses, bolometric luminosity, and spectrum of
  !   the composite stellar population at tage, normalized to 1 M_sun *formed*.
  
  use sps_vars, only: ntfull, nspec, time_full, tiny_number, sfhtab
  use sps_utils, only: locate, sfh_weight
  implicit none

  real(SP), intent(in), dimension(ntfull) :: mass_ssp, lbol_ssp
  real(SP), intent(in), dimension(nspec, ntfull) :: spec_ssp
  type(PARAMS), intent(in) :: pset
  real(SP), intent(in) :: tage
  
  real(SP), intent(out) :: mass_csp, csp_lbol
  real(SP), intent(in), dimension(nspec) :: spec_csp


  real(SP), dimension(ntfull) :: total_weights=0., w1=0., w2=0.
  integer :: j, imin=0, imax=ntfull
  type(SFHPARAMS) :: sfhpars
  real(SP) :: mass, m1, m2
  
  ! Build a structure containing useful units, numbers, and switches for the
  ! weight calculations.
  call convert_sfhparams(pset, tage, sfhpars)
     
  ! ----- Get SFH weights -----

  ! SSP.
  if (pset%sfh.eq.0) then
     ! Make sure to use SSP weighting scheme
     sfhpars%type = -1
     ! Use tage as the burst lookback time, instead of tage-tburst.
     sfhpars%tb = sfhpars%tage
     imin = min(max(locate(time_full, log10(sfhpars%tage)), 0), ntfull)
     imax = min(imin+1, ntfull)
     ! These come pre-normalized to 1 Msun
     total_weights = sfh_weight(sfhpars, imin, imax)
  endif
  
  ! Tau and delayed-tau.
  if ((pset%sfh.eq.1).or.(pset%sfh.eq.4)) then
     sfhpars%type = pset%sfh
     ! Only calculate SFH weights for SSPs up to tage (plus the next one).
     imax = min(max(locate(time_full, log10(sfhpars%tage)) + 2, 1), ntfull)
     total_weights = sfh_weight(sfhpars, 0, imax)
     ! Could save some loops by having proper normalization analytically from sfh_weight
     m1 = sum(total_weights)
     if (m1.lt.tiny_number) m1 = 1.0
     total_weights = total_weights / m1
  endif
  
  ! Add constant and burst weights for SFH=1,4
  if (((pset%sfh.eq.1).or.(pset%sfh.eq.4)).and.&
       ((pset%const.gt.0).or.(pset%fburst.gt.tiny_number))) then
     ! Constant
     sfhpars%type = 0
     w1 = sfh_weight(sfhpars, 0, ntfull)
     m1 = sum(w1)
     ! Burst.  These weights come pre-normalized to 1 Msun.
     sfhpars%type = -1
     w2 = sfh_weight(sfhpars, 0, ntfull)
     ! Sum with proper relative normalization.  Beware divide by zero.
     if (m1.lt.tiny_number) m1 = 1.0
     total_weights = (1 - pset%const - pset%fburst) * total_weights + &
                      pset%const * (w1 / m1) + &
                      pset%fburst * w2
  endif

  ! Simha
  if (pset%sfh.eq.5) then
     ! Only calculate SFH weights for SSPs up to tage (plus the next one).
     imax = min(max(locate(time_full, log10(sfhpars%tage)) + 2, 1), ntfull)
     ! Delayed-tau model portion.
     sfhpars%type = 4
     w1 = sfh_weight(sfhpars, 0, imax)
     m1 = sum(w1)
     ! Linear portion.  Need to set use_simha_limits flag to get correct limits.
     sfhpars%type = 5
     sfhpars%use_simha_limits = 1
     w2 = sfh_weight(sfhpars, 0, imax)
     sfhpars%use_simha_limits = 0
     m2 = sum(w2)
     ! Normalize and sum.  need to be careful of divide by zero here, if all
     ! linear or all delay-tau s.t. weights sum to zero for a component.
     if (m1.lt.tiny_number) m1 = 1.0
     if (m2.lt.tiny_number) m2 = 1.0
     total_weights = w1 / m1 + simha_norm(pset) * (w2 / m2)
  endif


  ! Tabular.  Time units in sfhtab are assumed to be linear years of lookback time.
  if (pset%sfh.eq.2.or.pset%sfh.eq.3) then
     total_weights = 0.
     call setup_tabular()
     ! Linearly interpolate in the bins.
     sfhpars%type = 5
     ! Loop over each bin.
     do i=1,ntabsfh-1
        ! mass formed in this bin assuming linear
        mass = (sfhtab(i,2) + sfhtab(i+1, 2)) * (sfhtab(i+1,1) - sfhtab(i, 1)) / 2
        ! min and max ssps to consider
        imin = min(max(locate(time_full, log10(sfhtab(i, 1))) - 1, 0), ntfull)
        imax = min(max(locate(time_full, log10(sfhtab(i+1, 1))) + 2, 0), ntfull)
        ! set integration limits
        sfhpars%tq = sfhtab(i, 1)
        sfhpars%tage = sfhtab(i+1, 1)
        sfhpars%sf_slope = (sfhtab(i, 2) - sfhtab(i+1, 2)) / (sfhtab(i+1, 1) - sfhtab(i, 1))
        
        ! get the weights for this bin in the tabulated sfh and add to the
        ! total weight, after normalizing
        w1 = sfh_weight(sfh_pars, imin, imax)
        m1 = sum(w1)
        if (m1.lt.tiny_number) m1 = 1.0
        total_weight = total_weight + w1 * (mass / m1)
     enddo
     imax = ntfull
  endif

  ! Now weight each SSP by `total_weight` and sum.
  ! This matrix multiply could probably be optimized!!!!
  do j=0, imax
     if (total_weight(j).gt.tiny_number) then
        spec_csp = spec_csp + total_weight(j) * spec_ssp(:, j)
     endif
  enddo
  mass_csp = sum(mass_ssp * total_weight)
  lbol_csp = sum(lbol_ssp * total_weight)
  
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
  !       - `tb` is the burst time, in lookback time.
  !
  use sps_vars, only: pset
  implicit none
  
  type(PARAMS), intent(in) :: pset
  real(SP), intent(in) :: tage
  
  type(SFHPARAMS), intent(inout) :: sfh

  ! Convert units from Gyr to yr
  sfh%tage = tage * 1e9
  sfh%tau = pset%tau * 1e9
  sfh%tburst = pset%tburst * 1e9
  sfh%sf_trunc = pset%sf_trunc * 1e9
  ! Note the sign flip here!
  sfh%sf_slope = -pset%sf_slope / 1e9

  ! convert tburst to lookback time
  sfh%tb = sfh%tage - sfh%tburst
  
  ! convert sf_trunc to to lookback time
  if ((sfh%sf_trunc.le.0).or.(sfh%sf_trunc.gt.sfh%tage)) then
     sfh%tq = 0.
  else
     sfh%tq = sfh%tage - sfh%sf_trunc
  endif
  ! For simha get zero crossing time (in lookback time), avoiding divison by zero
  ! Note that only positive slopes have a chance to hit zero SFR.
  if (sfh%sf_slope.gt.tiny_number) then
     sfh%t0 = sfh%tq - 1. / sfh%sf_slope
  else
     sfh%t0 = 0.
  endif
  ! If the zero crossing time is outside the linear regime, set it to zero.
  if ((sfh%t0.gt.sfh%tq).or.(sfh%t0.le.0)) then
     sfh%t0 = 0.
  endif

end subroutine convert_sfhparams


function simha_norm(pset)
  ! Get the relative normalization of the linear and delyed tau portions of the
  ! simha SFH
  !
  
  implicit none
  type(PARAMS), intent(in) :: pset

  real(SP) :: simha_norm

  real(SP) :: Tmax, Tz
  real(SP) :: mass_tau, mass_linear
  real(SP) :: m, sfr_q

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
  simha_norm = (sfr_q * mass_linear) / mass_tau
  
end function simha_norm
