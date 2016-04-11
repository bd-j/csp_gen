subroutine csp_gen(write_compsp,nzin,outfile,mass_ssp,&
     lbol_ssp,tspec_ssp,pset,ocompsp)
  !
  ! Return the spectrum of a composite stellar population
  !
  ! sfh=1: tau model
  ! sfh=2: tabulated SFH (from file)
  ! sfh=3: tabulated SFH (stored in sfhtab arr)
  ! sfh=4: delayed tau model
  ! sfh=5: custom SFH (see Simha et al. 2013)

  use sps_vars
  implicit none

  spec_ssp = tspec_ssp
  
  ! ------ Prepare SSPs ------
  ! Dust attenuation + emission, nebular emission etc.

  ! Add nebular emission
  IF (add_neb_emission.EQ.1) THEN
     IF (nzin.GT.1) THEN
        WRITE(*,*) 'COMPSP ERROR: cannot handle both nebular '//&
             'emission and mult-metallicity SSPs in compsp'
        STOP
     ENDIF
     CALL ADD_NEBULAR(pset,tspec_ssp(:,:,1),spec_ssp(:,:,1))
  ENDIF

  ! Add dust emission
  ! This is a wierd way to do this
  if ((pset%dust1.gt.tiny_number).or.(pset%dust2.gt.tiny_number)) then
     do i=1, ntfull
        csp1 = 0.0
        csp2 = 0.0
        IF (time_full(i).LT.pset%dust_tesc) THEN
           csp1 = spec_ssp(:,i,1)
        ELSE
           csp2 = spec_ssp(:,i,1)
        ENDIF
        !add dust and combine young and old csp
        CALL ADD_DUST(pset,csp1,csp2,spec_csp,mdust)
        spec_ssp(:,i,1) = spec_csp
     enddo
  endif
     
  ! ----- Get SFH weights -----

  ! Build a structure containing useful units and switches
  call convert_sfhparams(pset, sfhpars)
  
  ! Tau and delayed-tau
  if ((pset%sfh.EQ.1).OR.(pset%sfh.EQ.4)) THEN
     sfhpars%type = pset%sfh
     total_weights = sfh_weight(sfhpars, 1, ntfull)
     total_weights = total_weights / sum(total_weights)
  endif

  ! Simha
  if (pset%sfh.eq.5) then
     ! delayed-tau model portion
     sfhpars%type = 4
     w1 = sfh_weight(sfhpars, 1, ntfull)
     ! linear portion.  Need to set use_simha_limits flag to get corect limits
     sfhpars%type = 5
     sfhpars%use_simha_limits = 1
     w2 = sfh_weight(sfhpars, 1, ntfull)
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
     call setup_tabular()
     ! linearly interpolate in the bins
     sfhpars%type = 5
     ! loop over each bin
     do i=1,ntabsfh-1
        ! mass formed in this bin assuming linear
        mass = (sfhtab(i,2) + sfhtab(i+1, 2)) / (2 * (sfhtab(i+1,1) - sfhtab(i, 1)))
        ! min and max ssps to consider
        imin = locate(time_full, log10(sfhtab(i, 1))
        imax = locate(time_full, log10(sfhtab(i+1, 1))
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

  ! Add IGM absorption
  
end subroutine csp_gen


function simha_norm(pset)
  ! mass in the delayed tau portion
  if ((pset%sf_trunc.le.0).or.(pset%sf_trunc.gt.tage)) then
     Tmax = pset%tage
  else
     Tmax = pset%sf_trunc
  endif
  mass_tau = pset%tau * gammainc(2, Tmax/pset%tau)

  ! sfr at sf_trunc
  sfr_q = (Tmax/pset%tau) * exp(-Tmax/pset%tau)
  
  ! Mass in the linear portion.
  ! This is integral of (1 - m * (T - Tmax)) from Tmax to Tzero
  Tz = Tmax + 1/np.float64(sf_slope) ! need to deal with divide by zero
  if ((Tz.lt.Tmax).or.(Tz.gt.pset%tage).or.(pset%sf_slope.eq.0)) then
     Tz = pset%tage
  endif
  m = pset%sf_slope
  mass_linear = (Tz - Tmax) - m/2.*(Tz**2 + Tmax**2) + m*Tz*Tmax
  
  ! now compute normalization relative to 1 Msun in the delayed-tau portion
  simha_norm = mass_linear * sfr_q / mass_tau
  
end function simha_norm
