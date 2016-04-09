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

  ! Tau and delayed-tau
  if ((pset%sfh.EQ.1).OR.(pset%sfh.EQ.4)) THEN
     total_weights = sfh_weight(pset%sfh, 1, ntfull)
     total_weights = total_weights / sum(total_weights)
  endif

  ! Simha
  if (pset%sfh.eq.5) then
     ! tau model portion
     weights1 = sfh_weight(4, 1, ntfull)
     ! linear portion.  Need to set simha_linear flag to get corect limits
     weights2 = sfh_weight(5, 1, ntfull, 1)
     ! normalize and sum
     norm = 
     total_weights = weights1 / sum(weights1) + norm * (weights2 / sum(weights2))
  endif

  ! Add constant and burst weights
  if ((pset%const.gt.0).or.(pset%fburst.gt.tiny_number)) then
     ! Constant
     weights1 = sfh_weight(0, 1, ntfull)
     ! burst.  These weights come pre-normalized to 1 Msun
     weights2 = burst_weight(pset%tburst)
     ! sum with proper relative normalization
     total_weights = (1 - pset%const - pset%fburst) * total_weights + &
          pset%const * (weights1 / sum(weights1)) + &
          pset%fburst * weights2 )
  endif

  ! Tabular
  if (pset%sfh.EQ.2.OR.pset%sfh.EQ.3) then
     call setup_tabular()
  endif

end subroutine csp_gen
