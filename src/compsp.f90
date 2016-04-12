SUBROUTINE COMPSP(write_compsp,nzin,outfile,&
     mass_ssp,lbol_ssp,tspec_ssp,&
     pset,ocompsp)

  INTEGER, INTENT(in) :: write_compsp,nzin
  REAL(SP), INTENT(in), DIMENSION(ntfull,nzin) :: lbol_ssp,mass_ssp
  REAL(SP), INTENT(in), DIMENSION(nspec,ntfull,nzin) :: tspec_ssp
  REAL(SP), DIMENSION(nspec,ntfull,nzin) :: spec_ssp
  CHARACTER(100), INTENT(in) :: outfile

  ! ------ Various checks and setup ------
  
  IF (check_sps_setup.EQ.0) THEN
     WRITE(*,*) 'COMPSP ERROR: '//&
          'SPS_SETUP must be run before calling COMPSP. '
     STOP
  ENDIF
  !make sure various variables are set correctly
  CALL COMPSP_WARNING(maxtime,pset,nzin,write_compsp)
  !setup output files
  if (pset%tage.gt.0) then
     nage = 1
  else
     nage = ntfull
  endif
  IF (write_compsp.GT.0) &
       CALL COMPSP_SETUP_OUTPUT(write_compsp,pset,outfile,1,nage)

  ! Isochrone case just writes the CMDs and exits
  IF (write_compsp.EQ.5) THEN
     CALL WRITE_ISOCHRONE(outfile,pset)
     RETURN
  ENDIF

  ! ------ Prepare SSPs ------
  ! Dust attenuation + emission, nebular emission etc.
  ! We should probably only do this for ages up to tage, if it is set.
  ! Also we will operate on copies of the spectra

  spec_ssp = tspec_ssp
  
  ! Add nebular emission
  if (add_neb_emission.EQ.1) then
     if (nzin.GT.1) then
        WRITE(*,*) 'COMPSP ERROR: cannot handle both nebular '//&
             'emission and mult-metallicity SSPs in compsp'
        STOP
     endif
     call add_nebular(pset,tspec_ssp(:,:,1),spec_ssp(:,:,1))
  endif

  ! Add dust emission
  !
  ! This is a wierd way to do this. Also not strictly correct, since the age
  ! bin older than dust_tesc will include some contribution from young star
  ! dust due to the interpolation, and changing dust_tesc by values smaller
  ! than the ssp age grid will have no effect on the output.
  !
  ! The correct way is seprately calculate csp_spectra for the t<tesc and t >
  ! tesc portions and feed to add_dust as intended.  However, this requires
  ! significant extra logic (basically messing around with the integration
  ! limits in sfh_weight) for probably little gain in accuracy.
  !
  ! We could probably also do a little better by attenuating up to the SSP
  ! *nearest* in age to dust_tesc, instead of the oldest SSP still younger than
  ! dust_tesc
  if ((pset%dust1.gt.tiny_number).or.(pset%dust2.gt.tiny_number)) then
     do i=1, ntfull
        csp1 = 0.0
        csp2 = 0.0
        if (time_full(i).LT.pset%dust_tesc) then
           csp1 = spec_ssp(:,i,1)
        else
           csp2 = spec_ssp(:,i,1)
        endif
        !add dust and combine young and old csp
        call add_dust(pset,csp1,csp2,spec_dusty,mdust)
        spec_ssp(:,i,1) = spec_dusty
        mdust_ssp(i) = mdust
     enddo
  endif

  ! --- Get CSP spectra -------
  
  ! loop over output ages
  do i=1,ntfull
     if (pset%tage.gt.0) then
        ! A specific age was asked for, we will only compute one spectrum
        age = pset%tage
     else
        ! Otherwise we will calculate composite spectra for each SSP age.
        age = 10**(time_full(i)-9)
     endif

     ! -----
     ! Get the spectrum for this age.  Note this is always normalized to one
     ! solar mass formed, so we actually need to renormalize if computing all ages
     call csp_gen(mass_ssp,lbol_ssp,spec_ssp,mdust_ssp,&
          pset,age,&
          mass_csp,lbol_csp,spec_csp,mdust_csp)     

     call sfhinfo(pset, age, mass_frac, tsfr)
     
     ! -------
     ! Now do a bunch of stuff with the spectrum
     ! Smooth the spectrum
     if (pset%sigma_smooth.GT.0.0) then
        call smoothspec(spec_lambda,spec_csp,pset%sigma_smooth,&
             pset%min_wave_smooth,pset%max_wave_smooth)
     endif
     ! Add IGM absorption
     if (add_igm_absorption.EQ.1.AND.pset%zred.GT.tiny_number) then
        spec_csp = igm_absorb(spec_lambda,spec_csp,pset%zred,&
             pset%igm_factor)
     endif
     ! Compute spectral indices
     if (write_compsp.EQ.4) then
        call getindx(spec_lambda,spec_csp,indx)
     else
        indx=0.0
     endif
     ! Compute mags
     if (redshift_colors.EQ.0) then
        call getmags(pset%zred,spec_csp,mags,pset%mag_compute)
     else
        ! here we compute the redshift at the corresponding age
        zred = min(max(linterp(cosmospl(:,2),cosmospl(:,1),&
             10**time_full(i)/1E9),0.0),20.0)
        call getmags(zred,spec_csp,mags,pset%mag_compute)
     endif

     ! ---------
     ! Store the spectrum and write....
     call save_compsp(write_compsp,ocompsp(i),log10(age)+9,&
          mass_csp,lbol_csp,mags,tsfr,spec_csp,mdust_csp,indx)

     ! Terminate the loop if specific age was requested
     if (pset%tage.gt.0) then
        exit
     endif
     
  enddo

  if (write_compsp.EQ.1.OR.write_compsp.EQ.3) CLOSE(10)
  if (write_compsp.EQ.2.OR.write_compsp.EQ.3) CLOSE(20)

end subroutine compsp
