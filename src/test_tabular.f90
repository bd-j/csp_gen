PROGRAM TEST_TABULAR
  !set up modules
  USE sps_vars; USE sps_utils
  
  IMPLICIT NONE

  !NB: the various structure types are defined in sps_vars.f90
  !    variables not explicitly defined here are defined in sps_vars.f90
  INTEGER :: i
  !define variable for SSP spectrum
  REAL(SP), DIMENSION(ntfull,nspec)  :: spec_ssp
  !define variables for Mass and Lbol info
  REAL(SP), DIMENSION(ntfull)    :: mass_ssp,lbol_ssp
  CHARACTER(100) :: file2=''
  !structure containing all necessary parameters
  TYPE(PARAMS) :: pset
  !define structure for CSP spectrum
  TYPE(COMPSPOUT), DIMENSION(ntfull) :: ocompsp

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!
  
  ! Now we're going to show you how to use full  
  ! metallicity-dependent info                   
  
  imf_type = 0    ! Salpeter IMF
  pset%sfh = 0    ! compute SSP
  pset%zmet = 20 ! solar with padova/basel
  pset%mag_compute = 1
  !here we have to read in all the librarries
  CALL SPS_SETUP(pset%zmet)

  !compute all SSPs (i.e. at all Zs)
  !nz and the various *ssp_zz arrays are stored 
  !in the common block set up in sps_vars.f90
  !DO i=1,nz
  !   pset%zmet = i
  !   CALL SSP_GEN(pset,mass_ssp_zz(:,i),&
  !        lbol_ssp_zz(:,i),spec_ssp_zz(:,:,i))
  !ENDDO

  !compute the SSP
  CALL SSP_GEN(pset,mass_ssp,lbol_ssp,spec_ssp)

  !run compsp for a tabulated sfh with a single z (given by zmet)
  !NB: one must have setup all the SSPs, as was done in the DO-loop above
  pset%sfh = 2
  !pset%sfh_filename = 'test_sfh.dat'
  file2    = 'CSP_tabsfh_test.out'
  CALL COMPSP(1,1,file2,mass_ssp,lbol_ssp,&
          spec_ssp,pset,ocompsp)


END PROGRAM TEST_TABULAR
