  INTERFACE
     SUBROUTINE SFHINFO(pset, age, mfrac, sfr, frac_linear)
       USE sps_vars
       TYPE(PARAMS), INTENT(in) :: pset
       REAL(SP), INTENT(in) :: age
       REAL(SP), INTENT(out) :: mfrac, sfr, frac_linear
     END SUBROUTINE SFHINFO
  END INTERFACE

  INTERFACE
     SUBROUTINE CSP_GEN(mass_ssp, lbol_ssp, spec_ssp, mdust_ssp, pset, tage,&
                        mass_csp, lbol_csp, spec_csp, mdust_csp)
       USE sps_vars
       REAL(SP), DIMENSION(ntfull), INTENT(in) :: mass_ssp, lbol_ssp, mdust_ssp
       REAL(SP), DIMENSION(nspec, ntfull), INTENT(in) :: spec_ssp
       TYPE(PARAMS), intent(in) :: pset
       REAL(SP), INTENT(in)  :: tage
       REAL(SP), INTENT(out) :: mass_csp, lbol_csp, mdust_csp
       REAL(SP), INTENT(out), DIMENSION(nspec) :: spec_csp
     END SUBROUTINE CSP_GEN
  END INTERFACE

  INTERFACE
     FUNCTION SFH_WEIGHT(sfh, imin, imax)
       USE sps_vars
       TYPE(SFHPARAMS), INTENT(in) :: sfh
       INTEGER, INTENT(in) :: imin, imax
       REAL(SP), DIMENSION(ntfull) :: sfh_weight
     END FUNCTION SFH_WEIGHT
  END INTERFACE

  INTERFACE
     FUNCTION SFHLIMIT(tlim, sfh)
       USE sps_vars
       TYPE(SFHPARAMS), INTENT(in) :: sfh
       REAL(SP), INTENT(in) :: tlim
       REAL(SP) :: sfhlimit
     END FUNCTION SFHLIMIT
  END INTERFACE

  INTERFACE
     FUNCTION INTSFWGHT(sspind, logt, sfh)
       USE sps_vars
       TYPE(SFHPARAMS), INTENT(in) :: sfh
       INTEGER, INTENT(in) :: sspind
       REAL(SP), DIMENSION(2), INTENT(in) :: logt
       REAL(SP) :: intsfwght
     END FUNCTION INTSFWGHT
  END INTERFACE
