! this is what needs to be added to sps_vars
MODULE VARS
  implicit none
  save

  integer :: interpolation_type = 0
  real(SP) :: tiny_logt = -3
  
  TYPE SFHPARAMS
     REAL(SP) :: tau=1.0, tage=0., tburst=0., sf_trunc=0., sf_slope=0., tq=0., t0=0., tb=0.
     INTEGER :: type=0, use_simha_limits=0
  END TYPE SFHPARAMS

end module vars
