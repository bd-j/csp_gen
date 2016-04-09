! ------------------------------------
! Routines used to evaluate the indefinite integrals that arise when exactly
! integrating interpolated spectra weighted by particular SFHs.
!
! Adding new SFHs is as easy as adding new sfh%type cases in both sfhint_*
! functions below (and making sure the sfhlimits.f90 make sense)
! ------------------------------------

function sfhint(sspind, logt, sfh):
  ! Wrapper on the sfhint_* routines to choose the correct interpolation type
  
  use sps_vars, only: interpolation_type
  implicit none
  integer, intent(in) :: sspind
  real(SP), intent(in) :: logt
  type(SFHPARAMS), intent(in) :: sfh
  real(SP), intent(out) :: sfhint
  
  if (interpolation_type.eq.0) then
     sfhint = sfhint_log(sspind, logt, sfh)
  else
     sfhint = sfhint_lin(sspind, logt, sfh)
  endif
  
end function sfhint
   

function sfhint_log(sspind, logt, sfh)
  ! Indefinite integral of the interpolation weight in log time, weighted by
  ! the SFH, and evaluated at `logt`.  In detail, this function returns:
  !   (\int dt \, \mathrm{SFR}(t) \, (x - \log t) )|_{\texttt{logt}}
  ! where x = time_full(ssp_ind).
  !
  ! Inputs
  ! -------
  !
  ! sspind:
  !    The index of the SSP forming the bracket with the SSP you are getting a
  !    weight for.
  !
  ! logt:
  !    Where the integral is evaluated.
  !
  ! sfh:
  !    Structure containing the SFH parameters in units of years, including an
  !    integer specifiying the form of SFR(t).
  !
  ! Outputs
  !--------
  !  sfhint:
  !    The exact indefinite integral, evaluated at `logt`
  
  use sps_vars, only: time_full, pset
  use sps_utils, only: expi
  implicit none
  integer, intent(in) :: sspind
  real(SP), intent(in) :: logt
  type(SFHPARAMS), intent(in) :: sfh
  
  real(SP) :: loge
  real(SP) :: logage, tprime ! intermediate time variables
  real(SP) :: a, b, c !dummy variables used to break up long expressions
  
  real(SP), intent(out) :: sfhint_log

  logage = time_full(sspind)
  
  ! Convert from Gyr to yrs.
  ! This is stupid and time consuming to have in the inner loop. `sfh` should
  ! be a structure that contains the SFH variables (including type) already in
  ! the proper units
  ! sfh%tage = pset%tage * 1e9
  ! sfh%tau = pset%tau * 1e9
  ! sfh%sf_trunc = pset%sf_trunc * 1e9
  ! sfh%sf_slope = pset%sf_slope / 1e9
  loge = log10(exp(1))

  if (sfh%type.eq.0) then
     ! SFR = Constant ~ 1
     sfhint_log = 10**logt * (logage - logt + loge)
     
  else if (sfh%type.eq.1) then
     ! SFR = exponential ~ exp(-T/tau)
     tprime = 10**logt / sfh%tau
     sfhint_log = (logage - logt) * exp(tprime) + loge) * expi(tprime)
     
  else if (sfh%type.eq.4) then
     ! SFR = delayed exponential ~ T/tau exp(-T/tau)
     tprime = 10**logt / sfh%tau ! t/tau
     a = 10**logt * (logt - logage) + &
          sfh%tage * logage + &
          sfh%tau * (logage - loge)
     b = sfh%tau * (logt * exp(tprime) - loge * expi(tprime))
     sfhint_log = a * exp(tprime) - b * (tage / sfh%tau + 1)
     
  else if (sfh%type.eq.5) then
     !SFR = linear ~ (1 - sf_slope * (T - T_trunc)), T > T_trunc
     tprime = max(0, sfh%tage - sfh%sf_trunc) !t_q
     a = 1 - sfh%sf_slope * tprime
     b = a * 10**logt * (logage - logt + loge)
     c = sfh%sf_slope * (10**logt)**2 / 2 * (logage - logt + loge / 2)
     sfhint_log = b + c
     
  endif
  
end function sfhint_log


function sfhint_lin(sspind, t, sfh)
  ! Indefinite integral of the interpolation weight in linear time, weighted by
  ! the SFH, and evaluated at `logt`.  In detail, this function returns:
  !   (\int dt \, \mathrm{SFR}(t) \, (x - t) )|_{\texttt{t}}
  ! where x = 10**time_full(ssp_ind).
  !
  ! Inputs
  ! -------
  !
  ! sspind:
  !    The index of the SSP forming the bracket with the SSP you are getting a
  !    weight for.
  !
  ! logt:
  !    Where the integral is evaluated.
  !
  ! sfh:
  !    Structure containing the SFH parameters in units of years, including an
  !    integer specifiying the form of SFR(t).
  !
  ! Outputs
  !--------
  !  sfhint:
  !    The indefinite integral, evaluated at `logt`

  use sps_vars, only: time_full
  implicit none
  integer, intent(in) :: sspind 
  real(SP), intent(in) :: t
  type(SFHPARAMS), intent(in) :: sfh
  
  real(SP) :: age, tprime
  real(SP) :: loge, a

  real(SP), intent(out) :: sfhint_lin

  ! Convert from Gyr to yrs.
  ! This is stupid and time consuming to have in the inner loop. `sfh` should
  ! be a structure that contains the SFH variables (including type) already in
  ! the proper units
  ! sfh%tage = pset%tage * 1e9
  ! sfh%tau = pset%tau * 1e9
  ! sfh%sf_trunc = pset%sf_trunc * 1e9
  ! sfh%sf_slope = pset%sf_slope / 1e9
  loge = log10(exp(1.0))

  !convert from log(age_ssp) to age_ssp
  age = 10**time_full(sspind)

  if (sfh%type.eq.0) then
     ! SFR = Constant ~ 1
     sfhint_lin = age * t - t**2 / 2

  else if (sfh%type.eq.1) then
     ! SFR = exponential ~ 1/tau * exp(-T/tau)
     tprime = t / sfh%tau
     sfhint_lin = (age - t + sfh%tau) * exp(tprime)

  else if (sfh%type.eq.4) then
     ! SFR = delayed exponential ~ T/tau**2 * exp(-T/tau)
     tprime = t / sfh%tau
     a = sfh%tage * age - (sfh%tage + age) * (t - sfh%tau) + &
          t**2 - 2*t*sfh%tau + 2*sfh%tau**2
     sfhint_lin = a * exp(tprime)
     
  else if (sfh%type.eq.5) then
     ! SFR = linear ~ (1 - sf_slope * (T - T_trunc)), T > T_trunc
     tprime = max(0, sfh%tage - sfh%sf_trunc) !t_q
     a = 1 - sfh%sf_slope * tprime
     sfhint_lin = a * age * t + (sfh%sf_slope*age - a) * t**2 / 2 - sfh%sf_slope * t**3 / 3

end function sfhint_lin
