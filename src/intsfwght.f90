! ------------------------------------
! Routines used to evaluate the indefinite integrals that arise when exactly
! integrating interpolated spectra weighted by particular SFHs.
!
! Adding new SFHs is as easy as adding new sfh%type cases in both sfhint_*
! functions below (and making sure the sfhlimits.f90 make sense)
! ------------------------------------

function intsfwght(sspind, logt, sfh):
  ! Wrapper on the sfhint_* routines to choose the correct interpolation type
  ! and calculate the definite integral between the given limits
  
  use sps_vars, only: interpolation_type
  implicit none
  
  integer, intent(in) :: sspind
  real(SP), intent(in), dimension(2) :: logt
  type(SFHPARAMS), intent(in) :: sfh

  real(SP), intent(out) :: intsfwght

  !real(SP) :: dt

  if (interpolation_type.eq.0) then
     !dt = logt(2) - logt(1)
     intsfwght = (sfhint_log(sspind, logt(2), sfh) - sfhint_log(sspind, logt(1), sfh))
  else
     !dt = 10**logt(2) - 10**logt(1)
     intsfwght = (sfhint_lin(sspind, 10**logt(2), sfh) - sfhint_lin(sspind, 10**logt(1), sfh))
  endif
  !intsfwght = intsfwght / dt
  
end function sfhint
   

function sfhint_log(sspind, logt, sfh)
  ! Evaluates the indefinite integral of the interpolation weight in log time,
  ! weighted by the SFH, and evaluated at `logt`.  In detail, this function
  ! returns:
  !   (\int dt \, \mathrm{SFR}(t) \, (x - \log t) )|_{\texttt{logt}}
  ! where x = time_full(ssp_ind) and `t` is *lookback* time.
  !
  ! Inputs
  ! -------
  !
  ! sspind:
  !    The index of the SSP forming the bracket with the SSP you are getting a
  !    weight for.
  !
  ! logt:
  !    Where the indefinite integral is evaluated.
  !
  ! sfh:
  !    Structure containing the SFH parameters in units of years, including an
  !    integer `type` specifiying the form of SFR(t).
  !
  ! Outputs
  !--------
  !  sfhint:
  !    The exact indefinite integral, evaluated at `logt`
  
  use sps_vars, only: time_full, tiny_logt
  use sps_utils, only: expi
  implicit none
  
  integer, intent(in) :: sspind
  real(SP), intent(in) :: logt
  type(SFHPARAMS), intent(in) :: sfh

  real(SP), intent(out) :: sfhint_log

  real(SP) :: loge
  real(SP) :: logage, tprime ! intermediate time variables
  real(SP) :: a, b, c !dummy variables used to break up long expressions

  loge = log10(exp(1))
  ! zero index means use ~0 age
  if (sspind.gt.0) then
     logage = time_full(sspind)
  else
     logage = tiny_logt
  endif
     
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
  ! Evaluates the indefinite integral of the interpolation weight in linear
  ! time, weighted by the SFH, and evaluated at `logt`.  In detail, this
  ! function returns:
  !   (\int dt \, \mathrm{SFR}(t) \, (x - t) )|_{\texttt{t}}
  ! where x = 10**time_full(ssp_ind), and `t` is *lookback* time.
  !
  ! Inputs
  ! -------
  !
  ! sspind:
  !    The index of the SSP forming the bracket with the SSP you are getting a
  !    weight for.
  !
  ! t:
  !    Where the indefinite integral is evaluated.
  !
  ! sfh:
  !    Structure containing the SFH parameters in units of years, including an
  !    integer `type` specifiying the form of SFR(t).
  !
  ! Outputs
  !--------
  !  sfhint:
  !    The indefinite integral, evaluated at `t`

  use sps_vars, only: time_full, tiny_logt
  implicit none
  
  integer, intent(in) :: sspind 
  real(SP), intent(in) :: t
  type(SFHPARAMS), intent(in) :: sfh

  real(SP), intent(out) :: sfhint_lin
  
  real(SP) :: age, tprime
  real(SP) :: loge, a

  loge = log10(exp(1.0))

  ! convert from log(age_ssp) to age_ssp,
  ! accounting for the case sspind=0 (where age~0)
  if (sspind.gt.0) then
     age = 10**time_full(sspind)
  else
     age = 10**tiny_logt
  endif

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
