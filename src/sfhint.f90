! ------------------------------------
! Routines used to evaluate the indefinite integrals that arise when exactly
! integrating interpolated spectra weighted by particular SFHs.
! ------------------------------------


function sfhint_log(sspind, logt)
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
  ! Outputs
  !--------
  !  sfhint:
  !    The indefinite integral, evaluated at `logt`
  
  use sps_vars, only: ntfull, time_full, pset
  use sps_utils, only: expi
  implicit none
  integer, intent(in) :: sspind 
  real(SP), intent(in) :: logt
  real(SP) :: tage, tau, sft, sfslope, loge
  real(SP) :: logage, tprime
  real(SP), intent(out), dimension(ntfull) :: sfhint

  logage = time_full(sspind)
  
  ! Convert from Gyr to yrs.  This is stupid and time consuming to have in the inner loop
  tage = pset%tage * 1e9
  tau = pset%tau * 1e9
  sft = pset%sf_trunc * 1e9
  sf_slope = pset%sf_slope / 1e9
  loge = log10(exp(1))

  ! SFR = Constant ~ 1
  sfhint_log = 10**logt * (logage - logt + loge)

  ! SFR = exponential ~ exp(-T/tau)
  tprime = 10**logt / tau
  sfhint_log = (logage - logt) * exp(tprime) + loge) * expi(p1)

  ! SFR = delayed exponential ~ T/tau exp(-T/tau)
  tprime = 10**logt / tau ! t/tau
  br = 10**logt * (logt - logage) + &
       tage * logage + &
       tau * (logage - loge)
  h = tau * (logt * exp(tprime) - loge * expi(tprime))
  sfhint_log = br * exp(tprime) - (tage / tau + 1) * h

  !SFR = linear ~ (1 - sf_slope * (T - T_trunc)), T > T_trunc
  tprime = max(0, tage-sf_trunc) !t_q
  k = 1 - sf_slope * tprime
  term1 = k * 10**logt * (logage - logt + loge)
  term2 = sf_slope * (10**logt)**2 / 2 * (logage - logt + loge / 2)
  sfhint_log = term1 + term2

end function sfhint_log


function sfhint_lin(sspind, t)
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
  ! Outputs
  !--------
  !  sfhint:
  !    The indefinite integral, evaluated at `logt`

  use sps_vars, only: ntfull, time_full, pset
  use sps_utils, only: expi
  implicit none
  integer, intent(in) :: sspind 
  real(SP), intent(in) :: t
  real(SP) :: tage, tau, sft, sfslope, loge
  real(SP) :: tprime
  real(SP) :: age
  real(SP), intent(out), dimension(ntfull) :: sfhint

  ! Convert from Gyr to yrs
  ! Bad to have this in the inner loop
  tage = pset%tage * 1e9
  tau = pset%tau * 1e9
  sft = pset%sf_trunc * 1e9
  sf_slope = pset%sf_slope / 1e9
  loge = log10(exp(1))

  !convert from log(age_ssp) to age_ssp
  age = 10**time_full(sspind)
  
  ! SFR = Constant ~ 1
  sfhint_lin = age * t - t**2 / 2

  ! SFR = exponential ~ exp(-T/tau)
  tprime = t / tau
  sfhint_lin = (age - t + tau) * exp(tprime)

  ! SFR = delayed exponential ~ T/tau exp(-T/tau)
  tprime = t / tau
  bracket = tage * age - (tage + age) * (t - tau) + t**2 - 2*t*tau + 2*tau**2
  sfhint_lin bracket * exp(tprime)

  !SFR = linear ~ (1 - sf_slope * (T - T_trunc)), T > T_trunc
  tprime = max(0, tage-sf_trunc) !t_q
  k = 1 - sf_slope * tprime
  sfhint_lin = k * age * t + (sf_slope*age - k) * t**2 / 2 - sf_slope * t**3 / 3

end function sfhint_lin


