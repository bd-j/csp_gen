subroutine sfhinfo(pset, age, mfrac, sfr)

  mfrac = sfint(pset, age) / sfint(pset)
  sfr = get_sfr(pset, age)
  
end subroutine sfhinfo


function sfint(pset, tage)
  ! ugh.  This is a disaster
  ! Get the SFR integrated from T=0 to T=age
  implicit none
  type(PARAMS), intent(in) :: pset
  real(SP), intent(in), OPTIONAL :: tage
  
  real(SP) :: age, Tmax, Tz
  real(SP) :: mass_tau, mass_linear
  real(SP) :: power, m, sfr_q
  
  if (pset%sfh.eq.1) power = 1.
  if (pset%sfh.eq.4).or.(pset%sfh.eq.5) power = 2.
  
  ! Check for specified age
  if present(tage) then
     age = tage
  else
     age = pset%tage
  endif
  
  ! Get the mass in the tau portion.
  ! This is the integral of (T/tau)^power e^(-T/tau) from 0 to Tmax
  if ((pset%sf_trunc.le.0).or.(pset%sf_trunc.gt.age)) then
     Tmax = age - pset%sfstart
  else
     Tmax = pset%sf_trunc - pset%sfstart
  endif
  mass_tau = pset%tau * gammainc(power, Tmax/pset%tau)

  ! Add constant and burst
  if ((pset%sfh.eq.1).or.(pset%sfh.eq.4)) then
     
     sfint = (1. - pset%const - pset%fburst) * mass_tau + &
              pset%const*min/(pset%sf_trunc)
  endif
  
  ! Add the linear portion. for simha.
  ! This is integral of (1 - m * (T - Tmax)) from Tmax to Tzero

  if (pset%sfh.eq.5) then
     sfr_q = (Tmax/pset%tau) * exp(-Tmax/pset%tau)
     if (m.gt.0) then
        Tz = Tmax + 1.0/m
     else
        Tz = age
     endif
     if ((Tz.lt.Tmax).or.(Tz.gt.age)) then
        Tz = age
     endif
     mass_linear = sfr_q * ((Tz - Tmax) - m/2.*(Tz**2 + Tmax**2) + m*Tz*Tmax)
  else
     mass_linear = 0
  endif

  mass_par = mass_tau + mass_linear
  
  ! integral of 1 from 0 to Tmax
  
  ! add constant and burst
  mass = mass_tau + mass_linear
  if ((age.gt.pset%tburst).and.(pset%fburst.gt.tiny_number)) then
     use_burst = 1.
  else
     use_burst = 0.
  endif
  
  sfint = (1. - pset%const - pset%fburst) * (mass_tau + mass_linear) + &
       pset%const*min(age, sf_trunc)/(pset%sf_trunc) + &
       use_burst * pset%fburst

  
end function sfint


