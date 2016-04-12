subroutine sfhinfo(pset, age, mfrac, sfr)

  mfrac = sfint(pset, age) / sfint(pset, pset%tage)
  sfr = get_sfr(pset)
  
end subroutine sfhinfo


function sfint(pset, age)
  ! ugh.
  !
  implicit none
  type(PARAMS), intent(in) :: pset

  real(SP) :: Tmax, Tz
  real(SP) :: mass_tau, mass_linear
  real(SP) :: power, m, sfr_q
  
  if (pset%sfh.eq.1) power = 1.
  if (pset%sfh.eq.4).or.(pset%sfh.eq.5) power = 2.

  ! mass in the tau portion
  if ((pset%sf_trunc.le.0).or.(pset%sf_trunc.gt.age)) then
     Tmax = age
  else
     Tmax = pset%sf_trunc
  endif
  mass_tau = pset%tau * gammainc(power, Tmax/pset%tau)* (1 - pset%const)

  ! sfr at sf_trunc
  sfr_q = (Tmax/pset%tau) * exp(-Tmax/pset%tau) * (1 - pset%const) + pset%const
  
  ! Mass in the linear portion.
  ! This is integral of (1 - m * (T - Tmax)) from Tmax to Tzero
  ! we then multiply it by sfr_q
  if (m.gt.0) then
     Tz = Tmax + 1.0/m
  else
     Tz = age
  endif
  if ((Tz.lt.Tmax).or.(Tz.gt.age)) then
     Tz = age
  endif
  mass_linear = sfr_q * ((Tz - Tmax) - m/2.*(Tz**2 + Tmax**2) + m*Tz*Tmax)
  
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

