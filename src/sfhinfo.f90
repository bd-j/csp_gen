subroutine sfhinfo(pset, age, mfrac, sfr)
  ! Get the SFR integrated from T=0 to T=age, normalized by the SFR integrated
  ! from T=0 to T=Tmax, where Tmax is the maximum isochrone/SSP age.
  implicit none
  type(PARAMS), intent(in) :: pset
  real(SP), intent(in), OPTIONAL :: tage
  
  real(SP) :: age, Tmax, Tprime, Tz
  real(SP) :: mass_tau, mass_linear
  real(SP) :: power, m, sfr_q
  
  if (pset%sfh.eq.1) power = 1.
  if (pset%sfh.eq.4).or.(pset%sfh.eq.5) power = 2.

  ! Integration limits are from 0 to Tmax and 0 to Tprime, adjusted for sf_start
  ! The maximum age in the isochrone
  Tmax = 10**(max(time_full) - 9) - pset%sf_start
  ! The age to which we are integrating
  Tprime = age - pset%sf_start  
  ! Deal with truncation that happens after sf_start but before Tmax and/or Tprime
  if (pset%sf_trunc.gt.pset%sf_start) then
     Tmax = min(Tmax, pset%sf_trunc - pset%sf_start)
     Tprime = min(Tprime, pset%sf_trunc - pset%sf_start)
  endif
  
  ! Fraction of the tau-model mass formed by Tprime
  total_mass_tau = pset%tau * gammainc(power, Tmax/pset%tau)
  mass_tau = pset%tau * gammainc(power, Tprime/pset%tau)
  sfr_tau = (Tprime/pset%tau)**(power-1.) * exp(Tprime/pset%tau)

  ! Add the constant and burst or linear portions
  if ((pset%sfh.eq.1).or.(pset%sfh.eq.4)) then
     ! Fraction of the constant mass formed by Tprime
     total_mass_constant = Tmax
     mass_constant = Tprime
     sfr_const = 1.0

     !Fraction of the burst mass formed by Tprime
     if (Tprime.gt.(pset%tburst-pset%sf_start)) then
        mass_burst = 1.0
     else
        mass_burst = 0.0
     endif

     mfrac = (1. - pset%const - pset%fburst) * mass_tau / total_mass_tau + &
          pset%const * mass_constant / total_mass_constant + &
          pset%fburst * mass_burst
     ! N.B. for Tprime = tburst, sfr is infinite, but we ignore that case
     sfr = (1. - pset%const - pset%fburst) * sfr_tau / total_mass_tau + &
          pset%const * sfr_const / total_mass_const

  else if (pset%sfh.eq.5) then
   
  ! Add the linear portion. for simha.
  ! This is integral of (1 - m * (T - Tmax)) from Ttrunc to Tzero

     sfr_q = (Tmax/pset%tau) * exp(-Tmax/pset%tau)
     if (m.gt.0) then
        Tz = Tmax + 1.0 / m
     else
        Tz = age
     endif
     if ((Tz.lt.Tmax).or.(Tz.gt.age)) then
        Tz = age
     endif
     mass_linear = sfr_q * ((Tz - Tmax) - m/2.*(Tz**2 + Tmax**2) + m*Tz*Tmax)
  endif
  
 
end subroutine sfhinfo

function gammainc(power, arg)
  !
  ! Calculate incomplete gamma for a=1 or 2
  if (power.eq.2) then
     gammainc = 1.0 - exp(-arg) - arg * exp(-arg)
  else if (power.eq.1)
     gammainc = 1.0 - exp(-arg)
  endif
  
end function gammainc
