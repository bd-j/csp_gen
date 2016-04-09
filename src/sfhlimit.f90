! Deal with all the logic for computing the limits of integration in the weight
! calculations.
!
! This should return the input, clipped to valid limits, such that if no SFR
! occurs between t1 and t2 then sfhlimit(t1) = sfhlimit(t2)


function sfhlimit(tlim, sfh)
  !
  !
  !
  use sps_vars, only: interpolation_type, tiny_logt
  implicit none
  
  real(SP), intent(in) :: tlim
  type(SFHPARAMS), intent(in) :: sfh
  
  real(SP) :: tq, t0, tlo, thi

  real(SP), intent(out) :: sfhlimit
  
  ! convert sf_trunc to to t_lookback
  if ((sfh%sf_trunc.le.0).or.(sfh%sf_trunc.gt.sfh%tage)) then
     tq = 0
  else
     tq = sfh%tage - sfh%sf_trunc
  endif

  ! For the simha linear portion, we integrate from sf_trunc to tage or the
  ! zero crossing, whichever is smaller but still greater than sf_trunc.
  ! For everything else we integrate from 0 to tage or sf_trunc, whichever is
  ! smaller but non-zero.
  ! Of course, we are converting everything to lookback times as well!
  if (sfh%simha_limits.eq.1) then
     ! get zero crossing time (in lookback time), avoiding divison by zero
     if (abs(sfh%sf_slope).gt.tiny_number) then
        t0 = tq - 1. / sfh%sf_slope
     else
        t0 = 0
     endif
     ! if zero crossing time outside the linear regime, set it to zero
     if ((t0.gt.tq).or.(t0.le.0)) then
        t0 = 0
     endif
     tlo = t0
     thi = tq
  else
     tlo = tq
     thi = sfh%tage
  endif
  
  ! convert to log, taking care of possible zeros
  if (interpolation_type.eq.0) then
     tlo = log10(max(tlo, 10**tiny_logt))
     thi = log10(max(thi, 10**tiny_logt))
  endif

  ! Finally, clip limit to the upper and lower value
  sfhlimit = clip(tlim, tlo, thi)
  
end function sfhlimit

function clip(x, lo, hi):
  ! stupid function to clip x to the range (lo, hi)
  !
  real(SP), intent(in) :: x, lo, hi
  real(SP), intent(out) :: clip
  implicit none
  
  clip  = MIN(MAX(x, lo), hi)
  
end function clip

