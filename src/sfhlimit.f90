function sfhlimit(tlim, sfh)
  ! Deal with all the logic for computing the limits of integration in the
  ! weight calculations.  This should return the input `tlim`, clipped to valid
  ! limits, such that if no SFR occurs between t1 and t2 then sfhlimit(t1, sfh)
  ! = sfhlimit(t2, sfh)
  !
  ! Inputs
  ! ----------
  !
  ! tlim:
  !    The proposed limit of the integration, in log(lookback_time)
  !
  ! sfh:
  !    An sfhparams structure containing the relevant special lookback times
  !
  use sps_vars, only: tiny_logt
  implicit none
  
  real(SP), intent(in) :: tlim
  type(SFHPARAMS), intent(in) :: sfh

  real(SP), intent(out) :: sfhlimit
  
  real(SP) :: tlo, thi
  
  ! For the simha linear portion, we integrate from sf_trunc to tage or the
  ! zero crossing, whichever is smaller but still greater than sf_trunc.
  ! For everything else we integrate from 0 to tage or sf_trunc, whichever is
  ! smaller but non-zero.
  ! Of course, we are doing this using the lookback time versions!!!
  if (sfh%use_simha_limits.eq.1) then
     tlo = sfh%t0
     thi = sfh%tq
  else
     tlo = sfh%tq
     thi = sfh%tage
  endif
  
  ! convert to log, taking care of possible zeros
  tlo = log10(max(tlo, 10**tiny_logt))
  thi = log10(max(thi, 10**tiny_logt))
  
  ! Finally, clip proposed limit to the upper and lower value
  sfhlimit = clip(tlim, tlo, thi)
  
end function sfhlimit

function clip(x, lo, hi):
  ! stupid function to clip x to the range (lo, hi)
  !
  implicit none
  
  real(SP), intent(in) :: x, lo, hi
  real(SP), intent(out) :: clip
  
  clip  = MIN(MAX(x, lo), hi)
  
end function clip

