! Deal with all the logic for computing the limits of integration in the weight calculations.
! This should return the input, clipped to valid limits, such that if no SFR occurs between t1 and t2 then sfhlimit(t1) = sfhlimit(t2)


function sfhlimit(tlim)
  !
  use sps_vars, only: pset
  implicit none
  
  mint_logt = -3
  
  if ((sf_trunc.le.0).or.(sf_trunc.gt.tage)) then
     tq = 0
  else
     ! convert to t_lookback
     tq = tage - sf_trunc
  endif

  ! convert to log, taking care of possible zeros
  if (interp_type.eq.'logarithmic') then
     tq = max(tq, 10**mint_log)
     tage = max(tage, 10**mint_log)
     tq = log10(tq)
     tage = log10(tage)
  endif
  ! need to convert clip to ridiculous min(max(min logic, looped over values
  sfhlimit = clip(tlim, tq, tage)
  
end function sfhlimit



function sfhlimit_simha(tlim)
  
  use sps_vars, only: pset, tiny_number
  implicit none
  
  mint_logt = -3

  ! get truncation time in units of lookback time
  if ((sf_trunc.le.0).or.(sf_trunc.gt.tage)) then
     tq = 0
  else
     tq = tage - sf_trunc
  endif
  ! get zero crossing time, avoiding divison by zero
  if (abs(sf_slope).gt.tiny_number) then
     t0 = tq - 1. / sf_slope
  else
     t0 = 0
  endif
  ! if zero crossing time outside the linear regime, set it to zero
  if ((t0.gt.tq).or.(t0.le.0)) then
     t0 = 0
  endif

  ! convert to log, taking care of possible zeros
  if (interp_type.eq.'logarithmic') then
     t0 = max(t0, 10**mint_log)
     tq = max(tq, 10**mint_log)
     t0 = log10(t0)
     tq = log10(tq)
  endif
  
  ! need to convert np.clip to ridiculous min(max(min logic
  sfhlimit = clip(tlim, t0, tq)
  
