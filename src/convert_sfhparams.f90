subroutine convert_sfhparams(pset, sfh)
  ! Convert the pset values to yrs and pre-calculate some useful lookback
  ! times, and store in a structure.  Note that this subroutine uses but does
  ! not alter the pset values.
  !
  ! Outputs
  ! ----------
  ! sfh:
  !    An SFHPARAMS structure, containing SFH parameters in yrs, and some
  !    useful lookback times.
  !       - `tq` is the truncation time, in lookback time
  !       - `t0` is the zero crossing time for a linear SFH, in lookback time.
  !
  use sps_vars, only: pset
  implicit none
  
  type(SFHPARAMS), intent(in) :: pset
  
  type(SFHPARAMS), intent(inout) :: sfh

  ! Convert units from Gyr to yr
  sfh%tage = pset%tage * 1e9
  sfh%tau = pset%tau * 1e9
  sfh%sf_trunc = pset%sf_trunc * 1e9
  sfh%sf_slope = pset%sf_slope / 1e9

  ! convert sf_trunc to to t_lookback
  if ((sfh%sf_trunc.le.0).or.(sfh%sf_trunc.gt.sfh%tage)) then
     sfh%tq = 0
  else
     sfh%tq = sfh%tage - sfh%sf_trunc
  endif
  ! get zero crossing time (in lookback time), avoiding divison by zero
  if (abs(sfh%sf_slope).gt.tiny_number) then
     sfh%t0 = sfh%tq - 1. / sfh%sf_slope
  else
     sfh%t0 = 0
  endif
  ! If the zero crossing time is outside the linear regime, set it to zero.
  if ((sfh%t0.gt.sfh%tq).or.(sfh%t0.le.0)) then
     sfh%t0 = 0
  endif

end subroutine convert_sfhparams
