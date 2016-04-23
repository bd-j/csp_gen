function sfh_weight(sfh, imin, imax)
  ! Function to calculate the weights
  !
  !
  use sps_vars, only: ntfull, time_full, pset, interpolation_type, tiny_logt
  use sps_utils, only: sfhint, sfhlimits
  implicit none

  type(SFHPARAMS), intent(in) :: sfh
  integer, intent(in) :: imin, imax

  real(SP), intent(out), dimension(ntfull) ::sfh_weight=0.

  integer :: i, do_zero_bin=0
  real(SP) dimension(2) :: tlim
  real(SP) :: dt 
  real(SP), dimension(ntfull) :: left=0., right=0.

  ! Check if this is an SSP.  If so, do simple weights and return
  if (sfh%type.eq.-1) then
     sfh_weight = ssp_weight(sfh%tb)
     return
  endif

  ! Check if we need to do the zero bin, using imin=0 as a flag.
  ! If so, the zero bin is added at the end of this function.
  if (imin.eq.0) then
     do_zero_bin = 1
     imin = 1
  endif

  ! Loop over each SSP and calculate its weight in the given sfh
  do i=imin,imax
     if (i.gt.1) then
        ! There is a younger (`left`) bin, and we calculate its contribution to
        ! the weight.
        ! First calculate actual limits for the younger bin.  
        tlim(1) = sfhlimit(time_full(i-1), sfh)
        tlim(2) = sfhlimit(time_full(i), sfh)
        ! The elements of `tlim` will be equal if there is no valid SFR in the
        ! younger bin; only proceed if there is a non-zero sfr in the younger bin
        if (tlim(1).ne.tlim(2)) then
           dt = delta_time(time_full(i-1), time_full(i))
           ! Note sign flip here
           left(i) = 0. - intsfwght(i-1, tlim, sfh) / dt
        endif
     endif
     if (i.lt.ntfull) then
        ! There is an older (`right`) bin, we calculate its contribution to the weight
        tlim(1) = sfhlimit(time_full(i), sfh)
        tlim(2) = sfhlimit(time_full(i+1), sfh)
        ! The elements of `tlim` will be equal if there is no valid SFR in the
        ! younger bin; only proceed if there is a non-zero sfr in the older bin
        if (tlim(1).ne.tlim(2)) then
           dt = delta_time(time_full(i), time_full(i+1))
           right(i) = intsfwght(i+1, tlim, sfh) / dt
        endif
     endif
  enddo

  sfh_weight = right + left

  ! Do we need to add weights from the zero bin to the first SSP?
  !   We assume the t~0 spectrum is the same as the t=10**time_full(1) spectrum
  !   (i.e., nearest neighbor extrapolation), so the t~0 weight gets added to
  !   sfh_weight(1)
  if (do_zero_bin.eq.1) then
     tlim(1) = sfhlimit(tiny_logt, sfh)
     tlim(2) = sfhlimit(time_full(1), sfh)
     if (tlim(1).ne.tlim(2)) then
        dt = delta_time(tiny_logt, time_full(1))
        ! contribution of i=1 to younger bin
        sfh_weight(1) = sfh_weight(1) - intsfwght(0, tlim, sfh) / dt
        ! contribution of i=0 to older bin
        sfh_weight(1) = sfh_weight(1) + intsfwght(1, tlim, sfh) / dt
     endif
  endif

end function sfh_weight


function ssp_weight(tb)
  ! Quick function to calculate SSP weights for a single arbitrary age.
  !
  ! Inputs
  ! --------
  !
  ! tb:
  !   burst time, linear years of lookback time
  
  use sps_vars, only: time_full, ntfull
  implicit none

  real(SP), intent(in) :: tb

  real(SP), intent(out), dimension(ntfull) :: ssp_weight = 0.

  integer :: imin
  real(SP) :: log_tb, dt

  log_tb = log10(tb)
  imin = min(max(locate(time_full, log_tb), 1, ntfull-1))
  dt = delta_time(time_full(imin), time_full(imin+1))
  ssp_weight(imin) = delta_time(log_tb, timefull(imin+1)) / dt
  ssp_weight(imin+1) = delta_time(timefull(imin), log_tb) / dt

end function ssp_weight

  
function delta_time(logt1, logt2)
  ! Dumb function to properly calculate dt based on interpolation type.
  !
  ! Returns (logt2 - logt1), or (10**logt2 - 10**logt1)

  use sps_vars, only: interpolation_type
  implicit none

  real(SP), intent(in) :: logt1, logt2
  real(SP), intent(out) :: delta_time
  if (interpolation_type.eq.1) then
     delta_time = logt2 - logt1
  else
     delta_time = 10**logt2 - 10**logt1
  endif

end function delta_time
