function sfh_weight(sfh_type)

  use sps_vars, only: ntfull, time_full
  use sps_utils, only: sfhint, sfhlimits
  implicit none
  integer, intent(in) :: sfh_type
  integer :: i
  real(SP) dimension(2) :: tlim
  real(SP) :: dt 
  real(SP), dimension(ntfull) :: left=0., right=0.
  real(SP), intent(out), dimension(ntfull) ::sfh_weight=0.
  
  ! Loop over each SSP and calculate its weight in the given sfh
  do i=1,ntfull
     if (i.gt.1) then
        ! There is a younger (`left`) bin, and we calculate its contribution to
        ! the weight.
        ! First calculate actual limits for the younger bin.  
        tlim = sfhlimits(time_full(i-1), time_full(i))
        ! The elements of `tlim` will be equal if there is no valid SFR in the
        ! younger bin; only proceed if there is a non-zero sfr in the younger bin
        if (tlim(1).ne.tlim(2)) then
           dt = time_full(i) - time_full(i-1)
           ! Note sign flip here
           left(i) = 0 - dt * (sfhint(i-1, tlim(1)) - sfh_int(i-1, tlim(2)))
        endif
     endif
     if (i.lt.ntfull) then
        ! There is an older (`right`) bin, we calculate its contribution to the weight
        tlim = limits(time_full(i), time_full(i+1))
        if (tlim(1).ne.tlim(2)) then
           dt = time_full(i+1) - time_full(i)
           right(i) = dt * (sfhint(i+1, tlim(1)) - sfh_int(i+1, tlim(2)))
        endif
     endif
  enddo
  
  sfh_weight = right + left
