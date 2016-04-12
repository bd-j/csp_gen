function sfh_weight(sfh, imin, imax)
  ! Function to calculate the weights
  !
  !
  use sps_vars, only: ntfull, time_full, pset
  use sps_utils, only: sfhint, sfhlimits
  implicit none
  
  type(SFHPARAMS), intent(in) :: sfh
  integer, intent(in) :: imin, imax
  
  real(SP), intent(out), dimension(ntfull) ::sfh_weight=0.

  integer :: i
  real(SP) dimension(2) :: tlim
  real(SP) :: dt 
  real(SP), dimension(ntfull) :: left=0., right=0.
  
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
           dt = time_full(i) - time_full(i-1)
           ! Note sign flip here
           left(i) = 0 - dt * (sfhint(i-1, tlim(1), sfh) - sfhint(i-1, tlim(2), sfh))
        endif
     endif
     if (i.lt.ntfull) then
        ! There is an older (`right`) bin, we calculate its contribution to the weight
        tlim(1) = sfhlimit(time_full(i), sfh)
        tlim(2) = sfhlimit(time_full(i+1), sfh)
        ! The elements of `tlim` will be equal if there is no valid SFR in the
        ! younger bin; only proceed if there is a non-zero sfr in the older bin
        if (tlim(1).ne.tlim(2)) then
           dt = time_full(i+1) - time_full(i)
           right(i) = dt * (sfhint(i+1, tlim(1), sfh) - sfhint(i+1, tlim(2), sfh))
        endif
     endif
  enddo
  
  sfh_weight = right + left

  ! need to add weights from the zero bin
  
end function sfh_weight
