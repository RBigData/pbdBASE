! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2014, Schmidt


subroutine dmksym(triang, m, n, x)
  implicit none
  ! in/out
  character(len=1), intent(in) :: triang
  integer, intent(in) :: m, n
  double precision, intent(inout) :: x(m, n)
  ! local
  integer :: i, j
  
  
  ! copy to upper from lower
  if (triang == 'l') then
    do j = 1, n
      do i = j+1, m
        x(j, i) = x(i, j)
      end do
    end do
  
  ! copy to lower from upper
  else if (triang == 'u') then
    do j = 1, n
      do i = 1, j-1
        x(j, i) = x(i, j)
      end do
    end do 
  else
    return
  end if
  
  return
end


