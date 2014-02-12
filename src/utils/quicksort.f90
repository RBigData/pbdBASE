! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2013-2014, Schmidt


module quicksort_utils
  use swaps
  implicit none
  
  
  contains
  
  ! in-place partition function for quicksort
  !  inputs:
  !    x = vector from which median is to be chosen
  !    l = leftmost index
  !    r = rightmost index
  !    ind = index of pivot
  !  output:
  !    ret = number of elements in x <= pvt
  subroutine quicksort_partition(x, l, r, ind, ret)
    ! In/out
    integer, intent(in) :: l, r, ind
    double precision, intent(inout) :: x(*)
    integer, intent(out) :: ret
    ! local
    integer :: i
    double precision :: pvt
    
    
    ret = l
    
    pvt = x(ind)
    
    call swap(x, ind, r)
    
    do i = l, r
      if (x(i).lt.pvt) then
        call swap(x, i, ret)
        ret = ret + 1
      end if
    end do
    
    call swap(x, ret, r)
    
    ret = ret - 1 ! for "<=" return
    
    return
  end subroutine
  
  
  
  function quicksort_median_of_3(tmp) result(med)
    ! in/out
    double precision, intent(inout) :: tmp(3)
    double precision :: med
    ! local
    integer :: i
    
    
    ! bubble sort with a poor man's while loop
    1 continue
      do i = 1, 2
        if (tmp(i) > tmp(i+1)) then
            call swap(tmp, i, i+1)
            goto 1
        end if
      end do
    
    med = tmp(2)
    
    return
  end function
  
end module
