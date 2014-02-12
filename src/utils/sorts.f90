! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2013-2014, Schmidt


module sorts
  use swaps
  use quicksort_utils
  implicit none
  
  contains
  
  ! quick sort
  !  inputs:
  !    x = array to be sorted
  !    l = left most index of x
  !    r = right most index of x
  recursive subroutine quicksort(x, l, r)
    ! inputs
    integer, intent(in) :: l, r
    double precision, intent(inout) :: x(*)
    ! local
    integer :: itmp, ind
    double precision :: pvt
    double precision :: tmp(3)
    
    
    if (l < r) then
      ! choose median of l, r, and (l+r)/2 as pivot
      tmp = (/ l, (l+r)/2, r /)
      pvt = quicksort_median_of_3(tmp)
      
      ! partition x by pvt
      if (tmp(1) == pvt) then
        ind = l 
      else if (tmp(2) == pvt) then
        ind = (l+r)/2
      else
        ind = r
      end if
      
      call swap(x, ind, r)
      
      call quicksort_partition(x, l, r, ind, itmp)
      
      ! sort recursively
      call quicksort(x, l, itmp)
      call quicksort(x, itmp+2, r)
    end if
    
    return
  end subroutine
  
  
end module

