! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2013-2014, Schmidt


module quicksort_utils
  use swaps
  implicit none
  
  
  interface quicksort_median_of_3
    module procedure iquicksort_median_of_3, squicksort_median_of_3, dquicksort_median_of_3
  end interface
  
  interface quicksort_partition
    module procedure iquicksort_partition, squicksort_partition, dquicksort_partition
  end interface
  
  
  contains
  
  ! --------------------------------------------------------
  ! Median of a length 3 array; used in choosing pivot in quicksort
  ! --------------------------------------------------------
  function iquicksort_median_of_3(tmp) result(med)
    ! in/out
    integer, intent(inout) :: tmp(3)
    integer :: med
    ! local
    integer :: i
    
    include 'include/quicksort_median_of_3_generic.inc'
    
    return
  end function
  
  
  
  function squicksort_median_of_3(tmp) result(med)
    ! in/out
    real, intent(inout) :: tmp(3)
    real :: med
    ! local
    integer :: i
    
    include 'include/quicksort_median_of_3_generic.inc'
    
    return
  end function
  
  
  
  function dquicksort_median_of_3(tmp) result(med)
    ! in/out
    double precision, intent(inout) :: tmp(3)
    double precision :: med
    ! local
    integer :: i
    
    include 'include/quicksort_median_of_3_generic.inc'
    
    return
  end function
  
  
  ! --------------------------------------------------------
  ! in-place partition function for quicksort
  !  inputs:
  !    x = vector from which median is to be chosen
  !    l = leftmost index
  !    r = rightmost index
  !    ind = index of pivot
  !  output:
  !    ret = number of elements in x <= pvt
  subroutine iquicksort_partition(x, l, r, ind, ret)
    ! In/out
    integer, intent(in) :: l, r, ind
    integer, intent(inout) :: x(*)
    integer, intent(out) :: ret
    ! local
    integer :: i
    integer :: pvt
    
    include 'include/quicksort_partition_generic.inc'
    
    return
  end subroutine
  
  
  
  subroutine squicksort_partition(x, l, r, ind, ret)
    ! In/out
    integer, intent(in) :: l, r, ind
    real, intent(inout) :: x(*)
    integer, intent(out) :: ret
    ! local
    integer :: i
    real :: pvt
    
    include 'include/quicksort_partition_generic.inc'
    
    return
  end subroutine
  
  
  
  subroutine dquicksort_partition(x, l, r, ind, ret)
    ! In/out
    integer, intent(in) :: l, r, ind
    double precision, intent(inout) :: x(*)
    integer, intent(out) :: ret
    ! local
    integer :: i
    double precision :: pvt
    
    include 'include/quicksort_partition_generic.inc'
    
    return
  end subroutine
  
end module
