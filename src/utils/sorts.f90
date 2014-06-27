! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2013-2014, Schmidt


! I really wish I had c++ templates right now ;_;
module sorts
  use swaps
  use quicksort_utils
  implicit none
  
  
  interface insertionsort
    module procedure iinsertionsort, sinsertionsort, dinsertionsort
  end interface
  
  interface quicksort
    module procedure iquicksort, squicksort, dquicksort
  end interface
  
  
  contains
  
  
  ! --------------------------------------------------------
  ! Insertion sorts
  ! --------------------------------------------------------
  subroutine iinsertionsort(x, xlen)
    ! In/out
    integer, intent(in) :: xlen
    integer, intent(inout) :: x(*)
    ! local
    integer :: i, j
    integer :: tmp
    
    include 'include/insertionsort_generic.inc'
    
    return
  end subroutine
  
  
  
  subroutine sinsertionsort(x, xlen)
    ! In/out
    integer, intent(in) :: xlen
    real, intent(inout) :: x(*)
    ! local
    integer :: i, j
    real :: tmp
    
    include 'include/insertionsort_generic.inc'
    
    return
  end subroutine
  
  
  
  subroutine dinsertionsort(x, xlen)
    ! In/out
    integer, intent(in) :: xlen
    double precision, intent(inout) :: x(*)
    ! local
    integer :: i, j
    double precision :: tmp
    
    include 'include/insertionsort_generic.inc'
    
    return
  end subroutine  
  
  
  
  ! --------------------------------------------------------
  ! Quicksorts
  ! --------------------------------------------------------
  recursive subroutine iquicksort_r(x, l, r)
    ! inputs
    integer, intent(in) :: l, r
    integer, intent(inout) :: x(*)
    ! local
    integer :: itmp, ind
    integer :: pvt
    integer :: tmp(3)
    
    include 'include/quicksort_r_generic.inc'
      if (r-l <= 10) then
        call iinsertionsort(x(l), r-l+1)
      else
        ! sort recursively
        call iquicksort_r(x, l, itmp)
        call iquicksort_r(x, itmp+2, r)
      end if
    end if
    
    return
  end subroutine
  
  
  
  recursive subroutine squicksort_r(x, l, r)
    ! inputs
    integer, intent(in) :: l, r
    real, intent(inout) :: x(*)
    ! local
    integer :: itmp, ind
    real :: pvt
    real :: tmp(3)
    
    include 'include/quicksort_r_generic.inc'
      if (r-l <= 10) then
        call sinsertionsort(x(l), r-l+1)
      else
        ! sort recursively
        call squicksort_r(x, l, itmp)
        call squicksort_r(x, itmp+2, r)
      end if
    end if
    
    return
  end subroutine
  
  
  
  recursive subroutine dquicksort_r(x, l, r)
    ! inputs
    integer, intent(in) :: l, r
    double precision, intent(inout) :: x(*)
    ! local
    integer :: itmp, ind
    double precision :: pvt
    double precision :: tmp(3)
    
    include 'include/quicksort_r_generic.inc'
      if (r-l <= 10) then
        call dinsertionsort(x(l), r-l+1)
      else
        ! sort recursively
        call dquicksort_r(x, l, itmp)
        call dquicksort_r(x, itmp+2, r)
      end if
    end if
    
    return
  end subroutine
  
  
  
  subroutine iquicksort(x, xlen)
    ! inputs
    integer, intent(in) :: xlen
    integer, intent(inout) :: x(*)
    
    call iquicksort_r(x, 1, xlen)
    
  end subroutine
  
  
  
  subroutine squicksort(x, xlen)
    ! inputs
    integer, intent(in) :: xlen
    real, intent(inout) :: x(*)
    
    call squicksort_r(x, 1, xlen)
    
  end subroutine
  
  
  
  subroutine dquicksort(x, xlen)
    ! inputs
    integer, intent(in) :: xlen
    double precision, intent(inout) :: x(*)
    
    call dquicksort_r(x, 1, xlen)
    
  end subroutine
  
end module

