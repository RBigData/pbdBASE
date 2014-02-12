! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2014, Schmidt


function subsetter_find_indx_sign(arr, larr) result(indx_sign)
  use signs
  implicit none
  ! in/out
  integer, intent(in) :: larr
  integer, intent(in) :: arr(larr)
  integer :: indx_sign
  ! local
  integer :: i
  ! Functions
  
  
  if (larr > 0) then
    indx_sign = signof(arr(1))
    
    do i = 2, larr
      ! mixed signs
      if (signof(arr(i)) /= indx_sign) then
        indx_sign = -2
        return
      end if
    end do
  end if
  
  indx_sign = 0
  
  return
end function



function subsetter_check_indices(indx_sign, limit, arr, larr, errno) result(info)
  implicit none
  ! in/out
  integer, intent(in) :: limit, larr, indx_sign, errno
  integer, intent(in) :: arr(larr)
  integer :: info
  ! local
  integer :: i
  
  if (indx_sign == 1) then
    do i = 1, larr
      if (arr(i) > limit) then
        info = errno*10
        return
      end if
    end do
  else if (indx_sign == -1) then
    do i = 1, larr
      if (arr(i) < -limit) then
        info = errno*10
        return
      end if
    end do
  end if
  
  info = 0
  
  return
end function



! Subroutine to mimic R's `[` operator
!FIXME just doing positive indexing for now
! lrows = -1 : all rows are used and 'rows' is ignored
! lcols = -1 : all cols are used and 'cols' is ignored
subroutine subsetter(m, n, x, rows, lrows, cols, lcols, y, info)
  implicit none
  ! in/out
  integer, intent(in) :: m, n, lrows, lcols
  integer, intent(in) :: rows(*), cols(*)
  integer, intent(out) :: info
  double precision, intent(in) :: x(m, n)
  double precision, intent(out) :: y(m, n)
  ! local
  integer :: i, j
  integer :: row_indx_sign, col_indx_sign
  ! functions
  integer :: subsetter_find_indx_sign, subsetter_check_indices
  
  
  !!! Quick return if possible
  if (lrows == 0 .or. lcols == 0) then
    info = 0
    return
  else if (lrows == -1 .and. lcols == -1) then
    y = x
    info = 0
    return
  end if
  
  
  !!! Input checks
  if (lrows < -1) then
    info = -5
    return
  else if (lcols < -1) then
    info = -7
    return
  end if
  
  
  ! check that indices are all positive or all negative, and record sign
  row_indx_sign = subsetter_find_indx_sign(rows, lrows)
  if (row_indx_sign == -2) then
    info = -4
    return
  end if
  
  col_indx_sign = subsetter_find_indx_sign(cols, lcols)
  if (col_indx_sign == -2) then
    info = -6
    return
  end if
  
  
  ! check that indices are sane
  info = subsetter_check_indices(row_indx_sign, m, rows, lrows, -4)
  if (info /= 0) return
  info = subsetter_check_indices(col_indx_sign, n, cols, lcols, -6)
  if (info /= 0) return
  
  
  !!! Special cases
  if (lrows == -1) then
    if (col_indx_sign == 1) then
      call subsetter_allrows_posind(m, n, x, cols, lcols, y, info)
      return
    else if (col_indx_sign == -1) then
!      call subsetter_allrows_negind(m, n, x, cols, lcols, y, info)
      return
    end if
    return
  else if (lcols == -1) then
    if (row_indx_sign == 1) then
      call subsetter_allcols_posind(m, n, x, rows, lrows, y, info)
      return
    else if (row_indx_sign == -1) then
!      call subsetter_allcols_negind(m, n, x, rows, lrows, y, info)
      return
    end if
    return
  end if
  
  
  !!! General case
  
  
  info = 0
  return
end subroutine



subroutine subsetter_negind_to_posind(ret, lret, arr, larr)
  use sorts
  ! in/out
  integer, intent(in) :: lret, larr
  integer, intent(in) :: arr(larr)
  integer, intent(out) :: ret(lret)
  ! local
  integer :: i, j
  
  
  j = 1
  
!  call quicksort(arr, 1, larr)
  
!  do i = 1, maxlen
!    if () then
!      
!    else
!      
!    end if
!  end do
  
  
  return
end subroutine



subroutine subsetter_allcols_posind(m, n, x, rows, lrows, y, info)
  implicit none
  ! in/out
  integer, intent(in) :: m, n, lrows
  integer, intent(in) :: rows(*)
  integer, intent(out) :: info
  double precision, intent(in) :: x(m, n)
  double precision, intent(out) :: y(m, n)
  ! local
  integer :: i, j
  
  
  !!! Build y
  do j = 1, n
    do i = 1, lrows
      y(i, j) = x(rows(i), j)
    end do
  end do
  
  info = 0
  return
end subroutine



!subroutine subsetter_allcols_negind(m, n, x, rows, lrows, y, info)
!  implicit none
!  ! in/out
!  integer, intent(in) :: m, n, lrows
!  integer, intent(in) :: rows(*)
!  integer, intent(out) :: info
!  double precision, intent(in) :: x(m, n)
!  double precision, intent(out) :: y(m, n)
!  ! local
!  integer :: i, j
!  integer :: ltmp = m - lrows
!  double precision, allocatable :: tmp(:)
!  
!  
!  info = 0
!  
!  allocate(tmp(ltmp))
!  do i = 1, m
!    if () then
!      
!    else
!      
!    end if
!  end do
!  
!  call subsetter_allcols_posind(m, n, x, tmp, ltmp, y, info)
!  
!  
!  deallocate(tmp)
!  return
!end subroutine



subroutine subsetter_allrows_posind(m, n, x, cols, lcols, y, info)
  implicit none
  ! in/out
  integer, intent(in) :: m, n, lcols
  integer, intent(in) :: cols(*)
  integer, intent(out) :: info
  double precision, intent(in) :: x(m, n)
  double precision, intent(out) :: y(m, n)
  ! local
  integer :: i, j
  
  
  !!! Build y
  do j = 1, n
    do i = 1, lcols
      y(i, j) = x(i, cols(j))
    end do
  end do
  
  info = 0
  return
end subroutine


