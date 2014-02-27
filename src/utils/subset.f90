! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2014, Schmidt


! Subroutine to mimic R's `[` operator
! lrows = -1 : all rows are used and 'rows' is ignored
! lcols = -1 : all cols are used and 'cols' is ignored
module subset
  use subset_utils
  use subset_special
  use iso_c_binding
  
  implicit none
  
  
  contains
  
  subroutine subsetter(m, n, x, rows, lrows, cols, lcols, y, my, ny, info) &
  bind(c, name='subsetter_')
    ! in/out
    integer, intent(in) :: m, n, lrows, lcols
    integer, intent(in) :: my, ny
    integer, intent(in) :: rows(*), cols(*)
    integer, intent(out) :: info
    double precision, intent(in) :: x(m, n)
    double precision, intent(out) :: y(my, ny)
    ! local
    integer :: i, j
    integer :: row_indx_sign, col_indx_sign
    
    
    !!! Quick return if possible
    if (lrows == 0) then
      info = -5
      return
    else if (lcols == 0) then
      info = -7
      return
    else if (lrows == -1 .and. lcols == -1) then
      print *, y
!      allocate(y(m, n))
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
    
    
    !!! Determine dims and allocated
!    call subsetter_ydims(m, n, lrows, lcols, row_indx_sign, col_indx_sign, my, ny)
!    allocate(y(my,ny))
    
    
    !!! Special cases
    if (lrows == -1) then
      if (col_indx_sign == 1) then
        call subsetter_allrows_posind(m, n, x, cols, lcols, y, my, ny, info)
        return
      else if (col_indx_sign == -1) then
        call subsetter_allrows_negind(m, n, x, cols, lcols, y, my, ny, info)
        return
      end if
      return
    else if (lcols == -1) then
      if (row_indx_sign == 1) then
        call subsetter_allcols_posind(m, n, x, rows, lrows, y, my, ny, info)
        return
      else if (row_indx_sign == -1) then
        call subsetter_allcols_negind(m, n, x, rows, lrows, y, my, ny, info)
        return
      end if
      return
    end if
    
    
    !!! General case
    if (row_indx_sign == 1 .and. col_indx_sign == 1) then
      call subsetter_posind(m, n, x, rows, lrows, cols, lcols, y, my, ny, info)
!    else if 
      
    end if
    
    
    info = 0
    return
  end subroutine
  
  
  
  subroutine subsetter_posind(m, n, x, rows, lrows, cols, lcols, y, my, ny, info)
    implicit none
    ! in/out
    integer, intent(in) :: m, n, lrows, lcols, my, ny
    integer, intent(in) :: rows(*), cols(*)
    integer, intent(out) :: info
    double precision, intent(in) :: x(m, n)
    double precision, intent(out) :: y(my, ny)
    ! local
    integer :: i, j, col
    
    
    !!! Build y
    do j = 1, ny
     col = cols(j)
      do i = 1, my
        y(i, j) = x(rows(i), col)
      end do
    end do
    
    info = 0
    return
  end subroutine
  
end module

