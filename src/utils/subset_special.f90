! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2014, Schmidt


module subset_special
  implicit none
  
  
  contains
  
  subroutine subsetter_allrows_posind(m, n, x, cols, lcols, y, my, ny, info)
    ! in/out
    integer, intent(in) :: m, n, lcols, my, ny
    integer, intent(in) :: cols(*)
    integer, intent(out) :: info
    double precision, intent(in) :: x(m, n)
    double precision, intent(out) :: y(my, ny)
    ! local
    integer :: i, j, col
    
    
    !!! Build y
    do j = 1, ny
     col = cols(j)
      do i = 1, my
        y(i, j) = x(i, col)
      end do
    end do
    
    info = 0
    return
  end subroutine
  
  
  
  subroutine subsetter_allcols_posind(m, n, x, rows, lrows, y, my, ny, info)
    ! in/out
    integer, intent(in) :: m, n, lrows, my, ny
    integer, intent(in) :: rows(*)
    integer, intent(out) :: info
    double precision, intent(in) :: x(m, n)
    double precision, intent(out) :: y(my, ny)
    ! local
    integer :: i, j
    
    
    !!! Build y
    do j = 1, ny
      do i = 1, my
        y(i, j) = x(rows(i), j)
      end do
    end do
    
    info = 0
    return
  end subroutine
  
  
  
  subroutine subsetter_allrows_negind(m, n, x, cols, lcols, y, my, ny, info)
    use sorts
    ! in/out
    integer, intent(in) :: m, n, lcols, my, ny
    integer, intent(in) :: cols(*)
    integer, intent(out) :: info
    double precision, intent(in) :: x(m, n)
    double precision, intent(out) :: y(my, ny)
    ! local
    integer :: i, j, k
    integer :: col, ind, tmp
    integer, allocatable :: cols_cp(:)
    
    
    allocate(cols_cp(lcols))
    
    do i = 1, lcols
      cols_cp(i) = abs(cols(i))
    end do
    
    call quicksort(cols_cp, lcols)
    
    tmp = cols_cp(1)
    
    
    do j = 1, ny
      col = cols_cp(j)
      k = 1
      ind = 0
      do i = 1, my
        1 continue
          ind = ind + 1
          ! Fill positive indices
          if (ind /= tmp) then
            y(i, j) = x(i, j)
          ! Skip over negative ones
          else
            k = k + 1
            tmp = cols_cp(k)
            goto 1
          end if
      end do
    end do
    
    
    info = 0
    return
  end subroutine
  
  
  
  subroutine subsetter_allcols_negind(m, n, x, rows, lrows, y, my, ny, info)
    use sorts
    ! in/out
    integer, intent(in) :: m, n, lrows, my, ny
    integer, intent(in) :: rows(*)
    integer, intent(out) :: info
    double precision, intent(in) :: x(m, n)
    double precision, intent(out) :: y(my, ny)
    ! local
    integer :: i, j, k
    integer :: row, ind, tmp
    integer, allocatable :: rows_cp(:)
    
    
    allocate(rows_cp(lrows))
    
    do i = 1, lrows
      rows_cp(i) = abs(rows(i))
    end do
    
    call quicksort(rows_cp, lrows)
    
    tmp = rows_cp(1)
    
    
    do j = 1, ny
      k = 1
      ind = 0
      do i = 1, my
        row = rows_cp(i)
        1 continue
          ind = ind + 1
          ! Fill positive indices
          if (ind /= tmp) then
            y(i, j) = x(i, j)
          ! Skip over negative ones
          else
            k = k + 1
            tmp = rows_cp(k)
            goto 1
          end if
      end do
    end do
    
    
    info = 0
    return
  end subroutine
  
end module

