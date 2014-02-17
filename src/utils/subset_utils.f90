! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2014, Schmidt


module subset_utils
  implicit none
  
  
  contains
  
  function subsetter_find_indx_sign(arr, larr) &
  result(indx_sign)
    use signs
    ! in/out
    integer, intent(in) :: larr
    integer, intent(in) :: arr(larr)
    integer :: indx_sign
    ! local
    integer :: i
    
    
    if (larr > 0) then
      indx_sign = signof(arr(1))
      
      do i = 2, larr
        ! mixed signs
        if (signof(arr(i)) /= indx_sign) then
          indx_sign = -2
          return
        end if
      end do
    else
      indx_sign = 0
    end if
    
    return
  end function
  
  
  
  function subsetter_check_indices(indx_sign, limit, arr, larr, errno) &
  result(info)
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
  
  
  
  subroutine subsetter_ydims(m, n, lrows, lcols, row_indx_sign, col_indx_sign, my, ny)
    ! in/out
    integer, intent(in) :: m, n
    integer, intent(in) :: lrows, lcols
    integer, intent(in) :: row_indx_sign, col_indx_sign
    integer, intent(out) :: my, ny
    
    
    !!! Rows
    if (lrows == -1) then
      my = m
    else if (row_indx_sign == 1) then
      my = lrows
    else
      my = m - lrows
    end if
    
    !!! Columns
    if (lcols == -1) then
      ny = n
    else if (col_indx_sign == 1) then
      ny = lcols
    else
      ny = n - lcols
    end if
    
    return
  end subroutine
  
  
  
  subroutine subsetter_negind_to_posind(ret, lret, arr, larr)
    use sorts
    ! in/out
    integer, intent(in) :: lret, larr
    integer, intent(in) :: arr(larr)
    integer, intent(out) :: ret(lret)
    ! local
    integer :: i, j, ind, tmp
    integer, allocatable :: arr_cp(:)
    
    
    j = 1
    ind = 0
    
    allocate(arr_cp(larr))
    
    do i = 1, larr
      arr_cp(i) = abs(arr(i))
    end do
    
    call quicksort(arr_cp, larr)
    
    tmp = arr_cp(1)
    
    do i = 1, lret
      1 continue
        ind = ind + 1
        ! Fill positive indices
        if (ind /= tmp) then
          ret(i) = ind
        ! Skip over negative ones
        else
          j = j + 1
          tmp = arr_cp(j)
          goto 1
        end if
    end do
    
    deallocate(arr_cp)
    
    return
  end subroutine
  
end module
