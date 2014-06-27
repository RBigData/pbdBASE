! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2014, Schmidt


module string_tools
  implicit none
  
  contains
  
  
  ! Makes string all lowercase based on ascii encoding
  function toupper(str) &
  result(ret)
    implicit none
    ! in/out
    character(len=*) :: str
    character(len=len(str)) :: ret
    ! Local
    integer :: i, ascii_i
    integer :: ascii_low = iachar("a")
    integer :: ascii_high = iachar("z")
    
    do i = 1, len(str)
      ascii_i = iachar(str(i:i))
        if (ascii_low <= ascii_i .and. ascii_i <= ascii_high) then
          ret(i:i) = achar(iachar(str(i:i)) - 32)
        else
          ret(i:i) = str(i:i)
        end if
    end do
    
    return
  end function
  
  
  
  ! Makes string all uppercase based on ascii encoding
  function tolower(str) &
  result(ret)
    implicit none
    ! in/out
    character(len=*) :: str
    character(len=len(str)) :: ret
    ! Local
    integer :: i, ascii_i
    integer :: ascii_low = iachar("A")
    integer :: ascii_high = iachar("Z")
    
    do i = 1, len(str)
      ascii_i = iachar(str(i:i))
        if (ascii_low <= ascii_i .and. ascii_i <= ascii_high) then
          ret(i:i) = achar(iachar(str(i:i)) + 32)
        else
          ret(i:i) = str(i:i)
        end if
    end do
    
    return
  end function
  
  
end module


