module string_tools
  implicit none
  
  contains
  
  
  ! Makes string all lowercase based on ascii encoding
  function toupper(str) result(ret)
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
    
  end function
  
  
  
  ! Makes string all uppercase based on ascii encoding
  function tolower(str) result(ret)
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
    
  end function
  
  
end module


