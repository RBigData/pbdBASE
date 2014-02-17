! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2014, Schmidt


! A "not dumb" version of sign.
module signs
  implicit none
  
  interface signof
    module procedure isignof, ssignof, dsignof
  end interface
  
  contains
  
  
  function isignof(x) result(sgn)
    implicit none
    integer, intent(in) :: x
    integer :: sgn
    
    
    if (x > 0) then
      sgn = 1
    else if (x < 0) then
      sgn = -1
    else
      sgn = 0
    end if
    
    return
  end function
  
  
  
  function ssignof(x) result(sgn)
    implicit none
    real, intent(in) :: x
    integer :: sgn
    
    
    if (x > 0) then
      sgn = 1
    else if (x < 0) then
      sgn = -1
    else
      sgn = 0
    end if
    
    return
  end function
  
  
  
  function dsignof(x) result(sgn)
    implicit none
    double precision, intent(in) :: x
    integer :: sgn
    
    
    if (x > 0) then
      sgn = 1
    else if (x < 0) then
      sgn = -1
    else
      sgn = 0
    end if
    
    return
  end function
  
end module
