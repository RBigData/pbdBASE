! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2013, Schmidt


! This will not compile with ifortran, or probably anything else that isn't
! gfortran.  Used for profiling; you have no reason to ever use this.

! Call once with ONOFF = 1 to start, then call again with ONOFF = 0 to 
! halt and print the wallclock runtime.
  ! If ONOFF = 1, then TIME is input
  ! If ONOFF = 0, then TIME is output.
subroutine timer(time, onoff)
  implicit none
  ! in/out
  integer             time, onoff
  ! local
  double precision    timeold
  
  
  ! turn it on
  if (onoff >= 1) then
    call system_clock(time)
  ! turn it off
  else
    timeold = time
    call system_clock(time)
    time = time - timeold
!    if (blacs(4) == 0 .and. blacs(5) == 0) then
!      write (*,*) "elapsed: ", real(time)/1000.0
!    end if
  end if
  
  return
end subroutine


