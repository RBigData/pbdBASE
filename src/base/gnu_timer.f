! This will not compile with ifortran, or probably anything else that isn't
! gfortran.  Used for profiling; you have no reason to ever use this.
!
! Call once with ONOFF = 1 to start, then call again with ONOFF = 0 to 
! halt and print the wallclock runtime.
  ! If ONOFF = 1, then TIME is input
  ! If ONOFF = 0, then TIME is output.
      SUBROUTINE TIMER(TIME, ONOFF)
      IMPLICIT NONE
      ! IN/OUT
      INTEGER             TIME, ONOFF
      ! Local
      DOUBLE PRECISION    TIMEOLD
      
      
      ! Turn it on
      IF (ONOFF.GE.1) THEN
        CALL SYSTEM_CLOCK(TIME)
      ! Turn it off
      ELSE
        TIMEOLD = TIME
        CALL SYSTEM_CLOCK(TIME)
        TIME = TIME - TIMEOLD
!        IF (BLACS(4).EQ.0 .AND. BLACS(5).EQ.0) THEN
!          WRITE (*,*) "Elapsed: ", REAL(TIME)/1000.0
!        END IF
      END IF
      
      RETURN
      END
