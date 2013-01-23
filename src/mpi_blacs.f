! Initialize BLACS communicator, process grid
      SUBROUTINE mpi_blacs_initialize(NPROW, NPCOL, ICTXT, MYROW, MYCOL)
      ! Inputs
      INTEGER NPROW, NPCOL
      ! Outputs
      INTEGER ICTXT, MYROW, MYCOL
      ! External
      EXTERNAL SL_INIT, BLACS_GRIDINFO
      
      CALL SL_INIT( ICTXT, NPROW, NPCOL )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )

      RETURN
      END


!      SUBROUTINE DIMS(DESC, LDM, BLACS)
!      IMPLICIT NONE
!      ! IN/OUT
!      INTEGER             DESC( 9 ), LDM( 2 ), BLACS(4)
!      ! Functions
!      INTEGER             NUMROC
!      ! External
!      EXTERNAL            BLACS_GRIDINFO
!      
!      CALL BLACS_GRIDINFO(DESC(2), BLACS(1), BLACS(2), 
!     $                    BLACS(3), BLACS(4))
!      
!      WRITE (*,*) BLACS
!      
!      LDM(1) = NUMROC(DESC(3), DESC(5), BLACS(3), 
!     $                DESC(7), BLACS(1))
!      LDM(2) = NUMROC(DESC(4), DESC(6), BLACS(4), 
!     $                DESC(8), BLACS(2))
!      
!!      IF (LDM(1).LT.1) THEN
!!        LDM(1) = 0
!!      END IF
!!      
!!      IF (LDM(2).LT.1) THEN
!!        LDM(2) = 0
!!      END IF
!      
!      RETURN
!      END
