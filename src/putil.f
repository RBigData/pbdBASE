! Copyright 2013, Schmidt

! DESC = ScaLAPACK descriptor array
! DM = [DM1, DM2] = Global dimension
! LDM = [LDM1, LDM2] = Local dimension
! BLACS = [NPROCS, NPROW, NPCOL, MYPROW, MYPCOL]
      SUBROUTINE PDIMS(DESC, LDM, BLACS)
      IMPLICIT NONE
      ! IN/OUT
      INTEGER             DESC( 9 ), LDM( 2 ), BLACS( 5 )
      ! Functions
      INTEGER             NUMROC
      ! External
      EXTERNAL            BLACS_GRIDINFO
      
      
      CALL BLACS_GRIDINFO(DESC(2), BLACS(2), BLACS(3), 
     $                    BLACS(4), BLACS(5))
      
      IF (BLACS(2).EQ.(-1) .OR. BLACS(3).EQ.(-1)) THEN
        BLACS(1) = -1
      ELSE
        BLACS(1) = BLACS(2) * BLACS(3)
      END IF
      
      LDM(1) = NUMROC(DESC(3), DESC(5), BLACS(4), 
     $                DESC(7), BLACS(2))
      LDM(2) = NUMROC(DESC(4), DESC(6), BLACS(5), 
     $                DESC(8), BLACS(3))
      
      IF (LDM(1).LT.1) LDM(1) = 0
      IF (LDM(2).LT.1) LDM(2) = 0
      
      RETURN
      END

! Simple timing routine 
      SUBROUTINE TIMER(TIME, ONOFF)
      ! IN/OUT
      INTEGER             TIME, ONOFF
      ! Local
      DOUBLE PRECISION    TIMEOLD
      
      
      ! Turn it on
      IF (ONOFF.GE.1) THEN
        CALL SYSTEM_CLOCK(TIME)
      ELSE
        TIMEOLD = TIME
        CALL SYSTEM_CLOCK(TIME)
        TIME = TIME - TIMEOLD
!        IF (BLACS(4).EQ.0 .AND. BLACS(5).EQ.0) THEN
          WRITE (*,*) "Elapsed: ", REAL(TIME)/1000.0
!        END IF
      END IF
      
      RETURN
      END

! Local-to-global pair of indices; shorthand for calling INDXL2G twice.
      SUBROUTINE L2GPAIR(I, J, GI, GJ, DESC, BLACS)
      IMPLICIT NONE
      ! IN/OUT
      INTEGER             I, J, GI, GJ, DESC(9), BLACS(2)
      ! Functions
      INTEGER             INDXL2G
      
      
      GI = INDXL2G(I, DESC(5), BLACS(4), 0, BLACS(2))
      GJ = INDXL2G(J, DESC(6), BLACS(5), 0, BLACS(3))
      
      RETURN
      END


! Zero out the (global) triangle of a distributed matrix
      SUBROUTINE PTRI2ZERO(UPLO, DIAG, X, DESCX)
      IMPLICIT NONE
      ! IN/OUT
      INTEGER             IX, JX, DESCX(9)
      DOUBLE PRECISION    X(DESCX(9), *)
      CHARACTER*1         UPLO, DIAG
      ! Local
      INTEGER              M, N, I, J, GI, GJ, 
     $                    LDM(2), BLACS(5)
      ! Parameter
      DOUBLE PRECISION    ZERO
      PARAMETER ( ZERO = 0.0D0 )
      ! External
      EXTERNAL            PDIMS, L2GPAIR
      
      
      ! Get local and proc grid info
      CALL PDIMS(DESCX, LDM, BLACS)
      
      M = LDM(1)
      N = LDM(2)
      
      ! Only do work if we own any local pieces
      IF (M.GT.0 .AND. N.GT.0) THEN
      
      ! Zeroing both triangles
        IF (UPLO.EQ.'B') THEN
          ! Zero the diagonal as well
          IF (DIAG.EQ.'Y') THEN
            X(:,1:N) = ZERO
          ! Leave the diagonal alone
          ELSE IF (DIAG.EQ.'N') THEN
            DO J = 1, N
              DO I = 1, M
                CALL L2GPAIR(I, J, GI, GJ, DESCX, BLACS)
                IF (GI.NE.GJ) THEN
                  X(I,J) = ZERO
                END IF
              END DO 
            END DO
          ELSE
            STOP "Invalid argument 'DIAG'"
          END IF
        
        ! Zeroing the upper triangle
        ELSE IF (UPLO.EQ.'U') THEN
          ! Zero the diagonal as well
          IF (DIAG.EQ.'Y') THEN 
            DO J = 1, N
              DO I = 1, M
                CALL L2GPAIR(I, J, GI, GJ, DESCX, BLACS)
                IF (GI.GE.GJ) THEN
                  X(I,J) = ZERO
                END IF
              END DO 
            END DO
          ! Leave the diagonal alone
          ELSE IF (DIAG.EQ.'N') THEN 
            DO J = 1, N
              DO I = 1, M
                CALL L2GPAIR(I, J, GI, GJ, DESCX, BLACS)
                IF (GI.GT.GJ) THEN
                  X(I,J) = ZERO
                END IF
              END DO 
            END DO
          ELSE
            STOP "Invalid argument 'DIAG'"
          END IF
        
        ! Zeroing the lower triangle
        ELSE IF (UPLO.EQ.'L') THEN
          ! Zero the diagonal as well
          IF (DIAG.EQ.'Y') THEN 
            DO J = 1, N
              DO I = 1, M
                CALL L2GPAIR(I, J, GI, GJ, DESCX, BLACS)
                IF (GI.LE.GJ) THEN
                  X(I,J) = ZERO
                END IF
              END DO 
            END DO
          ! Leave the diagonal alone
          ELSE IF (DIAG.EQ.'N') THEN 
            DO J = 1, N
              DO I = 1, M
                CALL L2GPAIR(I, J, GI, GJ, DESCX, BLACS)
                IF (GI.LT.GJ) THEN
                  X(I,J) = ZERO
                END IF
              END DO 
            END DO
          ELSE
            STOP "Invalid argument 'DIAG'"
          END IF
        ELSE
          STOP "Invalid argument 'UPLO'"
        END IF
      
      END IF
      
      RETURN
      END 


! Make symmetric via copying from one (global) triangle to the other.
! UPLO = 'U', copy FROM upper
      SUBROUTINE PDMKSYM(UPLO, X, IX, JX, DESCX)
      IMPLICIT NONE
      ! IN/OUT
      INTEGER             IX, JX, DESCX(9)
      DOUBLE PRECISION    X(DESCX(9), *)
      CHARACTER*1         UPLO
      ! Local
      INTEGER             M, N, I, J, GI, GJ, LDM(2), BLACS(5)
      ! Parameter
      DOUBLE PRECISION    ZERO, ONE, HALF
      PARAMETER ( ZERO = 0.0D0, ONE = 1.0D0, HALF = 0.5D0 )
      ! External
      EXTERNAL            PDIMS, L2GPAIR, PDGEADD, PTRI2ZERO
      
      
      ! Get local and proc grid info
      CALL PDIMS(DESCX, LDM, BLACS)
      
      M = LDM(1)
      N = LDM(2)
      
      ! Only do work if we own any local pieces
      IF (M.GT.0 .AND. N.GT.0) THEN
        
        ! Zero Lower
        IF (UPLO.EQ.'U') THEN
          CALL PTRI2ZERO('U', 'N', X, DESCX)
        ! Zero Upper
        ELSE IF (UPLO.EQ.'L') THEN
          CALL PTRI2ZERO('L', 'N', X, DESCX)
        ! Otherwise...
        ELSE IF (UPLO.NE.'B') THEN
          STOP "Invalid argument 'UPLO'"
        END IF
        
        ! X = t(X) + X
        CALL PDGEADD('T', DESCX(3), DESCX(4), ONE, X, IX, JX, 
     $                DESCX, ONE, X, IX, JX, DESCX)
        
        ! Correct diagonal
        DO J = 1, N
          DO I = 1, M
            CALL L2GPAIR(I, J, GI, GJ, DESCX, BLACS)
            IF (GI.EQ.GJ) THEN
              X(I,J) = X(I,J) * HALF
            END IF
          END DO 
        END DO
        
      END IF
      
      RETURN
      END


!!!! Does what pdlacpy SHOULD do, namely copies the triangle of a global 
!!!! distributed matrix onto the triangle of another matrix.
!!!      SUBROUTINE PDLACPY2(UPLO, IP, A, IA, JA, DESCA, B, IB, JB, DESCB)
!!!      ! IN/OUT
!!!      INTEGER             IA, JA, IB, JB, DESCA(9), DESCB(9)
!!!      DOUBLE PRECISION    X(DESCX(9), *), Y(DESCY(9), *)
!!!      CHARACTER*1         UPLO, IP ! should copy of A onto B be in place?
!!!      ! Local
!!!      INTEGER             I, J, GI, GJ, LDM(2), BLACS(4), NPROCS
!!!      DOUBLE PRECISION, ALLOCATABLE, CPA(:)
!!!      ! Parameter
!!!      DOUBLE PRECISION    ZERO, ONE
!!!      PARAMETER ( ZERO = 0.0D0, ONE = 1.0D0 )
!!!      ! External
!!!      EXTERNAL            PDLACPY, PDIMS
!!!      INTEGER             INDXL2G
!!!      
!!!      
!!!      CALL PDIMS(DESCA, LDM, BLACS)
!!!      
!!!      ! Upper part to be copied
!!!      IF (UPLO.EQ.'U') THEN
!!!        
!!!        IF (IP.EQ.'N') THEN
!!!          ALLOCATE(CPA(LDM(1) * LDM(2)))
!!!          
!!!        END IF
!!!        
!!!        ! Zero upper of B
!!!        
!!!        
!!!        ! Add the two matrices
!!!        
!!!      ! Lower part to be copied
!!!      ELSE IF (UPLO.EQ.'L') THEN
!!!        
!!!        
!!!        
!!!        
!!!        
!!!      ! A = B
!!!      ELSE
!!!        CALL PDLACPY('B', M, N, A, IA, JA, DESCA, B, IB, JB, DESCB)
!!!      END IF
!!!      
!!!      
!!!      RETURN
!!!      END

