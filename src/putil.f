! Copyright 2013, Schmidt


! Optimal process grid shape
      SUBROUTINE OPTIMALGRID(NPROCS, NROWS, NCOLS)
      IMPLICIT NONE
      ! IN/OUT
      INTEGER             NPROCS, NROWS, NCOLS
      ! Local
      INTEGER             I, N
      
      
      N = INT(SQRT(REAL(NPROCS)))
      
      DO I = 0, N-1
        NROWS = N - I
        NCOLS = INT(MOD(NPROCS, NROWS))
        IF (NCOLS.EQ.0) EXIT
      END DO
      
      NCOLS = NPROCS / NROWS
      
      RETURN
      END


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
      IMPLICIT NONE
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
!          WRITE (*,*) "Elapsed: ", REAL(TIME)/1000.0
!        END IF
      END IF
      
      RETURN
      END


! Reductions using BLACS
      SUBROUTINE DALLREDUCE(M, N, X, OP, SCOPE, ICTXT)
      IMPLICIT NONE
      ! IN/OUT
      INTEGER             M, N, ICTXT
      DOUBLE PRECISION    X( * )
      CHARACTER           OP, SCOPE
      ! External
      EXTERNAL            DGSUM2D, DGAMX2D, DGMN2D
      
      
      IF (OP.EQ.'MIN') THEN
        CALL DGAMN2D(ICTXT, SCOPE, ' ', M, N, X, 1,-1,-1,-1,-1,-1)
      ELSE IF (OP.EQ.'MAX') THEN
        CALL DGAMX2D(ICTXT, SCOPE, ' ', M, N, X, 1,-1,-1,-1,-1,-1)
      ELSE
        CALL DGSUM2D(ICTXT, SCOPE, ' ', M, N, X, 1,-1,-1)
      END IF
      
      RETURN 
      END


      SUBROUTINE DREDUCE(M, N, X, OP, RDEST, CDEST, SCOPE, ICTXT)
      IMPLICIT NONE
      ! IN/OUT
      INTEGER             M, N, RDEST, CDEST, ICTXT
      DOUBLE PRECISION    X( * )
      CHARACTER           OP, SCOPE
      ! External
      EXTERNAL            DGSUM2D, DGAMX2D, DGMN2D
      
      
      IF (OP.EQ.'MIN') THEN
        CALL DGAMN2D(ICTXT, SCOPE, ' ', M, N, X, 1, -1, -1, -1, 
     $               RDEST, CDEST)
      ELSE IF (OP.EQ.'MAX') THEN
        CALL DGAMX2D(ICTXT, SCOPE, ' ', M, N, X, 1, -1, -1, -1,
     $               RDEST, CDEST)
      ELSE
        CALL DGSUM2D(ICTXT, SCOPE, ' ', M, N, X, 1, RDEST, CDEST)
      END IF
      
      RETURN 
      END


      SUBROUTINE IALLREDUCE(M, N, X, OP, SCOPE, ICTXT)
      IMPLICIT NONE
      ! IN/OUT
      INTEGER             M, N, ICTXT, X( * )
      CHARACTER           OP, SCOPE
      ! External
      EXTERNAL            IGSUM2D, IGAMX2D, IGMN2D
      
      
      IF (OP.EQ.'MIN') THEN
        CALL IGAMN2D(ICTXT, SCOPE, ' ', M, N, X, 1,-1,-1,-1,-1,-1)
      ELSE IF (OP.EQ.'MAX') THEN
        CALL IGAMX2D(ICTXT, SCOPE, ' ', M, N, X, 1,-1,-1,-1,-1,-1)
      ELSE
        CALL IGSUM2D(ICTXT, SCOPE, ' ', M, N, X, 1,-1,-1)
      END IF
      
      RETURN 
      END


      SUBROUTINE IREDUCE(M, N, X, OP, RDEST, CDEST, SCOPE, ICTXT)
      IMPLICIT NONE
      ! IN/OUT
      INTEGER             M, N, RDEST, CDEST, ICTXT
      DOUBLE PRECISION    X( * )
      CHARACTER           OP, SCOPE
      ! External
      EXTERNAL            IGSUM2D, IGAMX2D, IGMN2D
      
      
      IF (OP.EQ.'MIN') THEN
        CALL IGAMN2D(ICTXT, SCOPE, ' ', M, N, X, 1, -1, -1, -1, 
     $               RDEST, CDEST)
      ELSE IF (OP.EQ.'MAX') THEN
        CALL IGAMX2D(ICTXT, SCOPE, ' ', M, N, X, 1, -1, -1, -1,
     $               RDEST, CDEST)
      ELSE
        CALL IGSUM2D(ICTXT, SCOPE, ' ', M, N, X, 1, RDEST, CDEST)
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


! Global-to-local pair of indices; shorthand for calling INDXG2L twice.
      SUBROUTINE G2LPAIR(I, J, GI, GJ, DESC, BLACS)
      IMPLICIT NONE
      ! IN/OUT
      INTEGER             I, J, GI, GJ, DESC(9), BLACS(2)
      ! Local
      INTEGER             DUM
      ! Functions
      INTEGER             INDXG2L
      
      
      I = INDXG2L(GI, DESC(5), DUM, DUM, BLACS(2))
      J = INDXG2L(GI, DESC(6), DUM, DUM, BLACS(3))
      
      RETURN
      END


! Construct local submatrix from global matrix
      SUBROUTINE MKSUBMAT(GBLX, SUBX, DESCX)!, RSRC, CSRC)
      IMPLICIT NONE
      ! IN/OUT
      INTEGER             DESCX(9), RSRC, CSRC
      DOUBLE PRECISION    GBLX(DESCX(3), DESCX(4)), SUBX(DESCX(9), *)
      ! Local
      INTEGER             M, N, I, J, GI, GJ, 
     $                    LDM(2), BLACS(5)
      ! External
      EXTERNAL            PDIMS, L2GPAIR
      
      
      ! Get local and proc grid info
      CALL PDIMS(DESCX, LDM, BLACS)
      
      M = LDM(1)
      N = LDM(2)
      
      ! Only do work if we own any local pieces
      IF (M.GT.0 .AND. N.GT.0) THEN
      
        DO J = 1, N
          DO I = 1, M
            CALL L2GPAIR(I, J, GI, GJ, DESCX, BLACS)
            SUBX(I,J) = GBLX(GI,GJ)
          END DO 
        END DO
      
      END IF
      
      RETURN
      END


! Construct global matrix from local submatrix
      SUBROUTINE MKGBLMAT(GBLX, SUBX, DESCX, RDEST, CDEST)
      IMPLICIT NONE
      ! IN/OUT
      INTEGER             DESCX(9), RDEST, CDEST, PROC
      DOUBLE PRECISION    GBLX(DESCX(3), DESCX(4)), SUBX(DESCX(9), *)
      ! Local
      INTEGER             M, N, I, J, GI, GJ, 
     $                    LDM(2), BLACS(5)
      ! Parameter
      DOUBLE PRECISION    ZERO
      PARAMETER ( ZERO = 0.0D0 )
      ! External
      EXTERNAL            PDIMS, G2LPAIR, DGSUM2D, DALLREDUCE
      
      
      ! Get local and proc grid info
      CALL PDIMS(DESCX, LDM, BLACS)
      
      M = LDM(1)
      N = LDM(2)
      
      GBLX = ZERO
      
      DO J = 1, N
        DO I = 1, M
          CALL L2GPAIR(I, J, GI, GJ, DESCX, BLACS)
          GBLX(GI,GJ) = SUBX(I,J)
        END DO 
      END DO
      
      IF (RDEST.EQ.-1) THEN
        CALL DALLREDUCE(DESCX(3), DESCX(4), GBLX, 'S', 'All', DESCX(2))
      ELSE
        CALL DREDUCE(DESCX(3), DESCX(4), GBLX, 'S', RDEST, CDEST,
     $               'All', DESCX(2))
      END IF
      
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
            RETURN! "Invalid argument 'DIAG'"
          END IF
        
        ! Zeroing the upper triangle
        ELSE IF (UPLO.EQ.'U') THEN
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
            RETURN! "Invalid argument 'DIAG'"
          END IF
        
        ! Zeroing the lower triangle
        ELSE IF (UPLO.EQ.'L') THEN
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
            RETURN! "Invalid argument 'DIAG'"
          END IF
        ELSE
          RETURN! "Invalid argument 'UPLO'"
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
          CALL PTRI2ZERO('L', 'N', X, DESCX)
        ! Zero Upper
        ELSE IF (UPLO.EQ.'L') THEN
          CALL PTRI2ZERO('U', 'N', X, DESCX)
        ! Otherwise...
        ELSE IF (UPLO.NE.'B') THEN
          RETURN! "Invalid argument 'UPLO'"
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
!!!!
!!!! If TRANS = 'Y', then A^T will be used.  This requires communication.
!!!!
!!!! For simplicity, we assume DESCA = DESCB.
!!!! Not as efficient as it could be, but them's the breaks.
!!!!
!!!! PDMKSYM should be more efficient if A=B.
!!!      SUBROUTINE PDGELACPY2(TRANS, UPLO, DIAG, A, IA, JA, DESCA, 
!!!     $                      B, IB, JB)
!!!      IMPLICIT NONE
!!!      ! IN/OUT
!!!      INTEGER             IA, JA, DESCA(9), IB, JB
!!!      DOUBLE PRECISION    A(DESCA(9), *), B(DESCA(9), * )
!!!      CHARACTER*1         TRANS, UPLO, DIAG
!!!      ! Local
!!!      INTEGER             M, N, I, J, GI, GJ
!!!      DOUBLE PRECISION, ALLOCATABLE CPX(:,:)
!!!      ! Parameter
!!!      DOUBLE PRECISION    ONE
!!!      PARAMETER ( ONE = 1.0D0 )
!!!      ! External
!!!      EXTERNAL            PDIMS, PDLACPY, PTRI2ZERO, PDGEADD
!!!      
!!!      
!!!      ! Get local and proc grid info
!!!      CALL PDIMS(DESCA, LDM, BLACS)
!!!      
!!!      M = LDM(1)
!!!      N = LDM(2)
!!!      
!!!      ! Upper part to be copied
!!!      IF (TRANS.EQ.'N') THEN
!!!        ! Copy over the desired triangle
!!!        CALL PDLACPY(UPLO, M, N, A, IA, JA, DESCA, 
!!!     $               B, IB, JB, DESCB)
!!!        
!!!        ! Copy over the diagonal if needed
!!!        IF (DIAG.EQ.'Y') THEN
!!!          DO J = 1, N
!!!            DO I = 1, M
!!!              CALL L2GPAIR(I, J, GI, GJ, DESCX, BLACS)
!!!              IF (GI.EQ.GJ) THEN
!!!                B(I,J) = A(I,J)
!!!              END IF
!!!            END DO 
!!!          END DO
!!!        END IF
!!!      
!!!      ELSE IF (TRANS.EQ.'Y') THEN
!!!        IF (UPLO.EQ.'U' .OR. UPLO.EQ.'L') THEN
!!!          ! Zero out to-be unused parts of B
!!!          CALL PTRI2ZERO(UPLO, DIAG, B, DESCB) 
!!!          
!!!          ! Copy over via PDGEADD
!!!          CALL PDGEADD('N', DESCX(3), DESCX(4), ONE, X, IX, JX, 
!!!     $                  DESCX, ONE, X, IX, JX, DESCX)
!!!        ELSE
!!!          
!!!        END IF
!!!      ELSE
!!!        RETURN "Invalid argument 'TRANS'"
!!!      END IF
!!!      
!!!      RETURN
!!!      END

