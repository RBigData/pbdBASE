! Copyright 2013, Schmidt


! "Optimal" process grid shape; tries to make the grid as close to square as 
! possible.  If a square grid is not possible, e.g. NPROCS=10, then the 
! larger number goes to NROWS, e.g. NROWS=5, NCOLS=2.
! INPUTS
  ! NPROCS = Number of processors to be used.
! OUTPUTS
  ! NROWS = "Optimal" number of process rows
  ! NCOLS = "Optimal" number of process columns
      SUBROUTINE OPTIMALGRID(NPROCS, NROWS, NCOLS)
      IMPLICIT NONE
      ! IN/OUT
      INTEGER             NPROCS, NROWS, NCOLS
      ! Local
      INTEGER             I, N
      
      
      N = INT(SQRT(REAL(NPROCS)))
      
      DO I = 0, N-1
        NCOLS = N - I
        NROWS = INT(MOD(NPROCS, NCOLS))
        IF (NROWS.EQ.0) EXIT
      END DO
      
      NROWS = NPROCS / NCOLS
      
      RETURN
      END


! Takes a ScaLAPACK descriptor array and from it determines full (1) local
! dimension information (calling NUMROC) and (2) BLACS grid information.
! The local dimension LDIM is set to (/ 0, 0 /) if there is not actual 
! ownership of local data on that process.
! INPUTS
  ! DESC = ScaLAPACK descriptor array
! OUTPUTS
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
      
      IF (LDM(1).LT.1 .OR. LDM(2).LT.1) THEN
        LDM(1) = 0
        LDM(2) = 0
      END IF
      
      RETURN
      END


! Reductions using BLACS routines.  Bindings meant to mimic ALLREDUCE, which is
! much more friendly than the weird BLACS calls.
! INPUTS/OUTPUTS
  ! X = Array of data to be reduced
! INPUTS
  ! DESCX = Descriptor array of X
  ! OP = Operation: 'MIN', 'MAX', 'SUM'.  Default if you pass something else is 
    ! to sum.
  ! SCOPE = BLACS processor grid scope along which the reduction will take 
    ! place: 'All', 'Column', 'Row'
  ! RDEST/CDEST: Row/Column destination in the case of *REDUCE (rather than ALL)
      SUBROUTINE DALLREDUCE(X, DESCX, OP, SCOPE)
      IMPLICIT NONE
      ! IN/OUT
      INTEGER             DESCX(9)
      DOUBLE PRECISION    X( * )
      CHARACTER           OP, SCOPE
      ! Local
      INTEGER             M, N, LDA, ICTXT
      ! External
      EXTERNAL            DGSUM2D, DGAMX2D, DGMN2D
      
      
      M = DESCX(3)
      N = DESCX(4)
      LDA = DESCX(9)
      ICTXT = DESCX(2)
      
      IF (OP.EQ.'MIN') THEN
        CALL DGAMN2D(ICTXT,SCOPE,' ',M,N,X,LDA,-1,-1,-1,-1,-1)
      ELSE IF (OP.EQ.'MAX') THEN
        CALL DGAMX2D(ICTXT,SCOPE,' ',M,N,X,LDA,-1,-1,-1,-1,-1)
      ELSE ! default to sum
        CALL DGSUM2D(ICTXT, SCOPE, ' ', M, N, X, LDA, -1, -1)
      END IF
      
      RETURN 
      END


      SUBROUTINE DREDUCE(X, DESCX, OP, RDEST, CDEST, SCOPE)
      IMPLICIT NONE
      ! IN/OUT
      INTEGER             DESCX(9), RDEST, CDEST
      DOUBLE PRECISION    X( * )
      CHARACTER           OP, SCOPE
      ! Local
      INTEGER             M, N, LDA, ICTXT
      ! External
      EXTERNAL            DGSUM2D, DGAMX2D, DGMN2D
      
      
      M = DESCX(3)
      N = DESCX(4)
      LDA = DESCX(9)
      ICTXT = DESCX(2)
      
      IF (OP.EQ.'MIN') THEN
        CALL DGAMN2D(ICTXT, SCOPE, ' ', M, N, X, LDA, -1, -1, -1, 
     $               RDEST, CDEST)
      ELSE IF (OP.EQ.'MAX') THEN
        CALL DGAMX2D(ICTXT, SCOPE, ' ', M, N, X, LDA, -1, -1, -1,
     $               RDEST, CDEST)
      ELSE
        CALL DGSUM2D(ICTXT, SCOPE, ' ', M, N, X, LDA, RDEST, CDEST)
      END IF
      
      RETURN 
      END


      SUBROUTINE IALLREDUCE(X, DESCX, OP, SCOPE)
      IMPLICIT NONE
      ! IN/OUT
      INTEGER             DESCX(9), X( * )
      CHARACTER           OP, SCOPE
      ! Local
      INTEGER             M, N, LDA, ICTXT
      ! External
      EXTERNAL            IGSUM2D, IGAMX2D, IGMN2D
      
      
      M = DESCX(3)
      N = DESCX(4)
      LDA = DESCX(9)
      ICTXT = DESCX(2)
      
      IF (OP.EQ.'MIN') THEN
        CALL IGAMN2D(ICTXT, SCOPE, ' ', M, N, X, 1,-1,-1,-1,-1,-1)
      ELSE IF (OP.EQ.'MAX') THEN
        CALL IGAMX2D(ICTXT, SCOPE, ' ', M, N, X, 1,-1,-1,-1,-1,-1)
      ELSE
        CALL IGSUM2D(ICTXT, SCOPE, ' ', M, N, X, 1,-1,-1)
      END IF
      
      RETURN 
      END


      SUBROUTINE IREDUCE(X, DESCX, OP, RDEST, CDEST, SCOPE)
      IMPLICIT NONE
      ! IN/OUT
      INTEGER             DESCX(9), X( * ), RDEST, CDEST
      CHARACTER           OP, SCOPE
      ! Local
      INTEGER             M, N, LDA, ICTXT
      ! External
      EXTERNAL            IGSUM2D, IGAMX2D, IGMN2D
      
      
      M = DESCX(3)
      N = DESCX(4)
      LDA = DESCX(9)
      ICTXT = DESCX(2)
      
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
! INPUTS
  ! I/J = Local coordinates.
  ! DESC = BLACS descriptor array.
  ! BLACS = BLACS process grid information, taken from PDIMS.
! OUTPUTS
  ! GI/GJ = Global coordinates.
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
! INPUTS
  ! GI/GJ = Global coordinates.
  ! DESC = BLACS descriptor array.
  ! BLACS = BLACS process grid information, taken from PDIMS.
! OUTPUTS
  ! I/J = Local coordinates.
      SUBROUTINE G2LPAIR(I, J, GI, GJ, DESC, BLACS)
      IMPLICIT NONE
      ! IN/OUT
      INTEGER             I, J, GI, GJ, DESC(9), BLACS(2)
      ! Local
      INTEGER             DUM
      ! Functions
      INTEGER             INDXG2L
      
      
      I = INDXG2L(GI, DESC(5), DUM, DUM, BLACS(2))
      J = INDXG2L(GJ, DESC(6), DUM, DUM, BLACS(3))
      
      RETURN
      END


! Wrapper for pdgemr2d
! INPUTS
  ! X = Input submatrix.
  ! IX/JX = 
  ! DESCX = Descriptor array for X.
  ! IY/JY = 
  ! DESCY = Descriptor array for Y.
  ! CMNCTXT = Common BLACS context for X and Y.
! OUTPUTS
  ! Y = 
!!!      SUBROUTINE REDIST(X, IX, JX, DESCX, Y, IY, JY, DESCY, CMNCTXT)
!!!      IMPLICIT NONE
!!!      ! IN/OUT
!!!      INTEGER             IX, JX, DESCX(9), IY, JY, DESCY(9), CMNCTXT
!!!      DOUBLE PRECISION    X( * ), Y( * )
!!!      ! Local
!!!      INTEGER             M, N, MXLDM, DESCA(9),
!!!     $                    LDMX(2), LDMY(2), BLACSX(4), BLACSY(4)
!!!      ! External
!!!      EXTERNAL            PDGEMR2D
!!!      
!!!      
!!!      ! Get local and proc grid info
!!!      CALL PDIMS(DESCX, LDMX, BLACSX)
!!!      CALL PDIMS(DESCX, LDMY, BLACSY)
!!!      
!!!      M = DESCX(3)
!!!      N = DESCX(4)
!!!      
!!!      ! Adjust LDA since PDGEMR2D crashes all the time when LDA=1
!!!      DESCA(3) = 1
!!!      DESCA(4) = 1
!!!      DESCA(9) = 1
!!!      
!!!      MXLDM = MAX(LDMX)
!!!      DESCA(2) = DESCX(2)
!!!      CALL IALLREDUCE(MXLDM, DESCA, 'MAX', 'All')
!!!      IF (DESCX(9).EQ.1 .AND. DESCX(3).GT.1) DESCX(9) = MXLDM
!!!      
!!!      MXLDM = MAX(LDMY)
!!!      DESCA(2) = DESCY(2)
!!!      CALL IALLREDUCE(MXLDM, DESCA, 'MAX', 'All')
!!!      IF (DESCY(9).EQ.1 .AND. DESCY(3).GT.1) DESCY(9) = MXLDM
!!!      
!!!      ! Redistribute
!!!      CALL PDGEMR2D(M, N, X, IX, JX, DESCX,
!!!     $              Y, IY, JY, DESCY, CMNCTXT)
!!!      
!!!      RETURN
!!!      END


! Construct local submatrix from global matrix
! INPUTS
  ! GBLX = Global, non-distributed matrix.  Owned by which processor(s) depends 
    ! on R/CSRC values
  ! DESCX = ScaLAPACK descriptor array for SUBX (not a typo).
  ! RSRC/CSRC = Row/Column process value corresponding to BLACS grid for the
    ! value in DESCX(2) (the ICTXT) on which the data is stored.  If RSRC = -1,
    ! then CSRC is ignored and total ownership is assumed, i.e., GBLX is owned 
    ! by all processors.
! OUTPUTS
  ! SUBX = Local submatrix.
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


! Construct global matrix from local submatrix.
! INPUTS
  ! SUBX = Local submatrix.
  ! DESCX = ScaLAPACK descriptor array for SUBX.
  ! RDEST/CDEST = Row/Column process value corresponding to BLACS grid for the
    ! value in DESCX(2) (the ICTXT) on which the global matrix GBLX will be 
    ! stored.  If RDEST = -1, then CDEST is ignored and total ownership is 
    ! assumed, i.e., GBLX is given to all processors.
! OUTPUTS
  ! GBLX = Global, non-distributed matrix.
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
      EXTERNAL            PDIMS, L2GPAIR, DGSUM2D, DALLREDUCE
      
      
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
        CALL DALLREDUCE(GBLX, DESCX, 'S', 'All')
      ELSE
        CALL DREDUCE(GBLX, DESCX, 'S', RDEST, CDEST, 'All')
      END IF
      
      RETURN
      END


! Zero out the (global) triangle and/or diagonal of a distributed matrix
! INPUTS/OUTPUTS
  ! X = Array of data to be zeroed (in part or whole)
! INPUTS
  ! UPLO = Char specifying whether 'U'pper or 'L'ower triangle is to be zeroed,
    ! or whether 'B'oth should be zeroed, or whether 'N'either triangle should.
  ! DIAG = Char specifying whether 'Y'es, zero the diagonal too, or 'N'o, do not.
! EXAMPLE:  UPLO='B', DIAG='N' produces a diagonal matrix
!           UPLO='N', DIAG='Y' zeroes out only the diagonal
!           UPLO='L', DIAG='Y' zeros the diagonal and upper triangle
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
      EXTERNAL            PDIMS, L2GPAIR, PDLASET
      
      
      ! Get local and proc grid info
      CALL PDIMS(DESCX, LDM, BLACS)
      
      M = LDM(1)
      N = LDM(2)
      
      
      ! Let ScaLAPACK do the work if possible:
      IF ( (UPLO.EQ.'U' .OR. UPLO.EQ.'L') .AND. DIAG.EQ.'Y') THEN
          CALL PDLASET(UPLO, DESCX(3), DESCX(4), ZERO, ZERO, X, 
     $                 1, 1, DESCX)
      ! Only do work if we own any local pieces
      ELSE IF (M.GT.0 .AND. N.GT.0) THEN
        
        ! Diagonal only (zeroing neither triangle
        IF (UPLO.EQ.'N') THEN
          ! Quick return if possible
          IF (DIAG.EQ.'N') THEN
            RETURN
          ELSE IF (DIAG.EQ.'Y') THEN
            DO J = 1, N
              DO I = 1, M
                CALL L2GPAIR(I, J, GI, GJ, DESCX, BLACS)
                IF (GI.EQ.GJ) THEN
                  X(I,J) = ZERO
                END IF
              END DO 
            END DO
          ELSE
            RETURN! "Invalid argument 'DIAG'"
          END IF
        END IF
        
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
!          ! Zero the diagonal as well
!          IF (DIAG.EQ.'Y') THEN 
!            DO J = 1, N
!              DO I = 1, M
!                CALL L2GPAIR(I, J, GI, GJ, DESCX, BLACS)
!                IF (GI.LE.GJ) THEN
!                  X(I,J) = ZERO
!                END IF
!              END DO
!            END DO
          ! Leave the diagonal alone
          IF (DIAG.EQ.'N') THEN 
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
!          ! Zero the diagonal as well
!          IF (DIAG.EQ.'Y') THEN 
!            DO J = 1, N
!              DO I = 1, M
!                CALL L2GPAIR(I, J, GI, GJ, DESCX, BLACS)
!                IF (GI.GE.GJ) THEN
!                  X(I,J) = ZERO
!                END IF
!              END DO
!            END DO
          ! Leave the diagonal alone
          IF (DIAG.EQ.'N') THEN 
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


! Make the matrix symmetric via copying from one (global) triangle to the other.
! INPUTS/OUTPUTS
  ! X = Array of data to be zeroed (in part or whole)
! INPUTS
  ! UPLO = Copy FROM 'U'pper or copy FROM 'L'ower.  The modified data lives in 
    ! the opposite triangle from UPLO.
  ! IX/JX = 
  ! DESCX = ScaLAPACK descriptor array for X.
      SUBROUTINE PDMKSYM(UPLO, X, IX, JX, DESCX)
      IMPLICIT NONE
      ! IN/OUT
      INTEGER             IX, JX, DESCX(9)
      DOUBLE PRECISION    X(DESCX(9), *)
      CHARACTER*1         UPLO
      ! Local
      INTEGER             M, N, I, J, GI, GJ, LDM(2), BLACS(5)
      DOUBLE PRECISION, ALLOCATABLE :: CPX(:,:)
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
        
        ALLOCATE(CPX(M,N))
        CPX(1:M,1:N) = X(1:M,1:N)
        
        ! X = t(X) + X
        CALL PDGEADD('T', DESCX(3), DESCX(4), ONE, CPX, IX, JX, 
     $                DESCX, ONE, X, IX, JX, DESCX)
        
        DEALLOCATE(CPX)
        
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


!!! Does what pdlacpy SHOULD do, namely copies the triangle of a global 
!!! distributed matrix onto the triangle of another matrix.
!!!
!!! If TRANS = 'Y', then A^T will be used.  This requires communication.
!!!
!!! For simplicity, we assume DESCA = DESCB.
!!! Not as efficient as it could be, but them's the breaks.
!!!
!!! PDMKSYM should be more efficient if A=B.
!!      SUBROUTINE PDGELACPY2(TRANS, UPLO, DIAG, A, IA, JA, DESCA, 
!!     $                      B, IB, JB, INFO)
!!      IMPLICIT NONE
!!      ! IN/OUT
!!      INTEGER             IA, JA, DESCA(9), IB, JB, INFO
!!      DOUBLE PRECISION    A(DESCA(9), *), B(DESCA(9), * )
!!      CHARACTER*1         TRANS, UPLO, DIAG
!!      ! Local
!!      INTEGER             M, N, I, J, GI, GJ, LDIM(2), BLACS(5)
!!      DOUBLE PRECISION, ALLOCATABLE :: CPX(:,:)
!!      ! Parameter
!!      DOUBLE PRECISION    ZERO, ONE
!!      PARAMETER ( ZERO = 0.0D0, ONE = 1.0D0 )
!!      ! External
!!      EXTERNAL            PDIMS, PDLACPY, PTRI2ZERO, PDGEADD
!!      
!!      
!!      INFO = 0
!!      
!!      ! Get local and proc grid info
!!      CALL PDIMS(DESCA, LDIM, BLACS)
!!      
!!      M = LDIM(1)
!!      N = LDIM(2)
!!      
!!      ! Upper part to be copied
!!      IF (TRANS.EQ.'N') THEN
!!        IF (UPLO.EQ.'U') THEN
!!          ! Copy over the diagonal if needed
!!          IF (DIAG.EQ.'Y') THEN
!!            DO J = 1, N
!!              DO I = 1, M
!!                CALL L2GPAIR(I, J, GI, GJ, DESCA, BLACS)
!!                IF (GI.GE.GJ) THEN
!!                  B(I,J) = A(I,J)
!!                END IF
!!              END DO 
!!            END DO
!!          ELSE IF (DIAG.EQ.'N') THEN
!!            DO J = 1, N
!!              DO I = 1, M
!!                CALL L2GPAIR(I, J, GI, GJ, DESCA, BLACS)
!!                IF (GI.GT.GJ) THEN
!!                  B(I,J) = A(I,J)
!!                END IF
!!              END DO 
!!            END DO
!!          ELSE
!!            INFO = -3 ! Invalid argument 'DIAG'
!!            RETURN
!!          END IF
!!      
!!      ELSE IF (TRANS.EQ.'Y') THEN
!!        IF (UPLO.EQ.'U' .OR. UPLO.EQ.'L') THEN
!!          ! Zero out to-be unused parts of B
!!          CALL PTRI2ZERO(UPLO, DIAG, B, DESCA) 
!!          
!!          ! Copy over via PDGEADD
!!          CALL PDGEADD('N', DESCA(3), DESCA(4), ONE, A, IA, JA, 
!!     $                  DESCA, ONE, A, IA, JA, DESCA)
!!        ELSE
!!          
!!        END IF
!!      ELSE
!!        INFO = -1 ! Invalid argument 'TRANS'
!!        RETURN 
!!      END IF
!!      
!!      RETURN
!!      END




! Copyright 2013, Schmidt

! For internal use.
! In the case of Matrix-Vector operations where the vector is global and not
! necessarily of "appropriate" length, this is used to adjust which element
! of the vector is used with the matrix.  Enables the equivalent of doing,
! for example, something like matrix(1, nrow=5, ncol=2) + 1:3 in R.
! INPUTS
  ! I = Index
  ! M = Modulus
      INTEGER FUNCTION IND(I, M)
      IMPLICIT NONE
      INTEGER             I, M
      
      
      IND = MOD(I, M)
      IF (IND.EQ.0) THEN
        IND = M
      END IF
      
      RETURN
      END


! Convert matrix indexing to vector indexing, or vice versa
!      SUBROUTINE INDMAT2VEC(STORAGE, M, N, I, J, K)
!      IMPLICIT NONE
!      CHARACTER*1         STORAGE
!      INTEGER             M, N, I, J, K
!      ! Function
!      INTEGER             IND
!      
!      
!      IF (STORAGE .EQ. 'C') THEN
!        
!      ELSE IF (STORAGE .EQ. 'R') THEN
!        
!      ELSE
!        IND = -1
!      END IF
!      
!      RETURN
!      END


!      SUBROUTINE INDVEC2MAT(STORAGE, M, N, I, J, K)
!      IMPLICIT NONE
!      CHARACTER*1         STORAGE
!      INTEGER             M, N, I, J, K
!      ! Function
!      INTEGER             IND
!      
!      
!      IF (STORAGE .EQ. 'C') THEN
!        I = IND(K, M)
!        J = K/M + 1
!      ELSE IF (STORAGE .EQ. 'R') THEN
!        I = IND(K, N)
!        J = K/N + 1
!      ELSE
!        I = -1
!        J = -1
!      END IF
!      
!      RETURN
!      END



! SWEEP array out of distributed matrix
! INPUTS/OUTPUTS
  ! X = Submatrix of data which should globally be "swept"
! INPUTS
  ! IX/JX = 
  ! DESCX = Descriptor array for X
  ! VEC = Vector to "sweep" through X
  ! LVEC = Length of VEC
  ! MARGIN = 1 for row sweeping, 2 for column sweeping
  ! FUN = Char with 4 possibilities, describing the type of sweep to perform:
    ! "+", "-", "*", "/"
      SUBROUTINE PDSWEEP(X, IX, JX, DESCX, VEC, LVEC, MARGIN, FUN)
      IMPLICIT NONE
      ! IN/OUT
      INTEGER             IX, JX, DESCX(9), MARGIN, LVEC
      DOUBLE PRECISION    X(DESCX(9), *), VEC(LVEC)
      CHARACTER*1         FUN
      ! Local
      INTEGER             K, M, N, POS, I, J, GI, GJ, LDM(2), BLACS(5)
      ! External
      EXTERNAL            PDIMS
      ! Function
      INTEGER             IND
      
      
      ! Get local and proc grid info
      CALL PDIMS(DESCX, LDM, BLACS)
      
      M = LDM(1)
      N = LDM(2)
      
      ! Only do work if we own any local pieces
      IF (M.GT.0 .AND. N.GT.0) THEN
        ! Addition
        IF (FUN.EQ."+") THEN
          IF (MARGIN.EQ.1) THEN
            K = DESCX(3)
            DO J = 1, N
              DO I = 1, M
                CALL L2GPAIR(I, J, GI, GJ, DESCX, BLACS)
                POS = IND(GI + K*(GJ-1), LVEC)
                X(I, J) = X(I, J) + VEC(POS)
              END DO
            END DO
          ELSE IF (MARGIN.EQ.2) THEN
            K = DESCX(4)
            DO J = 1, N
              DO I = 1, M
                CALL L2GPAIR(I, J, GI, GJ, DESCX, BLACS)
                POS = IND(GJ + K*(GI-1), LVEC)
                X(I, J) = X(I, J) + VEC(POS)
              END DO
            END DO
          END IF
        ! Subtraction
        ELSE IF (FUN.EQ."-") THEN
          IF (MARGIN.EQ.1) THEN
            K = DESCX(3)
            DO J = 1, N
              DO I = 1, M
                CALL L2GPAIR(I, J, GI, GJ, DESCX, BLACS)
                POS = IND(GI + K*(GJ-1), LVEC)
                X(I, J) = X(I, J) - VEC(POS)
              END DO
            END DO
          ELSE IF (MARGIN.EQ.2) THEN
            K = DESCX(4)
            DO J = 1, N
              DO I = 1, M
                CALL L2GPAIR(I, J, GI, GJ, DESCX, BLACS)
                POS = IND(GJ + K*(GI-1), LVEC)
                X(I, J) = X(I, J) - VEC(POS)
              END DO
            END DO
          END IF
        ! Multiplication
        ELSE IF (FUN.EQ."*") THEN
          IF (MARGIN.EQ.1) THEN
            K = DESCX(3)
            DO J = 1, N
              DO I = 1, M
                CALL L2GPAIR(I, J, GI, GJ, DESCX, BLACS)
                POS = IND(GI + K*(GJ-1), LVEC)
                X(I, J) = X(I, J) * VEC(POS)
              END DO
            END DO
          ELSE IF (MARGIN.EQ.2) THEN
            K = DESCX(4)
            DO J = 1, N
              DO I = 1, M
                CALL L2GPAIR(I, J, GI, GJ, DESCX, BLACS)
                POS = IND(GJ + K*(GI-1), LVEC)
                X(I, J) = X(I, J) * VEC(POS)
              END DO
            END DO
          END IF
        ! Division
        ELSE IF (FUN.EQ."/") THEN
          IF (MARGIN.EQ.1) THEN
            K = DESCX(3)
            DO J = 1, N
              DO I = 1, M
                CALL L2GPAIR(I, J, GI, GJ, DESCX, BLACS)
                POS = IND(GI + K*(GJ-1), LVEC)
                X(I, J) = X(I, J) / VEC(POS)
              END DO
            END DO
          ELSE IF (MARGIN.EQ.2) THEN
            K = DESCX(4)
            DO J = 1, N
              DO I = 1, M
                CALL L2GPAIR(I, J, GI, GJ, DESCX, BLACS)
                POS = IND(GJ + K*(GI-1), LVEC)
                X(I, J) = X(I, J) / VEC(POS)
              END DO
            END DO
          END IF
        END IF
      END IF
      
      RETURN
      END


! Grab diagonal of matrix
! INPUTS
  ! X = Submatrix whose global diagonal is to be taken
  ! IX/JX = 
  ! DESCX = Descriptor array for X
! OUTPUTS
  ! DIAG = Diagonal of the global matrix for which X is the submatrix.
      SUBROUTINE PDGDGTK(X, IX, JX, DESCX, DIAG, RDEST, CDEST)
      IMPLICIT NONE
      ! IN/OUT
      INTEGER             IX, JX, DESCX(9), RDEST, CDEST
      DOUBLE PRECISION    X(DESCX(9), *), DIAG( * )
      CHARACTER*1         REDUCE
      ! Local
      INTEGER             K, M, N, I, J, GI, GJ, LDM(2), 
     $                    BLACS(5), DESC(9)
      ! Parameter
      DOUBLE PRECISION    ZERO
      PARAMETER ( ZERO = 0.0D0 )
      ! External
      EXTERNAL            PDIMS, DALLREDUCE
      
      
      ! Get local and proc grid info
      CALL PDIMS(DESCX, LDM, BLACS)
      
      M = LDM(1)
      N = LDM(2)
      
      K = MIN(DESCX(3), DESCX(4))
      
      DIAG(1:K) = ZERO
      
      DO J = 1, N
        DO I = 1, M
          CALL L2GPAIR(I, J, GI, GJ, DESCX, BLACS)
          IF (GI.EQ.GJ) THEN
            DIAG(GI) = X(I,J)
          END IF
        END DO
      END DO
      
      CALL DGSUM2D(DESCX(2), 'All', ' ', K, 1, DIAG, K, RDEST, CDEST)
      
      RETURN
      END


! Construct matrix containing diagonal DIAG
! INPUTS
  ! IX/JX = 
  ! DESCX = Descriptor array for X
  ! DIAG = Global vector which should be the diagonal of a distributed matrix
  ! LDIAG = Length of Diag.  Need not be equal to K = MIN(DESCX(3), DESCX(4)),
    ! though only the first K elements of DIAG will be used.
! OUTPUTS
  ! X = Submatrix of the global matrix with diagonal DIAG
      SUBROUTINE PDDIAGMK(X, IX, JX, DESCX, DIAG, LDIAG)
      IMPLICIT NONE
      ! IN/OUT
      INTEGER             IX, JX, DESCX(9), LDIAG
      DOUBLE PRECISION    X(DESCX(9), *), DIAG( * )
      ! Local
      INTEGER             M, N, I, J, GI, GJ, LDM(2), BLACS(5)
      ! Parameter
      DOUBLE PRECISION    ZERO
      PARAMETER ( ZERO = 0.0D0 )
      ! External
      EXTERNAL            PDIMS
      ! Function
      INTEGER             IND
      
      
      ! Get local and proc grid info
      CALL PDIMS(DESCX, LDM, BLACS)
      
      M = LDM(1)
      N = LDM(2)
      
      DO J = 1, N
        DO I = 1, M
          CALL L2GPAIR(I, J, GI, GJ, DESCX, BLACS)
          IF (GI.EQ.GJ) THEN
            X(I,J) = DIAG( IND(GI, LDIAG) )
          ELSE 
            X(I,J) = ZERO
          END IF
        END DO
      END DO
      
      RETURN
      END



