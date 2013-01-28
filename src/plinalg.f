! Copyright 2013, Schmidt

!!!! Distance D_METHOD(X,Y)
!!!! Methods:  'E' - Euclidean
!!!!           'S' - Supremum (infinity)
!!!!           'I' - Infimum (-infinity)
!!!!           'H' - Manhattan
!!!!           'C' - Canberra
!!!!           'K' - Minkowski
!!!      DOUBLE PRECISION FUNCTION DDIST(METHOD, N, X, INCX, Y, INCY, P)
!!!      IMPLICIT NONE
!!!      ! IN/OUT
!!!      INTEGER             N, P, INCX, INCY
!!!      DOUBLE PRECISION    X( * ), Y( * )
!!!      CHARACTER*1         METHOD
!!!      ! Local
!!!      INTEGER             I
!!!      DOUBLE PRECISION, ALLOCATABLE :: WORK(:)
!!!      ! Parameter
!!!      DOUBLE PRECISION    ONE
!!!      PARAMETER ( ONE = 1.0D0 )
!!!      ! Functions
!!!      EXTERNAL            DCOPY, DAXPY
!!!      INTEGER             IDAMAX, IDAMIN
!!!      DOUBLE PRECISION    DASUM, DNRM2, DNRM3
!!!      
!!!      ALLOCATE(WORK(N))
!!!      
!!!      ! ! WORK = X - Y
!!!      CALL DCOPY(N, X, INCX, WORK, INCY)
!!!      
!!!      CALL DAXPY(N, ONE, WORK, 1, Y, INCY)
!!!      
!!!      IF (METHOD.EQ.'E') THEN
!!!        DDIST = DNRM2(N, WORK, 1)
!!!      ELSE IF (METHOD.EQ.'S') THEN
!!!        DDIST = X( IDAMAX(N, X, INCX) )
!!!      ELSE IF (METHOD.EQ.'I') THEN
!!!        DDIST = X( IDAMIN(N, X, INCX) )
!!!      ELSE IF (METHOD.EQ.'H') THEN
!!!        DDIST = DASUM(N, WORK, 1)
!!!      ELSE IF (METHOD.EQ.'C') THEN
!!!        DDIST = 0
!!!        DO I = 1, N
!!!          DDIST = ABS(WORK(I)) / 
!!!     $       (ABS(X((I-1)*INCX + 1)) + ABS(Y((I-1)*INCY + 1)))
!!!        END DO
!!!      ELSE IF (METHOD.EQ.'K') THEN
!!!        DDIST = DNRM3(N, WORK, 1, P)
!!!      END IF
!!!      DEALLOCATE(WORK)
!!!      
!!!      RETURN
!!!      END


!!!! Mahalanobis distance d_M(X,Y)
!!!      DOUBLE PRECISION FUNCTION DMAHAL(N, X, INCX, Y, INCY, MN, COVINV)
!!!      IMPLICIT NONE
!!!      ! IN/OUT
!!!      INTEGER             N, INCX, INCY
!!!      DOUBLE PRECISION    X( * ), Y( * ), MN( * ), COVINV( * )
!!!      ! Local
!!!      INTEGER             I
!!!      DOUBLE PRECISION, ALLOCATABLE :: WORK1(:)
!!!      DOUBLE PRECISION, ALLOCATABLE :: WORK2(:)
!!!      ! Parameter
!!!      DOUBLE PRECISION    ZERO, ONE
!!!      PARAMETER ( ZERO = 0.0D0, ONE = 1.0D0 )
!!!      ! External
!!!      EXTERNAL            DCOPY, DAXPY, DGEMV
!!!      DOUBLE PRECISION    DDOT
!!!      
!!!      ALLOCATE(WORK1(N))
!!!      ALLOCATE(WORK2(N))
!!!      
!!!      ! WORK1 = X - Y
!!!      CALL DCOPY(N, X, INCX, WORK1, INCY)
!!!      
!!!      CALL DAXPY(N, ONE, WORK1, 1, Y, INCY)
!!!      
!!!      ! DMAHAL = WORK1^T * COVINV * WORK2
!!!      CALL DGEMV('N', N, N, ONE, COVINV, N, WORK1, 1, ZERO, WORK2, 1)
!!!      
!!!      DMAHAL = SQRT(DDOT(N, WORK1, 1, WORK2, 1))
!!!      
!!!      RETURN
!!!      END


! X^T * X or X * X^T
! TRANS = 'T' :  X^T*X, TRANS = 'N' : X*X^T
      SUBROUTINE PDCROSSPROD(TRANS, ALPHA, X, IX, JX,
     $                       DESCX, C, IC, JC, DESCC)
      IMPLICIT NONE
      ! IN/OUT
      INTEGER             IX, JX, DESCX(9), IC, JC, DESCC(9)
      DOUBLE PRECISION    X( * ), C( * ), ALPHA
      CHARACTER*1         TRANS
      ! Local
      INTEGER             LDX, LDC
      CHARACTER*1         NST
      ! Parameter
      DOUBLE PRECISION    ZERO
      PARAMETER ( ZERO = 0.0D0 )
      ! External
      EXTERNAL            DSYRK, DMKSYM
      
      
      IF (TRANS.EQ.'T') THEN
        NST = 'N'
        LDX = DESCX(4)
        LDC = DESCX(3)
      ELSE 
        NST = 'T'
        LDX = DESCX(3)
        LDC = DESCX(4)
      END IF
      
      ! Compute upper triangle of X^T*X or X^T*X
      CALL PDSYRK('U', NST, LDC, LDX, ALPHA, X, IX, JX, DESCX,
     $            ZERO, C, IC, JC, DESCC)
      
      ! Fill lower triangle (make symmetric)
      CALL PDMKSYM('U', C, IC, JC, DESCC)
      
      RETURN
      END 


!!!! Matrix norm; wrapper for DLANGE.  Handles allocation, etc.
!!!      SUBROUTINE DMATNORM(VALUE, NORM, M, N, X)
!!!      IMPLICIT NONE
!!!      ! In/Out
!!!      CHARACTER*1         NORM
!!!      INTEGER             M, N
!!!      DOUBLE PRECISION    VALUE, X( * )
!!!      ! Local
!!!      INTEGER             LWORK
!!!      DOUBLE PRECISION    LDX
!!!      DOUBLE PRECISION, ALLOCATABLE :: WORK(:)
!!!      ! Functions
!!!      DOUBLE PRECISION    DLANGE
!!!      
!!!      IF (NORM.EQ."I") THEN
!!!        LWORK = M
!!!      ELSE 
!!!        LWORK = 0
!!!      END IF
!!!      
!!!      LDX = MAX(1, M)
!!!      
!!!      ALLOCATE (WORK(LWORK))
!!!      
!!!      VALUE = DLANGE(NORM, M, N, X, LDX, WORK)
!!!      
!!!      DEALLOCATE (WORK)
!!!      
!!!      RETURN
!!!      END


!!!!     Condition number estimator for general matrix
!!!!       Step 1:  get matrix norm of X
!!!!       Step 2:  factor X=L
!!!!       Step 3:  Call dgecon
!!!      SUBROUTINE DCONDNUM(NORM, M, N, X, RCOND, INFO)
!!!      IMPLICIT NONE
!!!      ! #######################################################################
!!!      ! In/Out
!!!      CHARACTER           NORM
!!!      INTEGER             M, N, LDIM1, INFO
!!!      DOUBLE PRECISION    RCOND, X( * )
!!!      ! Local
!!!      INTEGER             LWORK, LIPIV
!!!      DOUBLE PRECISION    XNORM, TMP, LDX
!!!      ! Parameters
!!!      INTEGER             IN1
!!!      PARAMETER ( IN1 = -1 )
!!!      ! Dynamic
!!!      DOUBLE PRECISION, ALLOCATABLE :: WORK(:)
!!!      INTEGER, ALLOCATABLE :: IWORK(:), IPIV(:)
!!!      ! Functions
!!!      EXTERNAL            MATNORM, DGETRF, DGECON
!!!      ! #######################################################################
!!!      
!!!      ! Step 1:  get matrix norm of A
!!!      CALL DMATNORM(XNORM, NORM, N, N, X)
!!!      
!!!      ! Step 2:  factor A=LU
!!!      LIPIV = MAX(1, MIN(M, N))
!!!      ALLOCATE (IPIV(LIPIV))
!!!      
!!!      LDX = MAX(1, M)
!!!      
!!!      CALL DGETRF(M, N, X, LDX, INFO)
!!!      
!!!      IF (INFO.NE.0) THEN
!!!        RETURN
!!!      END IF
!!!      
!!!      LDX = MAX(1, N)
!!!      
!!!      ! Step 3:  Call pdgecon
!!!      CALL DGECON(NORM, N, X, LDX, XNORM, RCOND, TMP,
!!!     $            N, INFO)
!!!      
!!!      LWORK = INT(TMP)
!!!      ALLOCATE (WORK(LWORK))
!!!      ALLOCATE (IWORK(N))
!!!      
!!!      CALL DGECON(NORM, N, X, LDX, XNORM, RCOND, WORK,
!!!     $             IWORK, INFO)
!!!      
!!!      DEALLOCATE (IPIV)
!!!      DEALLOCATE (IWORK)
!!!      DEALLOCATE (WORK)
!!!      
!!!      RETURN
!!!      END


!!!! Compute matrix inverse without having to understand LAPACK/BLAS peculiarities
!!!! In place version (X is overwritten with X^-1)
!!!      SUBROUTINE DINVIP(N, X, INFO)
!!!      IMPLICIT NONE
!!!      ! IN/OUT
!!!      INTEGER             N
!!!      DOUBLE PRECISION    X(N, N)
!!!      ! Local
!!!      INTEGER             I, J, LWORK, INFO
!!!      DOUBLE PRECISION    TMP
!!!      INTEGER, ALLOCATABLE :: IPIV(:)
!!!      DOUBLE PRECISION, ALLOCATABLE :: WORK(:)
!!!      ! External
!!!      EXTERNAL           DGETRF, DGETRI
!!!      
!!!      ! Factor X=LU
!!!      ALLOCATE(IPIV(N))
!!!      CALL DGETRF(N, N, X, N, IPIV, INFO)
!!!      
!!!      IF (INFO.NE.0) RETURN
!!!      
!!!      ! Invert X
!!!      LWORK = -1
!!!      
!!!      CALL DGETRI(N, X, N, IPIV, TMP, LWORK, INFO)
!!!      
!!!      LWORK = INT(TMP)
!!!      ALLOCATE(WORK(LWORK))
!!!      
!!!      CALL DGETRI(N, X, N, IPIV, WORK, LWORK, INFO)
!!!      
!!!      DEALLOCATE(WORK)
!!!      DEALLOCATE(IPIV)
!!!      
!!!      RETURN
!!!      END


!!!! Compute matrix inverse without having to understand LAPACK/BLAS peculiarities
!!!! Non-in-place version (on return, INV = X^-1)
!!!      SUBROUTINE DINV(N, X, INV, INFO)
!!!      IMPLICIT NONE
!!!      ! IN/OUT
!!!      INTEGER             N
!!!      DOUBLE PRECISION    X(N, N), INV(N, N)
!!!      ! Local
!!!      INTEGER             I, J, LWORK, INFO
!!!      DOUBLE PRECISION    TMP
!!!      INTEGER, ALLOCATABLE :: IPIV(:)
!!!      DOUBLE PRECISION, ALLOCATABLE :: WORK(:)
!!!      ! External
!!!      EXTERNAL           DLACPY, DINVIP
!!!      
!!!      ! INV = X
!!!      CALL DLACPY('A', N, N, X, N, INV, N)
!!!      
!!!      CALL DINVIP(N, INV, INFO)
!!!      
!!!      RETURN
!!!      END


