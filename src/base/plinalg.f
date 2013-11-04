! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2013, Schmidt


! X^T * X or X * X^T
! TRANS = 'T' :  X^T*X, TRANS = 'N' : X*X^T
      SUBROUTINE PDCROSSPROD(UPLO, TRANS, ALPHA, X, IX, JX,
     $                       DESCX, C, IC, JC, DESCC)
      IMPLICIT NONE
      ! IN/OUT
      INTEGER             IX, JX, DESCX(9), IC, JC, DESCC(9)
      DOUBLE PRECISION    X( * ), C( * ), ALPHA
      CHARACTER*1         UPLO, TRANS
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
      CALL PDSYRK(UPLO, NST, LDC, LDX, ALPHA, X, IX, JX, DESCX,
     $            ZERO, C, IC, JC, DESCC)
      
      ! Fill lower triangle (make symmetric)
      CALL PDMKSYM(UPLO, C, IC, JC, DESCC)
      
      RETURN
      END 


! compute inverse of a cholesky
      SUBROUTINE PDCHTRI(UPLO, X, IX, JX, DESCX, C, IC, JC, DESCC, INFO)
      IMPLICIT NONE
      ! IN/OUT
      INTEGER             IX, JX, DESCX(9), IC, JC, DESCC(9), INFO
      DOUBLE PRECISION    X( * ), C( * )
      CHARACTER*1         UPLO
      ! Local
      CHARACTER*1         LOUP
      ! Parameter
      DOUBLE PRECISION    ONE
      PARAMETER ( ONE = 1.0D0 )
      ! External
      EXTERNAL            PTRI2ZERO, PDTRTRI, PDCROSSPROD
      
      
      IF (UPLO.EQ.'L') THEN
        LOUP = 'U'
      ELSE IF (UPLO.EQ.'U') THEN
        LOUP = 'L'
      ELSE 
        INFO = -1
        RETURN
      END IF
      
      ! Zero triangle opposite UPLO
      CALL PTRI2ZERO(LOUP, 'N', X, DESCX)
      
      ! Invert the UPLO triangle
      CALL PDTRTRI(UPLO, 'N', DESCX(4), X, IX, JX, DESCX, INFO)
      
      ! 
      CALL PDCROSSPROD(UPLO, 'T', ONE, X, IX, JX, DESCC, 
     $                 C, IC, JC, DESCC)
      
      RETURN
      END 


! Compute matrix inverse without having to understand ScaLAPACK peculiarities
! In place version (X is overwritten with X^-1)
      SUBROUTINE PDINVIP(X, IX, JX, DESCX, INFO)
      IMPLICIT NONE
      ! IN/OUT
      INTEGER             IX, JX, DESCX(9), INFO
      DOUBLE PRECISION    X( * )
      ! Local
      INTEGER             N, LWORK, LIWORK, ALLOCERR
      DOUBLE PRECISION    TMP
      INTEGER, ALLOCATABLE :: IPIV(:), IWORK(:)
      DOUBLE PRECISION, ALLOCATABLE :: WORK(:)
      ! External
      EXTERNAL           PDGETRF, PDGETRI
      
      
      ALLOCERR = 0
      N = DESCX(3)
      
      ! Factor X=LU
      ALLOCATE(IPIV(N + DESCX(6)), STAT=ALLOCERR)
      IF (ALLOCERR.NE.0) RETURN! "Out of memory"
      
      CALL PDGETRF(N, N, X, IX, JX, DESCX, IPIV, INFO)
      IF (INFO.NE.0) RETURN
      
      ! Invert X
      LWORK = -1
      LIWORK = -1
      
      CALL PDGETRI(N, X, IX, JX, DESCX, IPIV, TMP, LWORK, LIWORK, 
     $             LIWORK, INFO)
      IF (INFO.NE.0) RETURN
      
      LWORK = INT(TMP)
      ALLOCATE(WORK(LWORK), STAT=ALLOCERR)
      IF (ALLOCERR.NE.0) RETURN! "Out of memory"
      
      ALLOCATE(IWORK(LIWORK), STAT=ALLOCERR)
      IF (ALLOCERR.NE.0) RETURN! "Out of memory"
      
      CALL PDGETRI(N, X, IX, JX, DESCX, IPIV, WORK, LWORK, IWORK, 
     $             LIWORK, INFO)
      
      DEALLOCATE(IPIV)
      DEALLOCATE(WORK)
      DEALLOCATE(IWORK)
      
      RETURN
      END


! Non-in-place version of matrix inverse (on return, INV = X^-1)
      SUBROUTINE PDINV(X, IX, JX, DESCX, INV, INFO)
      IMPLICIT NONE
      ! IN/OUT
      INTEGER             IX, JX, DESCX(9), INFO
      DOUBLE PRECISION    X( * ), INV( * )
      ! External
      EXTERNAL           PDLACPY, PDINVIP
      
      ! INV = X
      CALL PDLACPY('B', DESCX(3), DESCX(4), X, IX, JX, DESCX, INV,
     $             IX, JX, DESCX)
      
      CALL PDINVIP(INV, IX, JX, DESCX, INFO)
      
      RETURN
      END



! Eigenvalues for non-symmetric matrix
! See http://www.netlib.org/lapack/lug/node50.html for explanation
      SUBROUTINE PDNEP(JOB, X, IX, JX, DESCX, W, INFO)
      IMPLICIT NONE
      ! IN/OUT
      CHARACTER*1         JOB
      INTEGER             IX, JX, DESCX(9), INFO
      DOUBLE PRECISION    X(DESCX(9), *), W(*)
      ! Local
      INTEGER             LTAU, LWORK, ALLOCERR
      INTEGER             I, J
      INTEGER             M, N, GN
      INTEGER             BLACS(5), LDM(2)
      DOUBLE PRECISION    TMP
      DOUBLE PRECISION, ALLOCATABLE :: WORK(:)
      DOUBLE PRECISION, ALLOCATABLE :: TAU(:)
      DOUBLE PRECISION, ALLOCATABLE :: WR(:), WI(:)
      DOUBLE COMPLEX, ALLOCATABLE :: ZX(:,:)
      ! External
      EXTERNAL           PDGEHRD, PDLAHQR
      
      
      ! JOB = 'E' for eigenvalues only, 'A' for vals+vecs
      
      ! Quick return if possible
      GN = DESCX(3)
      IF (N .NE. DESCX(4)) THEN
        INFO = -4
        RETURN
      END IF
      
      ALLOCERR = 0
      
      LTAU = MAX(1, JX+GN-2)
      
      
      
      ! Get local and proc grid info
      CALL PDIMS(DESCX, LDM, BLACS)
      
      M = LDM(1)
      N = LDM(2)
      
      
      
      ALLOCATE(TAU(LTAU), STAT=ALLOCERR)
      IF (ALLOCERR.NE.0) RETURN
      
      !!! Reduce to upper Hessenberg
      ! Workspace query
      CALL PDGEHRD(GN, 1, 1, X, IX, JX, DESCX, TAU, TMP, -1, INFO)
      
      ! Allocate workspace
      LWORK = INT(TMP)
      
      ALLOCATE(WORK(LWORK), STAT=ALLOCERR)
      IF (ALLOCERR.NE.0) RETURN
      
      ! Reduce
      CALL PDGEHRD(GN, 1, 1, X, IX, JX, DESCX, TAU, WORK, LWORK, INFO)
      
      IF (INFO .NE. 0) GOTO 1
      
      
      !!! Schur Decomposition
      ! Workspace query
      CALL PDLAHQR(.FALSE., .FALSE., GN, 1, N, X, DESCX, WR, WI, 0, 0, 
     $             0.0D0, DESCX, TMP, -1, 0, 0, INFO)
      
      ! Reallocate workspace as needed
      IF (INT(TMP) .GT. LWORK) THEN
        DEALLOCATE(WORK)
        LWORK = INT(TMP)
        ALLOCATE(WORK(LWORK), STAT=ALLOCERR)
        IF (ALLOCERR.NE.0) STOP "Out of memory"
      END IF
      
      ! Compute
      CALL PDLAHQR(.FALSE., .FALSE., GN, 1, N, X, DESCX, WR, WI, 0, 0, 
     $             0.0D0, DESCX, WORK, LWORK, 0, 0, INFO)
      
      
      !!! Store Eigenvalues from the Schur decomposition
      
      !!! Compute Eigenvectors
      
      IF (JOB.EQ.'A') THEN
        ALLOCATE(ZX(M, N))
        
        DO J = 1, N
          DO I = 1, N
            ZX(I, J) = DCMPLX(X(I, J))
          END DO
        END DO
        
!        CALL PZTREVC(SIDE, HOWMNY, SELECT, GN, T, DESCT, VL, DESCVL, 
!     $               VR, DESCVR, MM, M, WORK, RWORK, INFO )
      
      END IF
      
      !!! Deallocate workspace and exit
    1 CONTINUE
      DEALLOCATE(WORK)
      DEALLOCATE(TAU)
      
      RETURN
      END



