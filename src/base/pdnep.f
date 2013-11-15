! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2013, Schmidt

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

