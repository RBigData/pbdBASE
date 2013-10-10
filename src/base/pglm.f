! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2013, Schmidt


! Logistic regression using IRLS.
!
! X is NxP, locally NXxPX
! Y is Nx1, locally NYx1
! Beta is Px1, locally PBETAx1
!
! ETA = X*BETA_OLD          (Nx1) locally NYx1
! MU = LOGIT_LINKINV(ETA)   (Nx1) locally NYx1
!
! WT is Nx1 locally NYx1
! RESIDS is Nx1 locally NYx1
!
! Each iteration after the first, we fit the linear model:  Z ~ X_TW
! where
!      X_TW = X*WT, and
!      Z = SQRT(W) * (ETA + 1/W * (Y-MU))
!
!
! INCPT = 'Y', 'N', for whether intercept should be included in null model
!
      SUBROUTINE PGLM_FIT(FAMILY, LINK, INCPT, STOPRULE, N, P, 
     $                     X, DESCX, Y, DESCY, BETA, DESCBETA,
     $                     WT, RESIDS, MAXITER, INFO, TOL)
      IMPLICIT NONE
      ! IN/OUT
      CHARACTER*1         INCPT, STOPRULE
      CHARACTER*8         FAMILY, LINK
      INTEGER             N, P, MAXITER, INFO
      INTEGER             DESCX(9), DESCY(9), DESCBETA(9)
      DOUBLE PRECISION    TOL, DEV_OLD
      DOUBLE PRECISION    X(DESCX(9),*), Y(DESCY(9)), BETA(DESCBETA(9)), 
     $                    WT(DESCY(9)), RESIDS(DESCY(9)), 
      ! Local
      INTEGER             CONVERGED, I, J, ITER, ALLOCERR, 
     $                    RANK, LWORK
      INTEGER             LN, LP, NX, NY, PBETA, PX
      INTEGER             LDM(2), BLACS(5)
      DOUBLE PRECISION    AIC, DEV, NULLDEV, TMP
      DOUBLE PRECISION, ALLOCATABLE :: BETA_OLD(:)
      DOUBLE PRECISION, ALLOCATABLE :: ETA(:)
      DOUBLE PRECISION, ALLOCATABLE :: MU(:)
      DOUBLE PRECISION, ALLOCATABLE :: Z(:)
      DOUBLE PRECISION, ALLOCATABLE :: SQWT(:)
      DOUBLE PRECISION, ALLOCATABLE :: X_TW(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: WORK(:)
      ! Parameter
      DOUBLE PRECISION    ZERO, ONE
      PARAMETER ( ZERO = 0.0D0, ONE = 1.0D0 )
      ! External
      INTEGER             GLM_CONVERGENCE, GLM_CHECK_FAM_LINK, 
     $                    PGLM_CHECK_RESPONSE, PGLM_CHECK_MU
      DOUBLE PRECISION    GLM_LOGLIK, GLM_DEVIANCE
      EXTERNAL            PDIMS
      INTRINSIC           MIN, MAX, DBLE, DSQRT
      
      
      ! Quick return if possible
      IF (N.LT.1) THEN
        INFO = -1
        RETURN
      ELSE IF (P.LT.1) THEN
        INFO = -2
        RETURN
      END IF
      
      
      INFO = GLM_CHECK_FAM_LINK(FAMILY, LINK)
      IF (INFO.LT.0) RETURN
      
      INFO = PGLM_CHECK_RESPONSE(N, Y)
      IF (INFO.EQ.-8) RETURN
      
      
      ! Get local and proc grid info
      CALL PDIMS(DESCX, LDM, BLACS)
      NX = LDM(1)
      PX = LDM(2)
      
      NY = LDM(1)
      
      PBETA = DESCBETA(9)
      
      LN = MAX(1, NX)
      LP = MAX(1, PX)
      
!      NY = MAX(1, NY)
      
      
      ! Allocate local storage
      ALLOCERR = 0
      
      ALLOCATE(BETA_OLD(MAX(1, PBETA)), STAT=ALLOCERR)
      IF (ALLOCERR.NE.0) RETURN
      ALLOCATE(ETA(LN), STAT=ALLOCERR)
      IF (ALLOCERR.NE.0) RETURN
      ALLOCATE(SQWT(LN), STAT=ALLOCERR)
      IF (ALLOCERR.NE.0) RETURN
      ALLOCATE(X_TW(LN, LP), STAT=ALLOCERR)
      IF (ALLOCERR.NE.0) RETURN
      ALLOCATE(MU(LN), STAT=ALLOCERR)
      IF (ALLOCERR.NE.0) RETURN
      ALLOCATE(Z(LN), STAT=ALLOCERR)
      IF (ALLOCERR.NE.0) RETURN
      
      
      ! Allocate workspace for linear models
      CALL PDGELS('N', N, P, 1, X, 1, 1, DESCX, Y, 1, 1, DESCY, 
     $            TMP, -1, INFO)
      LWORK = INT(TMP)
      ALLOCATE(WORK(LWORK), STAT=ALLOCERR)
      IF (ALLOCERR.NE.0) RETURN
      
      
      ! Initialize
      DO I = 1, PBETA
        BETA_OLD = ZERO
      END DO
      
      DO I = 1, NY
        WT = ONE
      END DO
      
      DEV = ZERO
      
      
      !!! Main loop
      MAIN: DO ITER = 1, MAXITER
        
        ! Form ETA = X*BETA
        CALL PDGEMM('N', 'N', N, 1, P, ONE, X, 1, 1, DESCX, 
     $              BETA_OLD, 1, 1, DESCBETA, ZERO, ETA, 1, 1, DESCY)
        
        
        ! Compute MU = INVERSE_LINK( ETA )
        IF (ITER.EQ.1) THEN
          CALL GLM_INITIAL_MU(FAMILY, NY, Y, WT, MU)
        ELSE 
          CALL GLM_LINKINV(LINK, NY, ETA, MU)
        END IF
        
        
        ! check for bad fit in the MU's
        INFO = PGLM_CHECK_MU(FAMILY, NY, MU, TOL, INFO)
        IF (INFO.NE.0) GOTO 1
        
        
        ! Update WT = MU*(1-MU)
        CALL GLM_VARIANCE(FAMILY, NY, MU, WT)
        
        DO I = 1, NY
          SQWT(I) = DSQRT(WT(I))
        END DO
        
        IF (STOPRULE.EQ.'3') THEN
          DEV_OLD = DEV
          DEV = GLM_DEVIANCE(FAMILY, N, Y, MU, WT)
        END IF
        
        
        ! Prepare LHS:  X_TW = X*WT
        DO J = 1, LP
          DO I = 1, LN
            X_TW(I,J) = SQWT(I) * X(I,J)
          END DO
        END DO
        
        
        ! Prepare RHS:  Z = SQRT(WT) * (X*BETA + 1/WT*(Y-MU))
        !                 = SQWT * ETA + 1/SQWT*(Y-MU)
        DO I = 1, NY
          TMP = SQWT(I)
          Z(I) = TMP*ETA(I) + ONE/(TMP)*(Y(I)-MU(I))
        END DO
        
        
        ! Update Beta:  Fit Z ~ X_TW
        CALL PGLM_UPDATE_BETA(N, P, BETA, BETA_OLD, X_TW, Z, 
     $                   WORK, LWORK, INFO)
        
        
        ! Check for convergence
        IF (ITER.GT.1) THEN
          CONVERGED = GLM_CONVERGENCE(STOPRULE, P, BETA_OLD, BETA, 
     $                     DEV, DEV_OLD, TOL, ITER, MAXITER, INFO)
        END IF
        
        IF (CONVERGED.EQ.1) THEN
          GOTO 10 ! converged
        ELSE IF (CONVERGED.EQ.2) THEN
          GOTO 1 ! infinite parameter values detected
        END IF
        
      END DO MAIN
      
      
      !!! Success --- now do all the other stuff
 10   CONTINUE
      
      ! AIC, Deviance, NullDeviance
      CALL GLM_LOGLIK_STATS(FAMILY, LINK, INCPT, N, P, X, Y, 
     $                  ETA, MU, BETA, BETA_OLD, 
     $                  DEV, AIC, NULLDEV)
      
      ! Compute working residuals
      CALL GLM_RESIDUALS(LINK, N, Y, MU, ETA, RESIDS)
      
      ITER = ITER-1 ! to agree with R's iteration count
      
      GOTO 1
      
      ! Exit subroutine
 1    CONTINUE
      
      DEALLOCATE(BETA_OLD)
      DEALLOCATE(ETA)
      DEALLOCATE(X_TW)
      DEALLOCATE(MU)
      DEALLOCATE(Z)
      DEALLOCATE(WORK)
      
      
      RETURN
      END


