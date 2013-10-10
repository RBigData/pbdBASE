! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2013, Schmidt


! GLM_FIT
! 
! Purpose
! =======
!
! Generalized Linear Models (GLM) using IRLS.
!
! Implementation is based on "Generalized Linear Models, Second 
! Edition", by P. McCullagh and J. Nelder.
!
! Arguments
! =========
!
! FAMILY    (input) CHARACTER*8
!           The distribution of the response.  Choices are:
!           BINOMIAL, GAMMA, GAUSSIAN, POISSON.
! 
! LINK      (input) CHARACTER*8
!           The link function to use.  Choice depends on the
!           distribution.  Choices are:
!           (BINOMIAL) LOGIT, PROBIT, LOG, CLOGLOG,
!           (GAMMA) INVERSE, IDENTITY, LOG,
!           (GAUSSIAN) IDENTITY, LOG, INVERSE
!           (POISSON) LOG, IDENTITY, SQRT
! 
! INCPT     (input) CHARACTER*1
!           
! 
! STOPRULE  (input) CHARACTER*1
!           
! 
! N         (input) INTEGER
!           The number of rows of the data matrix X.  N>=0.
! 
! P         (input) INTEGER
!           The number of columns of the data matrix X.  P>=0.
! 
! X         (input/output) DOUBLE PRECISION array, dimension (N,P).
!           On entry, the input data matrix.  On exit, the output
!           of a LAPACK QR factorization is stored for X.
! 
! Y         (input) DOUBLE PRECISION array, dimension (N,1).
!           The response variable.
! 
! OFFSET    (input)
!           
! 
! BETA      (output) DOUBLE PRECISION array, dimension (P,1).
!           The model coefficients.
! 
! WT        (input/output) DOUBLE PRECISION array, dimension (N,1).
!           On input, a set of starting weights (initialize to 1.0
!           if you have nothing else in mind).  On output, the final
!           set of weights ("working weights") in the IRLS fit.
! 
! OFFSET    (input) DOUBLE PRECISION array, dimension (N,1).
!           
! 
! RESIDS    (output)
!           
! 
! MAXITER   (input) INTEGER
!           Maximum number of iterations.  20 should be sufficient
!           for convergence for most applications if it's going to 
!           happen at all.
! 
! TOL       (input) DOUBLE PRECISION
!           Tolerance for the QR decomposition.
! 
! INFO      (output) INTEGER
!           = 0: successful exit.
!           < 0: if INFO = -i, the i-th argument had an illegal value.


!
! ETA = X*BETA_OLD
! MU = LOGIT_LINKINV(ETA)
!
!
! Each iteration after the first, we fit the linear model:  Z ~ X_TW
! where
!      X_TW = X*WT, and
!      Z = SQRT(W) * (ETA + 1/W * (Y-MU))
!
!
! INCPT = 'Y', 'N', for whether intercept should be included in null model
!
      SUBROUTINE GLM_FIT(FAMILY, LINK, INCPT, STOPRULE, N, P, X, Y,
     $                BETA, WT, OFFSET, RESIDS, MAXITER, TOL, INFO)
      IMPLICIT NONE
      ! IN/OUT
      CHARACTER*8         FAMILY, LINK
      CHARACTER*1         INCPT
      INTEGER             STOPRULE, N, P, MAXITER, INFO
      DOUBLE PRECISION    X(N,P), Y(N), OFFSET(N), BETA(P), WT(N),
     $                    RESIDS(N), TOL
      ! Local
      INTEGER             CONVERGED, I, J, ITER, ALLOCERR, 
     $                    RANK, LWORK
      DOUBLE PRECISION    AIC, DEV, NULLDEV, TMP, DEV_OLD
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
      INTEGER             GLM_CONVERGENCE, CHECK_FAM_LINK, 
     $                    CHECK_RESPONSE
      DOUBLE PRECISION    GLM_LOGLIK, GLM_DEVIANCE
      EXTERNAL            GLM_UPDATE_BETA
      INTRINSIC           MIN, MAX, DBLE, DSQRT
      
      
      ! Quick return if possible
      INFO = CHECK_FAM_LINK(FAMILY, LINK)
      IF (INFO.LT.0) RETURN
      
      INFO = CHECK_RESPONSE(N, Y)
      IF (INFO.EQ.-8) RETURN
      
      ! Allocate local storage
      ALLOCERR = 0
      
      ALLOCATE(BETA_OLD(P), STAT=ALLOCERR)
      IF (ALLOCERR.NE.0) STOP "Out of memory"
      ALLOCATE(ETA(N), STAT=ALLOCERR)
      IF (ALLOCERR.NE.0) STOP "Out of memory"
      ALLOCATE(SQWT(N), STAT=ALLOCERR)
      IF (ALLOCERR.NE.0) STOP "Out of memory"
      ALLOCATE(X_TW(N, P), STAT=ALLOCERR)
      IF (ALLOCERR.NE.0) STOP "Out of memory"
      ALLOCATE(MU(N), STAT=ALLOCERR)
      IF (ALLOCERR.NE.0) STOP "Out of memory"
      ALLOCATE(Z(N), STAT=ALLOCERR)
      IF (ALLOCERR.NE.0) STOP "Out of memory"
      
      ! Allocate workspace for linear models
      LWORK = MIN(N, P) + MAX(1, N, P)
      ALLOCATE(WORK(LWORK), STAT=ALLOCERR)
      IF (ALLOCERR.NE.0) STOP "Out of memory"
      
      
      ! Empty model case
      IF (N.LT.1 .OR. P.LT.1) THEN
        DO I = 1, N
          ETA(I) = ZERO + OFFSET(I)
        END DO
        
        CALL GLM_LINKINV(LINK, N, ETA, MU)
        
        CALL GLM_CHECK_MU(FAMILY, N, MU, TOL, INFO)
        IF (INFO.NE.0) RETURN
        
        CALL GLM_VARIANCE(FAMILY, N, MU, WT)
        
        CALL GLM_RESIDUALS(LINK, N, Y, MU, ETA, RESIDS)
        
        ITER = 0
        
        RETURN
      END IF
      
      
      ! Initialize
      DO I = 1, P
        BETA_OLD(I) = ZERO
      END DO
      
      DO I = 1, N
        WT(I) = ONE
      END DO
      
      DEV = ZERO
      
      
      !!! Main loop
      MAIN: DO ITER = 0, MAXITER
        
        
        ! Compute ETA = X*BETA and MU = INVERSE_LINK( ETA )
        IF (ITER.EQ.0) THEN
          CALL GLM_INITIAL_MU(FAMILY, N, Y, WT, MU)
          CALL GLM_LINK(LINK, N, MU, ETA)
        ELSE 
          CALL DGEMM('N', 'N', N, 1, P, ONE, X, N, BETA, P,  ZERO, 
     $                ETA, N)
          CALL GLM_LINKINV(LINK, N, ETA, MU)
        END IF
        
        
        ! check for bad fit in the MU's
        CALL GLM_CHECK_MU(FAMILY, N, MU, TOL, INFO)
        IF (INFO.NE.0) RETURN
        
        
        ! Update WT
        CALL GLM_VARIANCE(FAMILY, N, MU, WT)
        
        DO I = 1, N
          SQWT(I) = DSQRT(WT(I))
        END DO
        
        IF (STOPRULE.EQ.3) THEN
          DEV_OLD = DEV
          DEV = GLM_DEVIANCE(FAMILY, N, Y, MU, WT)
        END IF
        
        
        ! Prepare LHS:  X_TW = X*WT
        DO J = 1, P
          DO I = 1, N
            X_TW(I,J) = SQWT(I) * X(I,J)
          END DO
        END DO
        
        
        ! Prepare RHS:  Z = SQRT(WT) * (X*BETA + 1/WT*(Y-MU))
        !                 = SQWT * ETA + 1/SQWT*(Y-MU)
        DO I = 1, N
          TMP = SQWT(I)
          Z(I) = TMP*ETA(I) + ONE/(TMP)*(Y(I)-MU(I))
        END DO
        
        
        ! Update Beta:  Fit Z ~ X_TW
        CALL GLM_UPDATE_BETA(N, P, BETA, BETA_OLD, X_TW, Z, 
     $                   WORK, LWORK, INFO)
        
        
        ! Check for convergence
        IF (ITER.GT.0) THEN
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
      
      WRITE (*,*) "iter=",ITER
      
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




! ==================================
! Internal functions used in GLM_FIT
!
! Consider these "self-documenting" aka I'm too lazy to explain 
! this mess

      SUBROUTINE GLM_UPDATE_BETA(N, P, BETA, BETA_OLD, X, Y, 
     $                       WORK, LWORK, INFO)
      IMPLICIT NONE
      ! IN/OUT
      INTEGER             N, P, LWORK, INFO
      DOUBLE PRECISION    BETA(P), BETA_OLD(P), X(N,P), Y(N), 
     $                    WORK(LWORK)
      ! Local
      INTEGER             K, I
      ! External
      EXTERNAL            DGELS
      INTRINSIC           MIN
      
      
      CALL DGELS('N', N, P, 1, X, N, Y, N, WORK, LWORK, INFO)
      
      K = MIN(N, P)
      
      DO I = 1, K
        BETA_OLD(I) = BETA(I)
        BETA(I) = Y(I)
      END DO
      
      RETURN
      END





      SUBROUTINE GLM_CHECK_MU(FAMILY, N, MU, TOL, INFO)
      IMPLICIT NONE
      ! IN/OUT
      CHARACTER*8         FAMILY
      INTEGER             N, INFO
      DOUBLE PRECISION    MU(*), TOL
      ! Local
      INTEGER             I
      DOUBLE PRECISION    TMP
      ! Parameters
      DOUBLE PRECISION    ONE, ZERO
      PARAMETER ( ONE = 1.0D0, ZERO = 0.0D0 )
      
      
      IF (FAMILY.EQ.'BINOMIAL') THEN
        TMP = ONE-TOL
        DO I = 1, N
          IF (MU(I).GT.TMP .OR. MU(I).LT.TOL) THEN
            INFO = -101
            RETURN
          END IF
        END DO
      
      ELSE IF (FAMILY.EQ.'POISSON' .OR. FAMILY.EQ.'GAMMA') THEN
        DO I = 1, N
          IF (MU(I).LT.TOL) THEN
            INFO = -101
            RETURN
          END IF
        END DO
      END IF
      
      
      RETURN
      END





! 1 - converged
! 2 - infinite params
! 3 - no improvement
      INTEGER FUNCTION GLM_CONVERGENCE(STOPRULE, P, BETA_OLD, BETA,
     $                    DEV, DEV_OLD, TOL, ITER, MAXITER, INFO)
      IMPLICIT NONE
      ! IN/OUT
      INTEGER             STOPRULE
      INTEGER             P, ITER, MAXITER, INFO
      DOUBLE PRECISION    BETA_OLD(*), BETA(*), DEV, DEV_OLD, TOL
      ! Local
      INTEGER             I
      DOUBLE PRECISION    TMP1, TMP2
      ! Intrinsic
      INTRINSIC           DABS
      LOGICAL             DISNAN
      
      
      GLM_CONVERGENCE = 0
      
      ! Check that all parameters are finite and 
      DO I = 1, P
        IF(DISNAN(BETA(I))) THEN
          GLM_CONVERGENCE = 2
          INFO = 1
          GOTO 1
!        ELSE IF (.NOT.IEEE_IS_FINITE(BETA(I))) THEN
!          GLM_CONVERGENCE = 2
!          GOTO 1
        END IF
      END DO
      
      IF (STOPRULE.EQ.1) THEN
        RETURN
      ELSE IF (STOPRULE.EQ.2) THEN
        DO I = 1, P
          TMP1 = DABS(BETA(I) - BETA_OLD(I))
          TMP2 = TOL*1.0D-6 + TOL*DABS(BETA_OLD(I))
          
          IF (TMP1.GT.TMP2) THEN
            GLM_CONVERGENCE = -1
            IF (ITER.EQ.MAXITER) INFO = MAXITER
            GOTO 1
          END IF
        END DO
        
        GLM_CONVERGENCE = 1
        
 1      CONTINUE
      ELSE IF (STOPRULE.EQ.3) THEN
        TMP1 = DABS(DEV-DEV_OLD)
        TMP2 = (0.1D0 + DABS(DEV))
        
        
        IF (TMP1/TMP2.LT.TOL) THEN
          GLM_CONVERGENCE = 1
        ELSE
          GLM_CONVERGENCE = -1
        END IF
      END IF
      
      RETURN
      END





      SUBROUTINE GLM_INITIAL_MU(FAMILY, N, Y, WT, MU)
      IMPLICIT NONE
      ! IN/OUT
      CHARACTER*8         FAMILY
      INTEGER             N, P, LWORK, INFO
      DOUBLE PRECISION    Y(*), WT(*), MU(*)
      ! Local
      INTEGER             I
      ! Parameter
      DOUBLE PRECISION    ONE, TENTH, HALF
      PARAMETER ( ONE = 1.0D0, TENTH = 0.1D0, HALF = 0.5D0 )
      ! External
      EXTERNAL            DGELS
      
      
      IF (FAMILY.EQ.'BINOMIAL') THEN
        DO I = 1, N
          MU(I) = (WT(I) * Y(I) + HALF) / (WT(I) + ONE)
        END DO
      
      ELSE IF (FAMILY.EQ.'GAMMA' .OR. FAMILY.EQ.'GAUSSIAN') THEN
        DO I = 1, N
          MU(I) = Y(I)
        END DO
      
      ELSE IF (FAMILY.EQ.'POISSON') THEN
        DO I = 1, N
          MU(I) = Y(I) + TENTH
        END DO
      END IF
      
      RETURN
      END





      DOUBLE PRECISION FUNCTION GLM_LOGLIK(FAMILY, N, Y, MU)
      IMPLICIT NONE
      ! IN/OUT
      CHARACTER*8         FAMILY
      INTEGER             N
      DOUBLE PRECISION    MU(N), Y(N)
      ! Local
      INTEGER             I
      DOUBLE PRECISION    TMP
      ! Parameter
      DOUBLE PRECISION    ZERO, ONE, TWO
      PARAMETER ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0 )
      
      
      GLM_LOGLIK = ZERO
      
      !!! Gaussian
      IF (FAMILY.EQ.'GAUSSIAN') THEN
        DO I = 1, N
          TMP = Y(I) - MU(I)
          GLM_LOGLIK = GLM_LOGLIK - TMP*TMP
        END DO
        
        GLM_LOGLIK = GLM_LOGLIK / TWO
      
      !!! Poisson
      ELSE IF (FAMILY.EQ.'POISSON') THEN
        DO I = 1, N
          IF (Y(I).GT.ZERO) THEN
            GLM_LOGLIK = GLM_LOGLIK - Y(I)*DLOG(Y(I)/MU(I)) 
          END IF
          GLM_LOGLIK = GLM_LOGLIK + Y(I) - MU(I)
        END DO
        
      !!! Binomial
      ELSE IF (FAMILY.EQ.'BINOMIAL') THEN
        DO I = 1, N
          IF (Y(I).GT.ZERO) THEN
            GLM_LOGLIK = GLM_LOGLIK + DLOG(MU(I))
          ELSE
            GLM_LOGLIK = GLM_LOGLIK + DLOG(ONE-MU(I))
          END IF
        END DO
        
      ELSE IF (FAMILY.EQ.'GAMMA') THEN
        DO I = 1, N
          IF (Y(I).GT.ZERO) THEN
            GLM_LOGLIK = GLM_LOGLIK + DLOG(Y(I)/MU(I))
          END IF
          GLM_LOGLIK = GLM_LOGLIK - (Y(I) - MU(I))/MU(I)
        END DO
        
      END IF
      
      RETURN
      END





      DOUBLE PRECISION FUNCTION GLM_DEVIANCE(FAMILY, N, Y, MU)
      IMPLICIT NONE
      ! IN/OUT
      CHARACTER*8         FAMILY
      INTEGER             N
      DOUBLE PRECISION    MU(N), Y(N)
      ! Local
      INTEGER             I
      DOUBLE PRECISION    TMP
      ! Parameter
      DOUBLE PRECISION    ZERO, ONE, TWO
      PARAMETER ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0 )
      
      
      GLM_DEVIANCE = ZERO
      
      !!! Normal
      IF (FAMILY.EQ.'GAUSSIAN') THEN
        DO I = 1, N
          TMP = Y(I) - MU(I)
          GLM_DEVIANCE = GLM_DEVIANCE + TMP*TMP
        END DO
        
      !!! Poisson
      ELSE IF (FAMILY.EQ.'POISSON') THEN
        DO I = 1, N
          IF (Y(I).GT.ZERO) THEN
            GLM_DEVIANCE = GLM_DEVIANCE + Y(I)*DLOG(Y(I)/MU(I)) 
          END IF
          GLM_DEVIANCE = GLM_DEVIANCE + MU(I) - Y(I)
        END DO
        
        GLM_DEVIANCE = TWO*GLM_DEVIANCE
        
      !!! Binomial
      ELSE IF (FAMILY.EQ.'BINOMIAL') THEN
        DO I = 1, N
          IF (Y(I).GT.ZERO) THEN
            GLM_DEVIANCE = GLM_DEVIANCE + DLOG(MU(I))
          ELSE
            GLM_DEVIANCE = GLM_DEVIANCE + DLOG(ONE-MU(I))
          END IF
        END DO
        
        GLM_DEVIANCE = -TWO*GLM_DEVIANCE
      
      !!! GAMMA
      ELSE IF (FAMILY.EQ.'GAMMA') THEN
        DO I = 1, N
          IF (Y(I).GT.ZERO) THEN
            GLM_DEVIANCE = GLM_DEVIANCE - DLOG(Y(I)/MU(I))
          END IF
          GLM_DEVIANCE = GLM_DEVIANCE + (Y(I) - MU(I))/MU(I)
        END DO
        
        GLM_DEVIANCE = -TWO*GLM_DEVIANCE
      
      END IF
      
      RETURN
      END





      DOUBLE PRECISION FUNCTION GLM_NULLDEV(FAMILY, N, Y, MU)
      IMPLICIT NONE
      ! IN/OUT
      CHARACTER*8         FAMILY
      INTEGER             N
      DOUBLE PRECISION    MU(*), Y(*)
      ! Local
      INTEGER             I
      DOUBLE PRECISION    TMP
      ! Parameter
      DOUBLE PRECISION    ZERO, ONE, TWO
      PARAMETER ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0 )
      ! Functions
      DOUBLE PRECISION    GLM_DEVIANCE
      
      
      TMP = ZERO
      DO I = 1, N
        TMP = TMP + Y(I)/N
      END DO
      
      DO I = 1, N
        MU(I) = TMP
      END DO
      
      GLM_NULLDEV = GLM_DEVIANCE(FAMILY, N, Y, MU)
        
      
      RETURN
      END





      SUBROUTINE GLM_LOGLIK_STATS(FAMILY, LINK, INCPT, N, P, X, Y,
     $                            ETA, MU, BETA, BETA_OLD, 
     $                            DEV, AIC, NULLDEV)
      IMPLICIT NONE
      ! IN/OUT
      CHARACTER*1         INCPT
      CHARACTER*8         FAMILY, LINK
      INTEGER             N, P
      DOUBLE PRECISION    X(*), Y(N), ETA(N), MU(*),
     $                    BETA(*), BETA_OLD(P),
     $                    DEV, AIC, NULLDEV
      ! Local
      INTEGER             I
      DOUBLE PRECISION    TMP
      ! Parameter
      DOUBLE PRECISION    ZERO, ONE, TWO, NEGTWO
      PARAMETER ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0, 
     $            NEGTWO = -2.0D0 )
      ! Intrinsic
      DOUBLE PRECISION    GLM_LOGLIK, GLM_DEVIANCE, GLM_NULLDEV
      INTRINSIC           DBLE, SUM
      
      
      ! Model deviance
!      DEV = NEGTWO * GLM_LOGLIK(FAMILY, N, Y, MU)
      DEV = GLM_DEVIANCE(FAMILY, N, Y, MU)
      
      ! Model AIC
      AIC = TWO * DBLE(P) + DEV
      
      !!! null deviance
      ! Null deviance for model with no intercept
      IF (INCPT.EQ.'N') THEN
        DO I = 1, N
          ETA(I) = ZERO
        END DO
        
        CALL GLM_LINKINV(LINK, N, ETA, MU)
        NULLDEV = NEGTWO * GLM_LOGLIK(FAMILY, N, Y, MU)
      
      ! Null deviance (deviance for intercept-only model)
      ELSE
        NULLDEV = GLM_NULLDEV(FAMILY, N, Y, MU)
        
      END IF
      
      
      RETURN
      END





      SUBROUTINE GLM_LINK(LINK, N, X, Y)
      IMPLICIT NONE
      ! IN/OUT
      CHARACTER*8         LINK
      INTEGER             N
      DOUBLE PRECISION    X(*), Y(*)
      ! Local
      INTEGER             I
      DOUBLE PRECISION    TMP
      ! Parameter
      DOUBLE PRECISION    ONE
      PARAMETER ( ONE = 1.0D0 )
      ! Intrinsic
      INTRINSIC           DLOG, DSQRT
      
      
      IF (LINK.EQ.'CLOGLOG') THEN
        DO I = 1, N
          Y(I) = DLOG(-DLOG(ONE-X(I)))
        END DO
      
      ELSE IF (LINK.EQ.'IDENTITY') THEN
        DO I = 1, N
          Y(I) = X(I)
        END DO
      
      ELSE IF (LINK.EQ.'INVERSE') THEN
        DO I = 1, N
          Y(I) = ONE/X(I)
        END DO
      
      ELSE IF (LINK.EQ.'LOG') THEN
        DO I = 1, N
          Y(I) = DLOG(X(I))
        END DO
      
      ELSE IF (LINK.EQ.'LOGIT') THEN
        DO I = 1, N
          TMP = X(I)
          Y(I) = DLOG(TMP / (ONE-TMP))
        END DO
      
      ELSE IF (LINK.EQ.'SQRT') THEN
        Y(I) = DSQRT(X(I))
      
      END IF
      
      RETURN
      END





      SUBROUTINE GLM_LINKINV(LINK, N, X, Y)
      IMPLICIT NONE
      ! IN/OUT
      CHARACTER*8         LINK
      INTEGER             N
      DOUBLE PRECISION    X(*), Y(*)
      ! Local
      INTEGER             I
      DOUBLE PRECISION    TMP
      ! Parameter
      DOUBLE PRECISION    ZERO, ONE
      PARAMETER ( ZERO = 0.0D0, ONE = 1.0D0 )
      ! Intrinsic
      INTRINSIC           DEXP, DSQRT
      
      
      IF (LINK.EQ.'CLOGLOG') THEN
        DO I = 1, N
          TMP = -DEXP(X(I))
          Y(I) = -DEXP(TMP)-ONE
        END DO
      
      ELSE IF (LINK.EQ.'IDENTITY') THEN
        DO I = 1, N
          Y(I) = X(I)
        END DO
      
      ELSE IF (LINK.EQ.'INVERSE') THEN
        DO I = 1, N
          IF (X(I).GT.ZERO) THEN
            Y(I) = ONE/X(I)
          ELSE
            Y(I) = ZERO
          END IF
        END DO
        
      ELSE IF (LINK.EQ.'LOG') THEN
          DO I = 1, N
            Y(I) = DEXP(X(I))
          END DO
      
      ELSE IF (LINK.EQ.'LOGIT') THEN
        DO I = 1, N
          TMP = DEXP(X(I))
          Y(I) = TMP / (ONE + TMP)
        END DO
      
      ELSE IF (LINK.EQ.'SQRT') THEN
        DO I = 1, N
          TMP = X(I)
          Y(I) = TMP*TMP
        END DO
      
      END IF
      
      RETURN
      END





      SUBROUTINE GLM_VARIANCE(FAMILY, N, MU, VAR)
      IMPLICIT NONE
      ! IN/OUT
      CHARACTER*8         FAMILY
      INTEGER             N
      DOUBLE PRECISION    MU(N), VAR(N)
      ! Local
      INTEGER             I
      DOUBLE PRECISION    TMP
      ! Parameter
      DOUBLE PRECISION    ONE
      PARAMETER ( ONE = 1.0D0 )
      
      
      IF (FAMILY.EQ.'BINOMIAL') THEN
          DO I = 1, N
            TMP = MU(I)
            VAR(I) = TMP * (ONE - TMP)
          END DO
      
      ELSE IF (FAMILY.EQ.'GAMMA') THEN
        DO I = 1, N
          TMP = MU(I)
          VAR(I) = TMP*TMP
        END DO
      
      ELSE IF (FAMILY.EQ.'GAUSSIAN') THEN
        DO I = 1, N
          VAR(I) = ONE
        END DO
      
      ELSE IF (FAMILY.EQ.'POISSON') THEN
        DO I = 1, N
          VAR(I) = MU(I)
        END DO
      END IF
      
      RETURN
      END





! 0 is ok, -1 means family not supported, -2 means link not supported
      INTEGER FUNCTION CHECK_FAM_LINK(FAMILY, LINK)
      ! IN/OUT
      CHARACTER*8         FAMILY, LINK
      ! Parameters
      INTEGER             BAD_FAM, BAD_LINK
      PARAMETER ( BAD_FAM = -1, BAD_LINK = -2 )
      
      
      IF (FAMILY.EQ.'BINOMIAL') THEN
        IF (LINK.NE.'CLOGLOG'  .AND. 
     $      LINK.NE.'LOG'      .AND.
     $      LINK.NE.'LOGIT')    THEN
          CHECK_FAM_LINK = BAD_LINK
        END IF
      
      ELSE IF (FAMILY.EQ.'GAMMA') THEN
        IF (LINK.NE.'IDENTITY'   .AND.
     $      LINK.NE.'LOG'        .AND. 
     $      LINK.NE.'INVERSE')   THEN
          CHECK_FAM_LINK = BAD_LINK
        END IF
      
      ELSE IF (FAMILY.EQ.'GAUSSIAN') THEN
        IF (LINK.NE.'IDENTITY'      .AND.
     $      LINK.NE.'LOG'           .AND. 
     $      LINK.NE.'INVERSE')       THEN
          CHECK_FAM_LINK = BAD_LINK
        END IF
      
      ELSE IF (FAMILY.EQ.'POISSON') THEN
        IF (LINK.NE.'IDENTITY'     .AND.
     $      LINK.NE.'LOG'          .AND. 
     $      LINK.NE.'SQRT')         THEN
          CHECK_FAM_LINK = BAD_LINK
        END IF
      
      ELSE
        ! Family not supported
        CHECK_FAM_LINK = BAD_FAM
      END IF
      
      RETURN
      END





      INTEGER FUNCTION CHECK_RESPONSE(FAMILY, N, Y)
      ! IN/OUT
      CHARACTER*8         FAMILY
      INTEGER             N
      DOUBLE PRECISION    Y(*)
      ! Local
      INTEGER             I
      ! Parameter
      INTEGER             FAIL
      DOUBLE PRECISION    ZERO, ONE
      PARAMETER ( FAIL = -8, ZERO = 0.0D0, ONE = 1.0D0 )
      
      
      IF (FAMILY.EQ.'BINOMIAL') THEN
        DO I = 1, N
          IF (Y(I).LT.ZERO .OR. Y(I).GT.ONE) THEN
            CHECK_RESPONSE = FAIL
            RETURN
          END IF
        END DO
      
      ELSE IF (FAMILY.EQ.'POISSON' .OR. FAMILY.EQ.'GAMMA') THEN
        DO I = 1, N
          IF (Y(I).LT.ZERO) THEN
            CHECK_RESPONSE = FAIL
            RETURN
          END IF
        END DO
      
      END IF
      
      CHECK_RESPONSE = 0
      
      RETURN
      END





      SUBROUTINE GLM_RESIDUALS(LINK, N, Y, MU, ETA, RESIDS)
      IMPLICIT NONE
      ! IN/OUT
      CHARACTER*8         LINK
      INTEGER             N
      DOUBLE PRECISION    Y(*), MU(*), ETA(*), RESIDS(*)
      ! Local
      INTEGER             I
      DOUBLE PRECISION    TMP
      ! Parameter
      DOUBLE PRECISION   ONE, TWO, NEGONE
      PARAMETER ( ONE = 1.0D0, TWO = 2.0D0, NEGONE = -1.0D0 )
      ! Intrinsic
      INTRINSIC           DEXP
      
      
      ! "Working" residuals
      
      IF (LINK.EQ.'CLOGLOG') THEN
        DO I = 1, N
          TMP = DEXP(ETA(I))
          RESIDS(I) = (Y(I) - MU(I)) / (TMP * DEXP(-TMP))
        END DO
      
      ELSE IF (LINK.EQ.'IDENTITY') THEN
        DO I = 1, N
          RESIDS(I) = Y(I) - MU(I)
        END DO
      
      ELSE IF (LINK.EQ.'INVERSE') THEN
        DO I = 1, N
          TMP = ETA(I)
          RESIDS(I) = NEGONE * TMP*TMP * (Y(I) - MU(I))
        END DO
      
      ELSE IF (LINK.EQ.'LOG') THEN
        DO I = 1, N
          RESIDS(I) = (Y(I) - MU(I)) / MU(I)
        END DO
      
      ELSE IF (LINK.EQ.'LOGIT') THEN
        DO I = 1, N
          TMP = MU(I)
          RESIDS(I) = (Y(I) - TMP) / TMP / (ONE - TMP)
        END DO
      
      ELSE IF (LINK.EQ.'SQRT') THEN
        DO I = 1, N
          RESIDS(I) = (Y(I) - MU(I)) / (TWO * ETA(I))
        END DO
      END IF
      
      RETURN
      END

