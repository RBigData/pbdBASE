! Update BETA_OLD = BETA and BETA = solution to Y~X*BETA
      SUBROUTINE PGLM_UPDATE_BETA(N, P, BETA, PBETA, BETA_OLD, 
     $                            X, DESCX, Y, NY, DESCY, 
     $                            WORK, LWORK, INFO)
      IMPLICIT NONE
      ! IN/OUT
      INTEGER             N, P, LWORK, INFO
      DOUBLE PRECISION    BETA(P), BETA_OLD(P), X(N,P), Y(N), 
     $                    WORK(LWORK)
      ! Local
      INTEGER             I, K
      ! External
      EXTERNAL            DGELS
      INTRINSIC           MIN
      
      
      CALL PDGELS('N', N, P, 1, X, 1, 1, DESCX, Y, 1, 1, DESCY, WORK, 
     $           LWORK, INFO)
      
      K = MIN(NY, PBETA)
      
      DO I = 1, PBETA
        BETA_OLD(I) = BETA(I)
        BETA(I) = Y(I)
      END DO
      
      
      RETURN
      END



! Check for convergence based on chosen stoprule.  Returns:
    ! 1 - converged
    ! 2 - infinite params
    ! 3 - no improvement
      INTEGER FUNCTION GLM_CONVERGENCE(STOPRULE, P, BETA_OLD, BETA,
     $                    DEV, DEV_OLD, TOL, ITER, MAXITER, INFO)
      IMPLICIT NONE
      ! IN/OUT
      CHARACTER*1         STOPRULE
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
      
      IF (STOPRULE.EQ.'1') THEN
        RETURN
      ELSE IF (STOPRULE.EQ.'2') THEN
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
      ELSE IF (STOPRULE.EQ.'3') THEN
        TMP1 = DABS(DEV-DEV_OLD)
        TMP2 = (0.1D0 + DABS(DEV))
        
        
!        IF (DEV.GT.0.0D0) THEN
          IF (TMP1/TMP2.LT.TOL) THEN
            GLM_CONVERGENCE = 1
          ELSE
            GLM_CONVERGENCE = -1
          END IF
!        END IF
      END IF
      
      RETURN
      END



! Initialize MU based on the error distribution
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


! Loglikelihood calculator
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



! Model deviance calculator
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



! Null deviance calculator
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



! Model statistics wrapper (AIC, deviance, and null deviance)
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



! Link function
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



! Inverse link function
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



! Variance
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



! Residuals
      SUBROUTINE GLM_RESIDUALS(LINK, N, Y, MU, ETA, RESID)
      IMPLICIT NONE
      ! IN/OUT
      CHARACTER*8         LINK
      INTEGER             N
      DOUBLE PRECISION    Y(*), MU(*), ETA(*), RESID(*)
      ! Local
      INTEGER             I
      DOUBLE PRECISION    TMP
      ! Parameter
      DOUBLE PRECISION   ONE, TWO, NEGONE
      PARAMETER ( ONE = 1.0D0, TWO = 2.0D0, NEGONE = -1.0D0 )
      ! Intrinsic
      INTRINSIC           DEXP
      
      
      !!! Working residuals
      IF (LINK.EQ.'LOGIT') THEN
        DO I = 1, N
          TMP = MU(I)
          RESID(I) = (Y(I) -TMP) / TMP / (ONE - TMP)
        END DO
      
      ELSE IF (LINK.EQ.'LOG') THEN
        DO I = 1, N
          RESID(I) = (Y(I) - MU(I)) / MU(I)
        END DO
      
      ELSE IF (LINK.EQ.'IDENTITY') THEN
        DO I = 1, N
          RESID(I) = Y(I) - MU(I)
        END DO
      
      ELSE IF (LINK.EQ.'INVERSE') THEN
        DO I = 1, N
          TMP = ETA(I)
          RESID(I) = NEGONE * TMP*TMP * (Y(I) - MU(I))
        END DO
      
      ELSE IF (LINK.EQ.'SQRT') THEN
        DO I = 1, N
          RESID(I) = (Y(I) - MU(I)) / (TWO * ETA(I))
        END DO
      END IF
      
      RETURN
      END


