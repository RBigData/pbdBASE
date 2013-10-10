! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2013, Schmidt


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


