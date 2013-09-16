! Check for valid MU (as appropriate)
      INTEGER FUNCTION GLM_CHECK_MU(FAMILY, N, MU, TOL)
      IMPLICIT NONE
      ! IN/OUT
      CHARACTER*8         FAMILY
      INTEGER             N, INFO
      DOUBLE PRECISION    MU(*), TOL
      ! Local
      INTEGER             I
      DOUBLE PRECISION    TMP
      ! Parameters
      DOUBLE PRECISION    ONE
      PARAMETER ( ONE = 1.0D0 )
      
      
      IF (FAMILY.EQ.'BINOMIAL') THEN
        TMP = ONE-TOL
        DO I = 1, N
          IF (MU(I).GT.TMP .OR. MU(I).LT.TOL) THEN
            GLM_CHECK_MU = -101
            RETURN
          END IF
        END DO
      
      ELSE IF (FAMILY.EQ.'POISSON' .OR. FAMILY.EQ.'GAMMA') THEN
        DO I = 1, N
          IF (MU(I).LT.TOL) THEN
            GLM_CHECK_MU = -101
            RETURN
          END IF
        END DO
      END IF
      
      
      RETURN
      END



! Check response for garbage values
      INTEGER FUNCTION PGLM_CHECK_RESPONSE(FAMILY, N, Y)
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
            PGLM_CHECK_RESPONSE = FAIL
            RETURN
          END IF
        END DO
      
      ELSE IF (FAMILY.EQ.'POISSON' .OR. FAMILY.EQ.'GAMMA') THEN
        DO I = 1, N
          IF (Y(I).LT.ZERO) THEN
            PGLM_CHECK_RESPONSE = FAIL
            RETURN
          END IF
        END DO
      
      END IF
      
      PGLM_CHECK_RESPONSE = 0
      
      RETURN
      END

