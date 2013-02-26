!! --------------------------------------------------------------------
!!     Calculate mean, column-wise of a distributed matrix object.
!!     Somewhat equivalent to colMeans(X) in R.
!! --------------------------------------------------------------------
      SUBROUTINE DDMATMN( X, M, LCM, LCN, ICTXT, MN)
      ! Inputs
      INTEGER             M, LCM, LCN, ICTXT
      DOUBLE PRECISION    X( LCM, LCN )
      ! Output
      DOUBLE PRECISION    MN( LCN )
      ! Local
      INTEGER             I, J
      DOUBLE PRECISION    MEAN
      ! Subroutines
      EXTERNAL            DGSUM2D
      
      ! --------------------------
      ! Mean calculator
      ! --------------------------
      
      DO I = 1, LCN
        MN(I) = 0
        DO J = 1, LCM
          MN(I) = MN(I) + (X(J,I) / M)
        END DO ! J
      END DO ! I
      
      CALL DGSUM2D( ICTXT, 'Column', ' ', LCN, 1, MN, 1, -1, -1 )
      
      RETURN
      END


!! --------------------------------------------------------------------
!!     Calculate variance, column-wise of a distributed matrix object.
!!     Somewhat equivalent to calling var(X) in R. Probably slower, but
!!     more memory efficient; possibly more numerically stable.
!! --------------------------------------------------------------------
      SUBROUTINE DDMATVAR( X, M, LCM, LCN, ICTXT, VAR)
      ! Inputs
      INTEGER             M, LCM, LCN
      DOUBLE PRECISION    X( LCM, LCN )
      ! Output
      DOUBLE PRECISION    VAR( LCN )
      ! Local
      INTEGER             I, ICTXT
      DOUBLE PRECISION    SCL,
     $                    MN( LCN )
      ! Subroutines
      EXTERNAL            DGSUM2D, DDMATMN
      
      ! --------------------------
      ! Variance calculator
      ! --------------------------
      
      ! mean
      CALL DDMATMN( X, M, LCM, LCN, DESCX, MN)
      
      ! variance
      DO I = 1, LCN
        VAR(I) = 0
        DO J = 1, LCM, 1
          VAR(I) = VAR(I) + ( X(J,I) ** 2 ) / (M-1)
        END DO
      END DO
      
      CALL DGSUM2D( ICTXT, 'Column', ' ', LCN, 1, VAR, 1, -1, -1 )
      
      SCL = 1.0D0 * M/(M-1)
!      VAR = VAR - SCL * (MN ** 2) ! loses precision?
      DO I = 1, LCN, 1
        VAR(I) = VAR(I) - SCL * (MN(I)**2)
      END DO
      
      RETURN
      END









