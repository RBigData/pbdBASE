!! Copyright 2013, Schmidt

! colmeans
      SUBROUTINE PDCLMN(X, DESCX, MN)
      IMPLICIT NONE
      ! IN/OUT
      INTEGER             DESCX( 9 )
      DOUBLE PRECISION    X(DESCX(9), *), MN(*)
      ! Local
      INTEGER             M, N, I, J, LDM(2), BLACS(5), NN
      DOUBLE PRECISION    DENOM
      !External
      EXTERNAL            PDIMS, DGSUM2D
      
      
      ! Get local and proc grid info
      CALL PDIMS(DESCX, LDM, BLACS)
      M = LDM(1)
      N = LDM(2)
      
      DENOM = DBLE(DESCX(3))
      
      DO J = 1, N
        MN(J) = 0.0D0
        DO I = 1, M
          MN(J) = MN(J) + X(I,J)
        END DO
        MN(J) = MN(J) / DENOM
      END DO
      
      NN = N
      CALL IGAMX2D(DESCX(2), 'Col', ' ', 1, 1, NN,1,-1,-1,-1,-1,-1)
      
      CALL DGSUM2D(DESCX(2), 'Col', ' ', NN, 1, MN, NN, -1, -1)
      
      RETURN
      END


! colVar
      SUBROUTINE PDCLVAR(X, DESCX, VAR)
      IMPLICIT NONE
      ! IN/OUT
      INTEGER             DESCX( 9 )
      DOUBLE PRECISION    X(DESCX(9), *), VAR(*)
      ! Local
      INTEGER             M, N, I, J, LDM(2), BLACS(5), ALLOCERR, NN
      DOUBLE PRECISION    SCL, DENOM
      DOUBLE PRECISION, ALLOCATABLE :: MN(:), WORK(:,:)
      ! Parameter
      DOUBLE PRECISION    ZERO
      PARAMETER ( ZERO = 0.0D0 )
      !External
      EXTERNAL            PDIMS, PDCLMN, IGAMX2D, DGSUM2D
      DOUBLE PRECISION    PDVAR
      
      
      ! Get local and proc grid info
      CALL PDIMS(DESCX, LDM, BLACS)
      M = LDM(1)
      N = LDM(2)
      
      DENOM = DBLE(DESCX(3))
      
      NN = N
      CALL IGAMX2D(DESCX(2), 'Col', ' ', 1, 1, NN,1,-1,-1,-1,-1,-1)
      
      ALLOCERR = 0
      ALLOCATE(MN(NN), STAT=ALLOCERR)
      IF (ALLOCERR.NE.0) RETURN ! STOP "Out of memory"
      ALLOCATE(WORK(2, NN), STAT=ALLOCERR)
      IF (ALLOCERR.NE.0) RETURN ! STOP "Out of memory"
      
      MN(1:N) = ZERO
      CALL PDCLMN(X, DESCX, MN)
      
      WORK(1:2, 1:NN) = ZERO
      
      DO J = 1, N
        DO I = 1, M
          WORK(1, J) = WORK(1, J) + (X(I,J) - MN(J))**2
          WORK(2, J) = WORK(2, J) + X(I,J) - MN(J)
        END DO
      END DO
      
      CALL DGSUM2D(DESCX(2), 'Col', ' ', 2*NN, 1, WORK, 2*NN, -1, -1)
      
      DO J = 1, N
        VAR(J) = (WORK(1, J) - (WORK(2, J)**2)/DENOM) / (DENOM - 1)
      END DO
      
      DEALLOCATE(MN)
      DEALLOCATE(WORK)
      
      RETURN
      END


