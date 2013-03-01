!! Copyright 2013, Schmidt

! colmeans
      SUBROUTINE PDCLMN(X, DESCX, MN)
      IMPLICIT NONE
      ! IN/OUT
      INTEGER             DESCX( 9 )
      DOUBLE PRECISION    X(DESCX(9), *), MN(*)
      ! Local
      INTEGER             M, N, I, J, LDM(2), BLACS(5)
      DOUBLE PRECISION    DENOM
      !External
      EXTERNAL            PDIMS, DGSUM2D
      
      
      ! Get local and proc grid info
      CALL PDIMS(DESCX, LDM, BLACS)
      M = LDM(1)
      N = LDM(2)
      
      DENOM = DBLE(DESCX(3))
      
      DO J = 1, N
        MN(J) = 0
        DO I = 1, M
          MN(J) = MN(J) + X(I,J)
        END DO
        MN(J) = MN(J) / DENOM
      END DO
      
      CALL DGSUM2D(DESCX(2), 'Column', ' ', N, 1, MN, N, -1, -1)
      
      RETURN
      END


! colVar
      SUBROUTINE PDCLVAR(X, DESCX, VAR)
      IMPLICIT NONE
      ! IN/OUT
      INTEGER             DESCX( 9 )
      DOUBLE PRECISION    X(DESCX(9), *), VAR(*)
      ! Local
      INTEGER             M, N, I, J, LDM(2), BLACS(5), ALLOCERR
      DOUBLE PRECISION    SCL, DENOM, WORK(2)
      DOUBLE PRECISION, ALLOCATABLE :: MN(:)
      !External
      EXTERNAL            PDIMS, PDCLMN, DGSUM2D
      DOUBLE PRECISION    PDVAR
      
      
      ! Get local and proc grid info
      CALL PDIMS(DESCX, LDM, BLACS)
      M = LDM(1)
      N = LDM(2)
      
      DENOM = DBLE(DESCX(3))
      
      ALLOCATE(MN(N), STAT=ALLOCERR)
      IF (ALLOCERR.NE.0) RETURN ! STOP "Out of memory"
      
      CALL PDCLMN(X, DESCX, MN)
      
      DO J = 1, N
        WORK(1:2) = 0.0D0
        DO I = 1, M
          WORK(1) = WORK(1) + (X(I,J) - MN(J))**2
          WORK(2) = WORK(2) + X(I,J) - MN(J)
        END DO
        CALL DGSUM2D(DESCX(2), 'Col', ' ', 2, 1, WORK, 2, -1, -1)
        VAR(J) = (WORK(1) - (WORK(2)**2)/DENOM) / (DENOM - 1)
      END DO
      
      DEALLOCATE(MN)
      
      RETURN
      END


