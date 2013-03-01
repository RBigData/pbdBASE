!! Copyright 2013, Schmidt


!     PDLANGE subroutine wrapper for use with F77_CALL in C
!     Also handles allocation of work vector
      SUBROUTINE MATNORM(VALUE, NORM, M, N, A, IA, JA, DESCA)
      ! #######################################################################
      ! In/Out
      CHARACTER           NORM
      INTEGER             M, N, IA, JA, DESCA( 9 )
      DOUBLE PRECISION    VALUE, A( * )
      ! Local
      INTEGER             LWORK, BL, IOFFA, IAMAR, NPROW, NPCOL, MYPROW, 
     $                    MYPCOL
      ! Dynamic
      DOUBLE PRECISION, ALLOCATABLE :: WORK(:)
      ! Parameters
      INTEGER             I0
      PARAMETER ( I0 = 0 )
      ! Functions
      DOUBLE PRECISION    PDLANGE
      INTEGER             INDXG2P, NUMROC
      ! #######################################################################
      
      IF (NORM.EQ."M" .OR. NORM.EQ."F") THEN
        LWORK = I0
      ELSE 
        CALL BLACS_GRIDINFO(DESCA(2), NPROW, NPCOL, MYPROW, MYPCOL)
        
        IF (NORM.EQ."O") THEN
          BL = DESCA(5)
          IOFFA = MOD(JA-1, BL)
          
          IAMAR = INDXG2P(IA, BL, MYPCOL, I0, NPCOL)
          
          LWORK = NUMROC(N+IOFFA, BL, MYPCOL, I0, NPCOL)
        ELSE IF (NORM.EQ."I") THEN
          BL = DESCA(6)
          IOFFA = MOD(IA-1, BL)
          
          IAMAR = INDXG2P(IA, BL, MYPROW, I0, NPROW)
          
          LWORK = NUMROC(M+IOFFA, BL, MYPROW, I0, NPROW)
        END IF
      END IF
      
      ALLOCATE (WORK(LWORK))
      VALUE = PDLANGE(NORM, M, N, A, IA, JA, DESCA, WORK)
      DEALLOCATE (WORK)
      
      RETURN
      END



!     Condition number estimator for general matrix
!       Step 1:  get matrix norm of A
!       Step 2:  factor A=LU
!       Step 3:  Call pdgecon
      SUBROUTINE CONDNUM(NORM, M, N, A, IA, JA, DESCA, RCOND, INFO)
      ! #######################################################################
      ! In/Out
      CHARACTER           NORM
      INTEGER             M, N, LDIM1, IA, JA, DESCA( 9 ), INFO
      DOUBLE PRECISION    RCOND, A( * )
      ! Local
      INTEGER             LWORK, LIWORK, LIPIV
      DOUBLE PRECISION    ANORM, TMP
      ! Parameters
      INTEGER             IN1
      PARAMETER ( IN1 = -1 )
      ! Dynamic
      DOUBLE PRECISION, ALLOCATABLE :: WORK(:)
      INTEGER, ALLOCATABLE :: IWORK(:), IPIV(:)
      ! Functions
      EXTERNAL            MATNORM, PDGETRF, PDGECON
      ! #######################################################################
      
      ! Step 1:  get matrix norm of A
      CALL MATNORM(ANORM, NORM, N, N, A, IA, JA, DESCA)
      
      ! Step 2:  factor A=LU
      LIPIV = DESCA(9) + DESCA(5) ! LOCr(m_a)+ mb_a
      ALLOCATE (IPIV(LIPIV))
      
      CALL PDGETRF(M, N, A, IA, JA, DESCA, IPIV, INFO)
      
      IF (INFO.NE.0) THEN
        RETURN
      END IF
      
      ! Step 3:  Call pdgecon
      CALL PDGECON(NORM, N, A, IA, JA, DESCA, ANORM, RCOND, 
     $             TMP, IN1, LIWORK, IN1, INFO)
      
      LWORK = INT(TMP)
      ALLOCATE (WORK(LWORK))
      ALLOCATE (IWORK(LIWORK))
      
      CALL PDGECON(NORM, N, A, IA, JA, DESCA, ANORM, RCOND, 
     $             WORK, LWORK, IWORK, LIWORK, INFO)
      
      DEALLOCATE (IPIV)
      DEALLOCATE (IWORK)
      DEALLOCATE (WORK)
      
      RETURN
      END


