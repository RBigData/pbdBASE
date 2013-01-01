! NOTE FIXME
! Generate R message:
! Error: number of cluster centres must lie between 1 and nrow(x)
! so assume that 1<k<nrow(X)=M

! mu = [mu1; mu2; ... ; mun]
! X = mxn
! mu = kxn
! Z = m-length vector

! Lloyd's Algorithm for K Means Clustering on block cyclic data
! K = number of clusters
! X = Data matrix to be clustered
! DESCX = ScaLAPACK descriptor array for X
! MU = Cluster centers, global, owned by all
! DESCMU = Descriptor array for MU
! Z = classification vector
! IMAX = Max number of iterations
      SUBROUTINE BCKMNCL1(K, X, DESCX, MU, DESCMU, Z, IMAX, LDIM1, LDIM2, 
     $                    INIMHD)
      ! #######################################################################
      ! In/Out
      INTEGER             K, IMAX, LDIM1, LDIM2, INIMHD
      INTEGER             DESCX( 9 ), 
     $                    DESCMU( 9 ), 
     $                    Z(LDIM1)
      DOUBLE PRECISION    X(LDIM1, LDIM2), 
     $                    MU(, )
      ! Local
      INTEGER             I, J, ITER, TESET, CHND,
     $                    NPROW, NPCOL, MYPROW, MYPCOL,
     $                    INFO,
     $                    DESCROW( 9 )
      DOUBLE PRECISION    ZOLD, TMP, NRM, LSTBST
      ! Constants
      INTEGER             IONE, INONE
      PARAMETER ( IONE = 1, INONE = -1 )
      DOUBLE PRECISION    ZERO, ONE, NEGONE
      PARAMETER ( ZERO = 0.0D0, ONE = 1.0D0, NEGONE = -1.0D0 )
      ! External Functions
      EXTERNAL           PDGEADD, PDNRM2, DGAMX2D, BLACS_GRIDINFO, DESCINIT
      ! #######################################################################
      
      ! Descriptor array for a row of the global distributed matrix
      CALL DESCINIT( DESCROW, DESCX(4), DESCX(5), DESCX(6), 1, 1, DESCX(2), 
     $               LDIM2, INFO )
      
      ! Initialize
      IF ( INIMHD.EQ.1 ) THEN
        CALL BCKMNIT1()
      ELSE IF ( INIMHD.EQ.2 ) THEN
        CALL BCKMNIT2()
      END IF
      
      ! Main
      DO 10, ITER = 2, IMAX
        ! Assign grouping by row x_n
        ! x_n in k whenever 
        DO 100, I = 1, DESCX(3) ! for each row of X
          TMP = X(I,:)
          TEST = 0
          DO 200, J = 1, K ! for each of the cluster centers
            TMPMU = MU(I,:)
            ! Compute euclidean norm (x(i,:)-mu(j,:))^T * (x(i,:)-mu(j,:))
            CALL PDGEADD('N', DESCX(4), 1, NEGONE, TMPMU, IONE, IONE, DESCROW, 
     $                    ONE, TMP, IONE, IONE, DESCROW)
            
            CALl PDNRM2(DESCX(4), NRM, TMP, IONE, IONE, DESCROW, IONE)
            
            IF ( I.EQ.1 ) THEN
              LSTBST = NRM
              Z(I) = 1
            ELSE
              IF ( LSTBST.GT.NRM ) THEN
                Z(I) = I
              END IF
            END IF
  200 CONTINUE
  100 CONTINUE
        
        ! Check if assignment has changed, act accordingly
        DO 20, I = 1, LDIM1
          ZOLD(I) = Z(I)
          IF ( Z(I).EQ.ZOLD(I) ) THEN
            CHND = 1
            GOTO 20
          END IF
          CHND = 0
   20 CONTINUE
        DO 21, J = I, LDIM1
          ZOLD(I) = Z(I)
   21 CONTINUE
        
        CALL DGAMX2D( DESCX(2), 'A', ' ', IONE, IONE, CHND, IONE, INONE, INONE, 
     $                INONE, INONE, INONE )
        
        IF ( CHND.EQ.1 ) THEN
          GOTO 10
        ELSE
          CALL BCKMNCUD(K, X, DESCX, MU, DESCMU, Z, M)
        END IF
        
   10 CONTINUE
      
      RETURN
      END



! Block-Cyclic K Means Centroid Updater
      SUBROUTINE BCKMNCUD(K, X, DESCX, MU, DESCMU, Z, M)
      ! #######################################################################
      ! In/Out
      INTEGER             K, M
      INTEGER             DESCX( 9 ), 
     $                    DESCMU( 9 ), 
     $                    Z( M )
      DOUBLE PRECISION    X(, ),
     $                    MU(, )
      ! Local
      INTEGER             COUNTS( K )
      ! #######################################################################
      
      ! MU_i = sum_j( Z(i)*X(j, ) ) / COUNTS(i)
      
      
      RETURN
      END



! Initialize Centroids --- Forgy's method
! Randomly choose K points as the centroids
      SUBROUTINE BCKMNIT1()
      ! #######################################################################
      ! In/Out
      INTEGER             K, M
      ! Local
      INTEGER             TIME,
     $                    ICTXT, NPROW, NPCOL, MYPROW, MYPCOL,
     $                    ROW( M ), MUROW( K )
      INTEGER, DIMENSION( M ) :: ROW = (/ (I, I=0, M-1) /)
      DOUBLE PRECISION    R, TMP
      ! External functions
      EXTERNAL            DGSUM2D, BLACS_GRIDINFO
      ! #######################################################################
      
      ! Get BLACS context information
      ICTXT = DESCX(2)
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYPROW, MYPCOL ) 
      
      IF ( MYPROW.EQ.0 .AND. MYPCOL.EQ.0 ) THEN
        CALL SYSTEM_CLOCK( TIME )
        CALL RANDOM_SEED( TIME )
        DO 10, I = 1, K
          CALL RANDOM_NUMBER( R )
          NUM = INT( M*R ) + 1
          TMP = ROW(NUM)
          ROW(NUM) = ROW(I)
          ROW(I) = TMP
          MUROW(I) = TMP
   10 CONTINUE
      
      CALL DGSUM2D( ICTXT, 'A', ' ', M, 1, MUROW, LDA, RDEST, CDEST )
      
      
      
      CALL DGSUM2D( ICTXT, 'A', ' ', K, N, MU, LDA, RDEST, CDEST )
      
      
      RETURN
      END






! Initialize Centroids --- Random partition
! Randomly assign each observation to a group
      SUBROUTINE BCKMNIT2()
      ! #######################################################################
      ! In/Out
      
      ! #######################################################################
      
      
      RETURN
      END




