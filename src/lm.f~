!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     RPDGEQPF:  Modified PDGEQPF to use the 'limited 
!     pivoting strategy' from R's dqrdc2
!     
!     Here, instead of pivoting on the column with greatest
!     norm, we instead pivot on the columns in order
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      SUBROUTINE RPDGEQPF( TOL, M, N, A, IA, JA, DESCA, IPIV, TAU, 
     $                    WORK, LWORK, RANK, INFO )
*
*     Originally:
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     March 14, 2000
*
*     Modified date October 9, 2012
*
*     .. Scalar Arguments ..
      INTEGER            IA, JA, INFO, LWORK, M, N, RANK
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), IPIV( * )
      DOUBLE PRECISION   A( * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PDGEQPF computes a QR factorization with column pivoting of a
*  M-by-N distributed matrix sub( A ) = A(IA:IA+M-1,JA:JA+N-1):
*
*                         sub( A ) * P = Q * R.
*
*  Notes
*  =====
*
*  Each global data object is described by an associated description
*  vector.  This vector stores the information required to establish
*  the mapping between an object element and its corresponding process
*  and memory location.
*
*  Let A be a generic term for any 2D block cyclicly distributed array.
*  Such a global array has an associated description vector DESCA.
*  In the following comments, the character _ should be read as
*  "of the global array".
*
*  NOTATION        STORED IN      EXPLANATION
*  --------------- -------------- --------------------------------------
*  DTYPE_A(global) DESCA( DTYPE_ )The descriptor type.  In this case,
*                                 DTYPE_A = 1.
*  CTXT_A (global) DESCA( CTXT_ ) The BLACS context handle, indicating
*                                 the BLACS process grid A is distribu-
*                                 ted over. The context itself is glo-
*                                 bal, but the handle (the integer
*                                 value) may vary.
*  M_A    (global) DESCA( M_ )    The number of rows in the global
*                                 array A.
*  N_A    (global) DESCA( N_ )    The number of columns in the global
*                                 array A.
*  MB_A   (global) DESCA( MB_ )   The blocking factor used to distribute
*                                 the rows of the array.
*  NB_A   (global) DESCA( NB_ )   The blocking factor used to distribute
*                                 the columns of the array.
*  RSRC_A (global) DESCA( RSRC_ ) The process row over which the first
*                                 row of the array A is distributed.
*  CSRC_A (global) DESCA( CSRC_ ) The process column over which the
*                                 first column of the array A is
*                                 distributed.
*  LLD_A  (local)  DESCA( LLD_ )  The leading dimension of the local
*                                 array.  LLD_A >= MAX(1,LOCr(M_A)).
*
*  Let K be the number of rows or columns of a distributed matrix,
*  and assume that its process grid has dimension p x q.
*  LOCr( K ) denotes the number of elements of K that a process
*  would receive if K were distributed over the p processes of its
*  process column.
*  Similarly, LOCc( K ) denotes the number of elements of K that a
*  process would receive if K were distributed over the q processes of
*  its process row.
*  The values of LOCr() and LOCc() may be determined via a call to the
*  ScaLAPACK tool function, NUMROC:
*          LOCr( M ) = NUMROC( M, MB_A, MYROW, RSRC_A, NPROW ),
*          LOCc( N ) = NUMROC( N, NB_A, MYCOL, CSRC_A, NPCOL ).
*  An upper bound for these quantities may be computed by:
*          LOCr( M ) <= ceil( ceil(M/MB_A)/NPROW )*MB_A
*          LOCc( N ) <= ceil( ceil(N/NB_A)/NPCOL )*NB_A
*
*  Arguments
*  =========
*
*  M       (global input) INTEGER
*          The number of rows to be operated on, i.e. the number of rows
*          of the distributed submatrix sub( A ). M >= 0.
*
*  N       (global input) INTEGER
*          The number of columns to be operated on, i.e. the number of
*          columns of the distributed submatrix sub( A ). N >= 0.
*
*  A       (local input/local output) DOUBLE PRECISION pointer into the
*          local memory to an array of dimension (LLD_A, LOCc(JA+N-1)).
*          On entry, the local pieces of the M-by-N distributed matrix
*          sub( A ) which is to be factored. On exit, the elements on
*          and above the diagonal of sub( A ) contain the min(M,N) by N
*          upper trapezoidal matrix R (R is upper triangular if M >= N);
*          the elements below the diagonal, with the array TAU, repre-
*          sent the orthogonal matrix Q as a product of elementary
*          reflectors (see Further Details).
*
*  IA      (global input) INTEGER
*          The row index in the global array A indicating the first
*          row of sub( A ).
*
*  JA      (global input) INTEGER
*          The column index in the global array A indicating the
*          first column of sub( A ).
*
*  DESCA   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix A.
*
*  IPIV    (local output) INTEGER array, dimension LOCc(JA+N-1).
*          On exit, if IPIV(I) = K, the local i-th column of sub( A )*P
*          was the global K-th column of sub( A ). IPIV is tied to the
*          distributed matrix A.
*
*  TAU     (local output) DOUBLE PRECISION array, dimension
*          LOCc(JA+MIN(M,N)-1). This array contains the scalar factors
*          TAU of the elementary reflectors. TAU is tied to the
*          distributed matrix A.
*
*  WORK    (local workspace/local output) DOUBLE PRECISION array,
*                                                   dimension (LWORK)
*          On exit, WORK(1) returns the minimal and optimal LWORK.
*
*  LWORK   (local or global input) INTEGER
*          The dimension of the array WORK.
*          LWORK is local input and must be at least
*          LWORK >= MAX(3,Mp0 + Nq0) + LOCc(JA+N-1)+Nq0.
*
*          IROFF = MOD( IA-1, MB_A ), ICOFF = MOD( JA-1, NB_A ),
*          IAROW = INDXG2P( IA, MB_A, MYROW, RSRC_A, NPROW ),
*          IACOL = INDXG2P( JA, NB_A, MYCOL, CSRC_A, NPCOL ),
*          Mp0   = NUMROC( M+IROFF, MB_A, MYROW, IAROW, NPROW ),
*          Nq0   = NUMROC( N+ICOFF, NB_A, MYCOL, IACOL, NPCOL ),
*          LOCc(JA+N-1) = NUMROC( JA+N-1, NB_A, MYCOL, CSRC_A, NPCOL )
*
*          and NUMROC, INDXG2P are ScaLAPACK tool functions;
*          MYROW, MYCOL, NPROW and NPCOL can be determined by calling
*          the subroutine BLACS_GRIDINFO.
*
*          If LWORK = -1, then LWORK is global input and a workspace
*          query is assumed; the routine only calculates the minimum
*          and optimal size for all work arrays. Each of these
*          values is returned in the first entry of the corresponding
*          work array, and no error message is issued by PXERBLA.
*
*
*  INFO    (global output) INTEGER
*          = 0:  successful exit
*          < 0:  If the i-th argument is an array and the j-entry had
*                an illegal value, then INFO = -(i*100+j), if the i-th
*                argument is a scalar and had an illegal value, then
*                INFO = -i.
*
*  Further Details
*  ===============
*
*  The matrix Q is represented as a product of elementary reflectors
*
*     Q = H(1) H(2) . . . H(n)
*
*  Each H(i) has the form
*
*     H = I - tau * v * v'
*
*  where tau is a real scalar, and v is a real vector with v(1:i-1) = 0
*  and v(i) = 1; v(i+1:m) is stored on exit in A(ia+i-1:ia+m-1,ja+i-1).
*
*  The matrix P is represented in jpvt as follows: If
*     jpvt(j) = i
*  then the jth column of P is the ith canonical unit vector.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      DOUBLE PRECISION   ONE, ZERO, TOL
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IACOL, IAROW, ICOFF, ICTXT, ICURROW,
     $                   ICURCOL, II, IIA, IOFFA, IPN, IPCOL, IPW,
     $                   IROFF, ITEMP, J, JB, JJ, JJA, JJPVT, JN, KB,
     $                   K, KK, KSTART, KSTEP, LDA, LL, LWMIN, MN, MP,
     $                   MYCOL, MYROW, NPCOL, NPROW, NQ, NQ0, PVT,
     $                   IJ
      DOUBLE PRECISION   AJJ, ALPHA, TEMP, TEMP2, TEMP3
*     ..
*     .. Local Arrays ..
      INTEGER            DESCN( DLEN_ ), IDUM1( 1 ), IDUM2( 1 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, DCOPY, DESCSET,
     $                   DGEBR2D, DGEBS2D, DGERV2D,
     $                   DGESD2D, DLARFG, DSWAP, IGERV2D,
     $                   IGESD2D, INFOG1L, INFOG2L, PCHK1MAT, PDAMAX,
     $                   PDELSET, PDLARF, PDLARFG, PDNRM2,
     $                   PXERBLA, PDELGET
*     ..
*     .. External Functions ..
      INTEGER            ICEIL, INDXG2P, NUMROC
      EXTERNAL           ICEIL, INDXG2P, NUMROC
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, IDINT, MAX, MIN, MOD, SQRT
*     ..
*     .. Executable Statements ..
*
*     Get grid parameters
*
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
*     Test the input parameters
*
      INFO = 0
      IF( NPROW.EQ.-1 ) THEN
         INFO = -(600+CTXT_)
      ELSE
         CALL CHK1MAT( M, 1, N, 2, IA, JA, DESCA, 6, INFO )
         IF( INFO.EQ.0 ) THEN
            IROFF = MOD( IA-1, DESCA( MB_ ) )
            ICOFF = MOD( JA-1, DESCA( NB_ ) )
            IAROW = INDXG2P( IA, DESCA( MB_ ), MYROW, DESCA( RSRC_ ),
     $                       NPROW )
            IACOL = INDXG2P( JA, DESCA( NB_ ), MYCOL, DESCA( CSRC_ ),
     $                       NPCOL )
            MP = NUMROC( M+IROFF, DESCA( MB_ ), MYROW, IAROW, NPROW )
            NQ = NUMROC( N+ICOFF, DESCA( NB_ ), MYCOL, IACOL, NPCOL )
            NQ0 = NUMROC( JA+N-1, DESCA( NB_ ), MYCOL, DESCA( CSRC_ ),
     $                    NPCOL )
            LWMIN = MAX( 3, MP + NQ ) + NQ0 + NQ
*
            WORK( 1 ) = DBLE( LWMIN )
            LQUERY = ( LWORK.EQ.-1 )
            IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY )
     $         INFO = -10
         END IF
         IF( LWORK.EQ.-1 ) THEN
            IDUM1( 1 ) = -1
         ELSE
            IDUM1( 1 ) = 1
         END IF
         IDUM2( 1 ) = 10
         CALL PCHK1MAT( M, 1, N, 2, IA, JA, DESCA, 6, 1, IDUM1, IDUM2,
     $                  INFO )
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PDGEQPF', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
      CALL INFOG2L( IA, JA, DESCA, NPROW, NPCOL, MYROW, MYCOL, IIA, JJA,
     $              IAROW, IACOL )
      IF( MYROW.EQ.IAROW )
     $   MP = MP - IROFF
      IF( MYCOL.EQ.IACOL )
     $   NQ = NQ - ICOFF
      MN = MIN( M, N )
*
*     Initialize the array of pivots
*
      LDA = DESCA( LLD_ )
      JN = MIN( ICEIL( JA, DESCA( NB_ ) ) * DESCA( NB_ ), JA+N-1 )
      KSTEP  = NPCOL * DESCA( NB_ )
*
      IF( MYCOL.EQ.IACOL ) THEN
*
*        Handle first block separately
*
         JB = JN - JA + 1
         DO 10 LL = JJA, JJA+JB-1
            IPIV( LL ) = JA + LL - JJA
   10    CONTINUE
         KSTART = JN + KSTEP - DESCA( NB_ )
*
*        Loop over remaining block of columns
*
         DO 30 KK = JJA+JB, JJA+NQ-1, DESCA( NB_ )
            KB = MIN( JJA+NQ-KK, DESCA( NB_ ) )
            DO 20 LL = KK, KK+KB-1
               IPIV( LL ) = KSTART+LL-KK+1
   20       CONTINUE
            KSTART = KSTART + KSTEP
   30    CONTINUE
      ELSE
         KSTART = JN + ( MOD( MYCOL-IACOL+NPCOL, NPCOL )-1 )*
     $                        DESCA( NB_ )
         DO 50 KK = JJA, JJA+NQ-1, DESCA( NB_ )
            KB = MIN( JJA+NQ-KK, DESCA( NB_ ) )
            DO 40 LL = KK, KK+KB-1
               IPIV( LL ) = KSTART+LL-KK+1
   40       CONTINUE
            KSTART = KSTART + KSTEP
   50    CONTINUE
      END IF
*
*     Initialize partial column norms, handle first block separately
*
      CALL DESCSET( DESCN, 1, DESCA( N_ ), 1, DESCA( NB_ ), MYROW,
     $              DESCA( CSRC_ ), ICTXT, 1 )
*
      IPN = 1
      IPW = IPN + NQ0 + NQ
      JJ = IPN + JJA - 1
      IF( MYCOL.EQ.IACOL ) THEN
         DO 60 KK = 0, JB-1
            CALL PDNRM2( M, WORK( JJ+KK ), A, IA, JA+KK, DESCA, 1 )
            WORK( NQ+JJ+KK ) = WORK( JJ+KK )
   60    CONTINUE
         JJ = JJ + JB
      END IF
      ICURCOL = MOD( IACOL+1, NPCOL )
*
*     Loop over the remaining blocks of columns
*
      DO 80 J = JN+1, JA+N-1, DESCA( NB_ )
         JB = MIN( JA+N-J, DESCA( NB_ ) )
*
         IF( MYCOL.EQ.ICURCOL ) THEN
            DO 70 KK = 0, JB-1
               CALL PDNRM2( M, WORK( JJ+KK ), A, IA, J+KK, DESCA, 1 )
               WORK( NQ+JJ+KK ) = WORK( JJ+KK )
   70       CONTINUE
            JJ = JJ + JB
         END IF
         ICURCOL = MOD( ICURCOL+1, NPCOL )
   80 CONTINUE
*
*     Compute factorization
*
      TEMP3 = 0
      DO 120 J = JA, JA+MN-1
         I = IA + J - JA
*
         CALL INFOG1L( J, DESCA( NB_ ), NPCOL, MYCOL, DESCA( CSRC_ ),
     $                 JJ, ICURCOL )
         K = JA + N - J 
         IF( K.GT.1 ) THEN
!            CALL PDAMAX( K, TEMP, PVT, WORK( IPN ), 1, J, DESCN,
!     $                   DESCN( M_ ) )
            DO 11 IJ = J, N, 1
              CALL PDELGET( 'A', 'T', TEMP, WORK(IPN), 1, IJ, DESCN )
              IF( TEMP.GT.TOL ) THEN
                PVT = IJ
                TEMP3 = IJ - J ! for tracking "N - RANK"
                EXIT
              END IF
   11 CONTINUE
         ELSE
            PVT = J
         END IF
         IF( J.NE.PVT ) THEN
            CALL INFOG1L( PVT, DESCA( NB_ ), NPCOL, MYCOL,
     $                    DESCA( CSRC_ ), JJPVT, IPCOL )
            IF( ICURCOL.EQ.IPCOL ) THEN
               IF( MYCOL.EQ.ICURCOL ) THEN
                  CALL DSWAP( MP, A( IIA+(JJ-1)*LDA ), 1,
     $                        A( IIA+(JJPVT-1)*LDA ), 1 )
                  ITEMP = IPIV( JJPVT )
                  IPIV( JJPVT ) = IPIV( JJ )
                  IPIV( JJ ) = ITEMP
                  WORK( IPN+JJPVT-1 ) = WORK( IPN+JJ-1 )
                  WORK( IPN+NQ+JJPVT-1 ) = WORK( IPN+NQ+JJ-1 )
               END IF
            ELSE
               IF( MYCOL.EQ.ICURCOL ) THEN
*
                  CALL DGESD2D( ICTXT, MP, 1, A( IIA+(JJ-1)*LDA ), LDA,
     $                          MYROW, IPCOL )
                  WORK( IPW )   = DBLE( IPIV( JJ ) )
                  WORK( IPW+1 ) = WORK( IPN + JJ - 1 )
                  WORK( IPW+2 ) = WORK( IPN + NQ + JJ - 1 )
                  CALL DGESD2D( ICTXT, 3, 1, WORK( IPW ), 3, MYROW,
     $                          IPCOL )
*
                  CALL DGERV2D( ICTXT, MP, 1, A( IIA+(JJ-1)*LDA ), LDA,
     $                          MYROW, IPCOL )
                  CALL IGERV2D( ICTXT, 1, 1, IPIV( JJ ), 1, MYROW,
     $                          IPCOL )
*
               ELSE IF( MYCOL.EQ.IPCOL ) THEN
*
                  CALL DGESD2D( ICTXT, MP, 1, A( IIA+(JJPVT-1)*LDA ),
     $                          LDA, MYROW, ICURCOL )
                  CALL IGESD2D( ICTXT, 1, 1, IPIV( JJPVT ), 1, MYROW,
     $                          ICURCOL )
*
                  CALL DGERV2D( ICTXT, MP, 1, A( IIA+(JJPVT-1)*LDA ),
     $                          LDA, MYROW, ICURCOL )
                  CALL DGERV2D( ICTXT, 3, 1, WORK( IPW ), 3, MYROW,
     $                          ICURCOL )
                  IPIV( JJPVT ) = IDINT( WORK( IPW ) )
                  WORK( IPN+JJPVT-1 ) = WORK( IPW+1 )
                  WORK( IPN+NQ+JJPVT-1 ) = WORK( IPW+2 )
*
               END IF
*
            END IF
*
         END IF
*
*        Generate elementary reflector H(i)
*
         CALL INFOG1L( I, DESCA( MB_ ), NPROW, MYROW, DESCA( RSRC_ ),
     $                 II, ICURROW )
         IF( DESCA( M_ ).EQ.1 ) THEN
            IF( MYROW.EQ.ICURROW ) THEN
               IF( MYCOL.EQ.ICURCOL ) THEN
                  IOFFA = II+(JJ-1)*DESCA( LLD_ )
                  AJJ = A( IOFFA )
                  CALL DLARFG( 1, AJJ, A( IOFFA ), 1, TAU( JJ ) )
                  IF( N.GT.1 ) THEN
                     ALPHA = ONE - TAU( JJ )
                     CALL DGEBS2D( ICTXT, 'Rowwise', ' ', 1, 1, ALPHA,
     $                                  1 )
                     CALL DSCAL( NQ-JJ, ALPHA, A( IOFFA+DESCA( LLD_ ) ),
     $                           DESCA( LLD_ ) )
                  END IF
                  CALL DGEBS2D( ICTXT, 'Columnwise', ' ', 1, 1,
     $                          TAU( JJ ), 1 )
                  A( IOFFA ) = AJJ
               ELSE
                  IF( N.GT.1 ) THEN
                     CALL DGEBR2D( ICTXT, 'Rowwise', ' ', 1, 1, ALPHA,
     $                             1, ICURROW, ICURCOL )
                     CALL DSCAL( NQ-JJ+1, ALPHA, A( I ), DESCA( LLD_ ) )
                  END IF
               END IF
            ELSE IF( MYCOL.EQ.ICURCOL ) THEN
               CALL DGEBR2D( ICTXT, 'Columnwise', ' ', 1, 1, TAU( JJ ),
     $                       1, ICURROW, ICURCOL )
            END IF
*
         ELSE
*
            CALL PDLARFG( M-J+JA, AJJ, I, J, A, MIN( I+1, IA+M-1 ), J,
     $                    DESCA, 1, TAU )
            IF( J.LT.JA+N-1 ) THEN
*
*              Apply H(i) to A(ia+j-ja:ia+m-1,j+1:ja+n-1) from the left
*
               CALL PDELSET( A, I, J, DESCA, ONE )
               CALL PDLARF( 'Left', M-J+JA, JA+N-1-J, A, I, J, DESCA,
     $                      1, TAU, A, I, J+1, DESCA, WORK( IPW ) )
            END IF
            CALL PDELSET( A, I, J, DESCA, AJJ )
*
         END IF
*
*        Update partial columns norms
*
         IF( MYCOL.EQ.ICURCOL )
     $      JJ = JJ + 1
         IF( MOD( J, DESCA( NB_ ) ).EQ.0 )
     $      ICURCOL = MOD( ICURCOL+1, NPCOL )
         IF( (JJA+NQ-JJ).GT.0 ) THEN
            IF( MYROW.EQ.ICURROW ) THEN
               CALL DGEBS2D( ICTXT, 'Columnwise', ' ', 1, JJA+NQ-JJ,
     $                       A( II+( MIN( JJA+NQ-1, JJ )-1 )*LDA ),
     $                       LDA )
               CALL DCOPY( JJA+NQ-JJ, A( II+( MIN( JJA+NQ-1, JJ )
     $                     -1)*LDA ), LDA, WORK( IPW+MIN( JJA+NQ-1,
     $                    JJ )-1 ), 1 )
            ELSE
               CALL DGEBR2D( ICTXT, 'Columnwise', ' ', JJA+NQ-JJ, 1,
     $                       WORK( IPW+MIN( JJA+NQ-1, JJ )-1 ),
     $                       MAX( 1, NQ ), ICURROW, MYCOL )
            END IF
         END IF
*
         JN = MIN( ICEIL( J+1, DESCA( NB_ ) ) * DESCA( NB_ ),
     $                    JA + N - 1 )
         IF( MYCOL.EQ.ICURCOL ) THEN
            DO 90 LL = JJ-1, JJ + JN - J - 2
               IF( WORK( IPN+LL ).NE.ZERO ) THEN
                  TEMP = ONE-( ABS( WORK( IPW+LL ) ) /
     $                         WORK( IPN+LL ) )**2
                  TEMP = MAX( TEMP, ZERO )
                  TEMP2 = ONE + 0.05D+0*TEMP*
     $                    ( WORK( IPN+LL ) / WORK( IPN+NQ+LL ) )**2
                  IF( TEMP2.EQ.ONE ) THEN
                     IF( IA+M-1.GT.I ) THEN
                        CALL PDNRM2( IA+M-I-1, WORK( IPN+LL ), A, I+1,
     $                               J+LL-JJ+2, DESCA, 1 )
                        WORK( IPN+NQ+LL ) = WORK( IPN+LL )
                     ELSE
                        WORK( IPN+LL ) = ZERO
                        WORK( IPN+NQ+LL ) = ZERO
                     END IF
                  ELSE
                     WORK( IPN+LL ) = WORK( IPN+LL ) * SQRT( TEMP )
                  END IF
               END IF
   90       CONTINUE
            JJ = JJ + JN - J
         END IF
         ICURCOL = MOD( ICURCOL+1, NPCOL )
*
         DO 110 K = JN+1, JA+N-1, DESCA( NB_ )
            KB = MIN( JA+N-K, DESCA( NB_ ) )
*
            IF( MYCOL.EQ.ICURCOL ) THEN
               DO 100 LL = JJ-1, JJ+KB-2
                  IF( WORK( IPN+LL ).NE.ZERO ) THEN
                     TEMP = ONE-( ABS( WORK( IPW+LL ) ) /
     $                            WORK( IPN+LL ) )**2
                     TEMP = MAX( TEMP, ZERO )
                     TEMP2 = ONE + 0.05D+0*TEMP*
     $                     ( WORK( IPN+LL ) / WORK( IPN+NQ+LL ) )**2
                     IF( TEMP2.EQ.ONE ) THEN
                        IF( IA+M-1.GT.I ) THEN
                           CALL PDNRM2( IA+M-I-1, WORK( IPN+LL ), A,
     $                                  I+1, K+LL-JJ+1, DESCA, 1 )
                           WORK( IPN+NQ+LL ) = WORK( IPN+LL )
                        ELSE
                           WORK( IPN+LL ) = ZERO
                           WORK( IPN+NQ+LL ) = ZERO
                        END IF
                     ELSE
                        WORK( IPN+LL ) = WORK( IPN+LL ) * SQRT( TEMP )
                     END IF
                  END IF
  100          CONTINUE
               JJ = JJ + KB
            END IF
            ICURCOL = MOD( ICURCOL+1, NPCOL )
*
  110    CONTINUE
*
  120 CONTINUE
  
      RANK = N - TEMP3
*
      WORK( 1 ) = DBLE( LWMIN )
*
      RETURN
*
*     End of PDGEQPF
*
      END
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     RPDGELS2:  Modified PDGELS to use custom RPDGEQPF
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      SUBROUTINE RPDGELS( TOL, TRANS, M, N, NRHS, A, IA, JA, DESCA, 
     $                    B, IB, JB, DESCB, WORK, LWORK, IPIV,
     $                    RANK, INFO )
*
*     Originally:
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     Modified date October 9, 2012
*
*     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            IA, IB, INFO, JA, JB, LWORK, M, N, NRHS, 
     $                   RANK
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( 9 ), DESCB( 9 ), IPIV( * )
      DOUBLE PRECISION   A( * ), B( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PDGELS solves overdetermined or underdetermined real linear
*  systems involving an M-by-N matrix sub( A ) = A(IA:IA+M-1,JA:JA+N-1),
*  or its transpose, using a QR or LQ factorization of sub( A ).  It is
*  assumed that sub( A ) has full rank.
*
*  The following options are provided:
*
*  1. If TRANS = 'N' and m >= n:  find the least squares solution of
*     an overdetermined system, i.e., solve the least squares problem
*                  minimize || sub( B ) - sub( A )*X ||.
*
*  2. If TRANS = 'N' and m < n:  find the minimum norm solution of
*     an underdetermined system sub( A ) * X = sub( B ).
*
*  3. If TRANS = 'T' and m >= n:  find the minimum norm solution of
*     an undetermined system sub( A )**T * X = sub( B ).
*
*  4. If TRANS = 'T' and m < n:  find the least squares solution of
*     an overdetermined system, i.e., solve the least squares problem
*                  minimize || sub( B ) - sub( A )**T * X ||.
*
*  where sub( B ) denotes B( IB:IB+M-1, JB:JB+NRHS-1 ) when TRANS = 'N'
*  and B( IB:IB+N-1, JB:JB+NRHS-1 ) otherwise. Several right hand side
*  vectors b and solution vectors x can be handled in a single call;
*  When TRANS = 'N', the solution vectors are stored as the columns of
*  the N-by-NRHS right hand side matrix sub( B ) and the M-by-NRHS
*  right hand side matrix sub( B ) otherwise.
*
*  Notes
*  =====
*
*  Each global data object is described by an associated description
*  vector.  This vector stores the information required to establish
*  the mapping between an object element and its corresponding process
*  and memory location.
*
*  Let A be a generic term for any 2D block cyclicly distributed array.
*  Such a global array has an associated description vector DESCA.
*  In the following comments, the character _ should be read as
*  "of the global array".
*
*  NOTATION        STORED IN      EXPLANATION
*  --------------- -------------- --------------------------------------
*  DTYPE_A(global) DESCA( DTYPE_ )The descriptor type.  In this case,
*                                 DTYPE_A = 1.
*  CTXT_A (global) DESCA( CTXT_ ) The BLACS context handle, indicating
*                                 the BLACS process grid A is distribu-
*                                 ted over. The context itself is glo-
*                                 bal, but the handle (the integer
*                                 value) may vary.
*  M_A    (global) DESCA( M_ )    The number of rows in the global
*                                 array A.
*  N_A    (global) DESCA( N_ )    The number of columns in the global
*                                 array A.
*  MB_A   (global) DESCA( MB_ )   The blocking factor used to distribute
*                                 the rows of the array.
*  NB_A   (global) DESCA( NB_ )   The blocking factor used to distribute
*                                 the columns of the array.
*  RSRC_A (global) DESCA( RSRC_ ) The process row over which the first
*                                 row of the array A is distributed.
*  CSRC_A (global) DESCA( CSRC_ ) The process column over which the
*                                 first column of the array A is
*                                 distributed.
*  LLD_A  (local)  DESCA( LLD_ )  The leading dimension of the local
*                                 array.  LLD_A >= MAX(1,LOCr(M_A)).
*
*  Let K be the number of rows or columns of a distributed matrix,
*  and assume that its process grid has dimension p x q.
*  LOCr( K ) denotes the number of elements of K that a process
*  would receive if K were distributed over the p processes of its
*  process column.
*  Similarly, LOCc( K ) denotes the number of elements of K that a
*  process would receive if K were distributed over the q processes of
*  its process row.
*  The values of LOCr() and LOCc() may be determined via a call to the
*  ScaLAPACK tool function, NUMROC:
*          LOCr( M ) = NUMROC( M, MB_A, MYROW, RSRC_A, NPROW ),
*          LOCc( N ) = NUMROC( N, NB_A, MYCOL, CSRC_A, NPCOL ).
*  An upper bound for these quantities may be computed by:
*          LOCr( M ) <= ceil( ceil(M/MB_A)/NPROW )*MB_A
*          LOCc( N ) <= ceil( ceil(N/NB_A)/NPCOL )*NB_A
*
*  Arguments
*  =========
*
*  TRANS   (global input) CHARACTER
*          = 'N': the linear system involves sub( A );
*          = 'T': the linear system involves sub( A )**T.
*
*  M       (global input) INTEGER
*          The number of rows to be operated on, i.e. the number of
*          rows of the distributed submatrix sub( A ). M >= 0.
*
*  N       (global input) INTEGER
*          The number of columns to be operated on, i.e. the number of
*          columns of the distributed submatrix sub( A ). N >= 0.
*
*  NRHS    (global input) INTEGER
*          The number of right hand sides, i.e. the number of columns
*          of the distributed submatrices sub( B ) and X.  NRHS >= 0.
*
*  A       (local input/local output) DOUBLE PRECISION pointer into the
*          local memory to an array of local dimension
*          ( LLD_A, LOCc(JA+N-1) ).  On entry, the M-by-N matrix A.
*          if M >= N, sub( A ) is overwritten by details of its QR
*            factorization as returned by PDGEQRF;
*          if M <  N, sub( A ) is overwritten by details of its LQ
*            factorization as returned by PDGELQF.
*
*  IA      (global input) INTEGER
*          The row index in the global array A indicating the first
*          row of sub( A ).
*
*  JA      (global input) INTEGER
*          The column index in the global array A indicating the
*          first column of sub( A ).
*
*  DESCA   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix A.
*
*  B       (local input/local output) DOUBLE PRECISION pointer into the
*          local memory to an array of local dimension
*          (LLD_B, LOCc(JB+NRHS-1)).  On entry, this array contains the
*          local pieces of the distributed matrix B of right hand side
*          vectors, stored columnwise;
*          sub( B ) is M-by-NRHS if TRANS='N', and N-by-NRHS otherwise.
*          On exit, sub( B ) is overwritten by the solution vectors,
*          stored columnwise:  if TRANS = 'N' and M >= N, rows 1 to N
*          of sub( B ) contain the least squares solution vectors; the
*          residual sum of squares for the solution in each column is
*          given by the sum of squares of elements N+1 to M in that
*          column; if TRANS = 'N' and M < N, rows 1 to N of sub( B )
*          contain the minimum norm solution vectors; if TRANS = 'T'
*          and M >= N, rows 1 to M of sub( B ) contain the minimum norm
*          solution vectors; if TRANS = 'T' and M < N, rows 1 to M of
*          sub( B ) contain the least squares solution vectors; the
*          residual sum of squares for the solution in each column is
*          given by the sum of squares of elements M+1 to N in that
*          column.
*
*  IB      (global input) INTEGER
*          The row index in the global array B indicating the first
*          row of sub( B ).
*
*  JB      (global input) INTEGER
*          The column index in the global array B indicating the
*          first column of sub( B ).
*
*  DESCB   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix B.
*
*  WORK    (local workspace/local output) DOUBLE PRECISION array,
*                                                  dimension (LWORK)
*          On exit, WORK(1) returns the minimal and optimal LWORK.
*
*  LWORK   (local or global input) INTEGER
*          The dimension of the array WORK.
*          LWORK is local input and must be at least
*          LWORK >= LTAU + MAX( LWF, LWS ) where
*          If M >= N, then
*            LTAU = NUMROC( JA+MIN(M,N)-1, NB_A, MYCOL, CSRC_A, NPCOL ),
*            LWF  = NB_A * ( MpA0 + NqA0 + NB_A )
*            LWS  = MAX( (NB_A*(NB_A-1))/2, (NRHSqB0 + MpB0)*NB_A ) +
*                   NB_A * NB_A
*          Else
*            LTAU = NUMROC( IA+MIN(M,N)-1, MB_A, MYROW, RSRC_A, NPROW ),
*            LWF  = MB_A * ( MpA0 + NqA0 + MB_A )
*            LWS  = MAX( (MB_A*(MB_A-1))/2, ( NpB0 + MAX( NqA0 +
*                   NUMROC( NUMROC( N+IROFFB, MB_A, 0, 0, NPROW ),
*                   MB_A, 0, 0, LCMP ), NRHSqB0 ) )*MB_A ) +
*                   MB_A * MB_A
*          End if
*
*          where LCMP = LCM / NPROW with LCM = ILCM( NPROW, NPCOL ),
*
*          IROFFA = MOD( IA-1, MB_A ), ICOFFA = MOD( JA-1, NB_A ),
*          IAROW = INDXG2P( IA, MB_A, MYROW, RSRC_A, NPROW ),
*          IACOL = INDXG2P( JA, NB_A, MYCOL, CSRC_A, NPCOL ),
*          MpA0 = NUMROC( M+IROFFA, MB_A, MYROW, IAROW, NPROW ),
*          NqA0 = NUMROC( N+ICOFFA, NB_A, MYCOL, IACOL, NPCOL ),
*
*          IROFFB = MOD( IB-1, MB_B ), ICOFFB = MOD( JB-1, NB_B ),
*          IBROW = INDXG2P( IB, MB_B, MYROW, RSRC_B, NPROW ),
*          IBCOL = INDXG2P( JB, NB_B, MYCOL, CSRC_B, NPCOL ),
*          MpB0 = NUMROC( M+IROFFB, MB_B, MYROW, IBROW, NPROW ),
*          NpB0 = NUMROC( N+IROFFB, MB_B, MYROW, IBROW, NPROW ),
*          NRHSqB0 = NUMROC( NRHS+ICOFFB, NB_B, MYCOL, IBCOL, NPCOL ),
*
*          ILCM, INDXG2P and NUMROC are ScaLAPACK tool functions;
*          MYROW, MYCOL, NPROW and NPCOL can be determined by calling
*          the subroutine BLACS_GRIDINFO.
*
*          If LWORK = -1, then LWORK is global input and a workspace
*          query is assumed; the routine only calculates the minimum
*          and optimal size for all work arrays. Each of these
*          values is returned in the first entry of the corresponding
*          work array, and no error message is issued by PXERBLA.
*
*  INFO    (global output) INTEGER
*          = 0:  successful exit
*          < 0:  If the i-th argument is an array and the j-entry had
*                an illegal value, then INFO = -(i*100+j), if the i-th
*                argument is a scalar and had an illegal value, then
*                INFO = -i.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY, TPSD
      INTEGER            BROW, IACOL, IAROW, IASCL, IBCOL, IBROW, IBSCL,
     $                   ICOFFA, ICOFFB, ICTXT, IPW, IROFFA, IROFFB,
     $                   LCM, LCMP, LTAU, LWF, LWMIN, LWS, MPA0, MPB0,
     $                   MYCOL, MYROW, NPB0, NPCOL, NPROW, NQA0,
     $                   NRHSQB0, SCLLEN
      DOUBLE PRECISION   ANRM, BIGNUM, BNRM, SMLNUM
*     ..
*     .. Local Arrays ..
      INTEGER            IDUM1( 2 ), IDUM2( 2 )
      DOUBLE PRECISION   RWORK( 1 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILCM
      INTEGER            INDXG2P, NUMROC
      DOUBLE PRECISION   PDLAMCH, PDLANGE
      EXTERNAL           ILCM, INDXG2P, LSAME, NUMROC, PDLAMCH,
     $                   PDLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, PCHK2MAT, PDGELQF,
     $                   PDGEQRF, PDLABAD, PDLASCL, PDLASET,
     $                   PDORMLQ, PDORMQR, PDTRSM, PXERBLA,
     $                   PDTZRZF, PDGEQPF, PDORGQR
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, ICHAR, MAX, MIN, MOD
*     ..
*     .. Executable Statements ..
*
*     Get grid parameters
*
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
*     Test the input parameters
*
      INFO = 0
      IF( NPROW.EQ.-1 ) THEN
         INFO = -( 800 + CTXT_ )
      ELSE
         CALL CHK1MAT( M, 2, N, 3, IA, JA, DESCA, 8, INFO )
         IF ( M .GE. N ) THEN
            CALL CHK1MAT( M, 2, NRHS, 4, IB, JB, DESCB, 12, INFO )
         ELSE
            CALL CHK1MAT( N, 3, NRHS, 4, IB, JB, DESCB, 12, INFO )
         ENDIF
         IF( INFO.EQ.0 ) THEN
            IROFFA = MOD( IA-1, DESCA( MB_ ) )
            ICOFFA = MOD( JA-1, DESCA( NB_ ) )
            IROFFB = MOD( IB-1, DESCB( MB_ ) )
            ICOFFB = MOD( JB-1, DESCB( NB_ ) )
            IAROW = INDXG2P( IA, DESCA( MB_ ), MYROW, DESCA( RSRC_ ),
     $                       NPROW )
            IACOL = INDXG2P( IA, DESCA( NB_ ), MYCOL, DESCA( CSRC_ ),
     $                       NPCOL )
            MPA0 = NUMROC( M+IROFFA, DESCA( MB_ ), MYROW, IAROW, NPROW )
            NQA0 = NUMROC( N+ICOFFA, DESCA( NB_ ), MYCOL, IACOL, NPCOL )
*
            IBROW = INDXG2P( IB, DESCB( MB_ ), MYROW, DESCB( RSRC_ ),
     $                       NPROW )
            IBCOL = INDXG2P( IB, DESCB( NB_ ), MYCOL, DESCB( CSRC_ ),
     $                       NPCOL )
            NRHSQB0 = NUMROC( NRHS+ICOFFB, DESCB( NB_ ), MYCOL, IBCOL,
     $                        NPCOL )
            IF( M.GE.N ) THEN
               MPB0 = NUMROC( M+IROFFB, DESCB( MB_ ), MYROW, IBROW,
     $                        NPROW )
               LTAU = NUMROC( JA+MIN(M,N)-1, DESCA( NB_ ), MYCOL,
     $                        DESCA( CSRC_ ), NPCOL )
               LWF  = DESCA( NB_ ) * ( MPA0 + NQA0 + DESCA( NB_ ) )
               LWS = MAX( ( DESCA( NB_ )*( DESCA( NB_ ) - 1 ) ) / 2,
     $               ( MPB0 + NRHSQB0 ) * DESCA( NB_ ) ) +
     $               DESCA( NB_ )*DESCA( NB_ )
            ELSE
               LCM = ILCM( NPROW, NPCOL )
               LCMP = LCM / NPROW
               NPB0 = NUMROC( N+IROFFB, DESCB( MB_ ), MYROW, IBROW,
     $                        NPROW )
               LTAU = NUMROC( IA+MIN(M,N)-1, DESCA( MB_ ), MYROW,
     $                        DESCA( RSRC_ ), NPROW )
               LWF  = DESCA( MB_ ) * ( MPA0 + NQA0 + DESCA( MB_ ) )
               LWS  = MAX( ( DESCA( MB_ )*( DESCA( MB_ ) - 1 ) ) / 2,
     $                ( NPB0 + MAX( NQA0 + NUMROC( NUMROC( N+IROFFB,
     $                DESCA( MB_ ), 0, 0, NPROW ), DESCA( MB_ ), 0, 0,
     $                LCMP ), NRHSQB0 ) )*DESCA( MB_ ) ) +
     $                DESCA( MB_ ) * DESCA( MB_ )
            END IF
            LWMIN = LTAU + MAX( LWF, LWS )
            WORK( 1 ) = DBLE( LWMIN )
            LQUERY = ( LWORK.EQ.-1 )
*
            TPSD = .TRUE.
            IF( LSAME( TRANS, 'N' ) )
     $         TPSD = .FALSE.
*
            IF( .NOT.( LSAME( TRANS, 'N' ) .OR.
     $          LSAME( TRANS, 'T' ) ) ) THEN
               INFO = -1
            ELSE IF( M.LT.0 ) THEN
               INFO = -2
            ELSE IF( N.LT.0 ) THEN
               INFO = -3
            ELSE IF( NRHS.LT.0 ) THEN
               INFO = -4
            ELSE IF( M.GE.N .AND. IROFFA.NE.IROFFB ) THEN
               INFO = -10
            ELSE IF( M.GE.N .AND. IAROW.NE.IBROW ) THEN
               INFO = -10
            ELSE IF( M.LT.N .AND. ICOFFA.NE.IROFFB ) THEN
               INFO = -10
            ELSE IF( M.GE.N .AND. DESCA( MB_ ).NE.DESCB( MB_ ) ) THEN
               INFO = -( 1200 + MB_ )
            ELSE IF( M.LT.N .AND. DESCA( NB_ ).NE.DESCB( MB_ ) ) THEN
               INFO = -( 1200 + MB_ )
            ELSE IF( ICTXT.NE.DESCB( CTXT_ ) ) THEN
               INFO = -( 1200 + CTXT_ )
            ELSE IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
               INFO = -14
            END IF
         END IF
*
         IF( .NOT.TPSD ) THEN
            IDUM1( 1 ) = ICHAR( 'N' )
         ELSE
            IDUM1( 1 ) = ICHAR( 'T' )
         END IF
         IDUM2( 1 ) = 1
         IF( LWORK.EQ.-1 ) THEN
            IDUM1( 2 ) = -1
         ELSE
            IDUM1( 2 ) = 1
         END IF
         IDUM2( 2 ) = 14
         CALL PCHK2MAT( M, 2, N, 3, IA, JA, DESCA, 8, N, 3, NRHS, 4,
     $                  IB, JB, DESCB, 12, 2, IDUM1, IDUM2, INFO )
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PDGELS', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( MIN( M, N, NRHS ).EQ.0 ) THEN
         CALL PDLASET( 'Full', MAX( M, N ), NRHS, ZERO, ZERO, B,
     $                 IB, JB, DESCB )
         RETURN
      END IF
*
*     Get machine parameters
*
      SMLNUM = PDLAMCH( ICTXT, 'S' )
      SMLNUM = SMLNUM / PDLAMCH( ICTXT, 'P' )
      BIGNUM = ONE / SMLNUM
      CALL PDLABAD( ICTXT, SMLNUM, BIGNUM )
*
*     Scale A, B if max entry outside range [SMLNUM,BIGNUM]
*
      ANRM = PDLANGE( 'M', M, N, A, IA, JA, DESCA, RWORK )
      IASCL = 0
      IF( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) THEN
*
*        Scale matrix norm up to SMLNUM
*
         CALL PDLASCL( 'G', ANRM, SMLNUM, M, N, A, IA, JA, DESCA,
     $                 INFO )
         IASCL = 1
      ELSE IF( ANRM.GT.BIGNUM ) THEN
*
*        Scale matrix norm down to BIGNUM
*
         CALL PDLASCL( 'G', ANRM, BIGNUM, M, N, A, IA, JA, DESCA,
     $                 INFO )
         IASCL = 2
      ELSE IF( ANRM.EQ.ZERO ) THEN
*
*        Matrix all zero. Return zero solution.
*
         CALL PDLASET( 'F', MAX( M, N ), NRHS, ZERO, ZERO, B, IB, JB,
     $                 DESCB )
         GO TO 10
      END IF
*
      BROW = M
      IF( TPSD )
     $   BROW = N
*
      BNRM = PDLANGE( 'M', BROW, NRHS, B, IB, JB, DESCB, RWORK )
*
      IBSCL = 0
      IF( BNRM.GT.ZERO .AND. BNRM.LT.SMLNUM ) THEN
*
*        Scale matrix norm up to SMLNUM
*
         CALL PDLASCL( 'G', BNRM, SMLNUM, BROW, NRHS, B, IB, JB,
     $                 DESCB, INFO )
         IBSCL = 1
      ELSE IF( BNRM.GT.BIGNUM ) THEN
*
*        Scale matrix norm down to BIGNUM
*
         CALL PDLASCL( 'G', BNRM, BIGNUM, BROW, NRHS, B, IB, JB,
     $                 DESCB, INFO )
         IBSCL = 2
      END IF
*
      IPW = LTAU + 1
*
      IF( M.GE.N ) THEN
*
*        compute QR factorization of A
*
!
!
!
!
!         CALL PDGEQRF( M, N, A, IA, JA, DESCA, WORK, WORK( IPW ),
!     $                 LWORK-LTAU, INFO )
         CALL RPDGEQPF( TOL, M, N, A, IA, JA, DESCA, IPIV, 
     $                  WORK, WORK( IPW ), LWORK-LTAU, RANK, INFO)
!         CALL PDTZRZF( M, N, A, IA, JA, DESCA, TAU, WORK, LWORK, 
!     $                 INFO )
!         CALL PDORGQR( M, N, RANK, A, IA, JA, DESCA, WORK, WORK( IPW ), 
!     $                 LWORK-LTAU, INFO)
!
!
!
!
!
*
*        workspace at least N, optimally N*NB
*
         IF( .NOT.TPSD ) THEN
*
*           Least-Squares Problem min || A * X - B ||
*
*           B(IB:IB+M-1,JB:JB+NRHS-1) := Q' * B(IB:IB+M-1,JB:JB+NRHS-1)
*
            CALL PDORMQR( 'Left', 'Transpose', M, NRHS, N, A, IA, JA,
     $                    DESCA, WORK, B, IB, JB, DESCB, WORK( IPW ),
     $                    LWORK-LTAU, INFO )
*
*           workspace at least NRHS, optimally NRHS*NB
*
*           B(IB:IB+N-1,JB:JB+NRHS-1) := inv(R) *
*                                        B(IB:IB+N-1,JB:JB+NRHS-1)
*
            CALL PDTRSM( 'Left', 'Upper', 'No transpose', 'Non-unit', N,
     $                   NRHS, ONE, A, IA, JA, DESCA, B, IB, JB, DESCB )
*
            SCLLEN = N
*
         ELSE
*
*           Overdetermined system of equations sub( A )' * X = sub( B )
*
*           sub( B ) := inv(R') * sub( B )
*
            CALL PDTRSM( 'Left', 'Upper', 'Transpose', 'Non-unit', N,
     $                   NRHS, ONE, A, IA, JA, DESCA, B, IB, JB, DESCB )
*
*           B(IB+N:IB+M-1,JB:JB+NRHS-1) = ZERO
*
            CALL PDLASET( 'All', M-N, NRHS, ZERO, ZERO, B, IB+N, JB,
     $                    DESCB )
*
*           B(IB:IB+M-1,JB:JB+NRHS-1) := Q(1:N,:) *
*                                        B(IB:IB+N-1,JB:JB+NRHS-1)
*
            CALL PDORMQR( 'Left', 'No transpose', M, NRHS, N, A, IA, JA,
     $                    DESCA, WORK, B, IB, JB, DESCB, WORK( IPW ),
     $                    LWORK-LTAU, INFO )
*
*           workspace at least NRHS, optimally NRHS*NB
*
            SCLLEN = M
*
         END IF
*
      ELSE
*
*        Compute LQ factorization of sub( A )
*
         CALL PDGELQF( M, N, A, IA, JA, DESCA, WORK, WORK( IPW ),
     $                 LWORK-LTAU, INFO )
*
*        workspace at least M, optimally M*NB.
*
         IF( .NOT.TPSD ) THEN
*
*           underdetermined system of equations sub( A ) * X = sub( B )
*
*           B(IB:IB+M-1,JB:JB+NRHS-1) := inv(L) *
*                                        B(IB:IB+M-1,JB:JB+NRHS-1)
*
            CALL PDTRSM( 'Left', 'Lower', 'No transpose', 'Non-unit', M,
     $                   NRHS, ONE, A, IA, JA, DESCA, B, IB, JB, DESCB )
*
*           B(IB+M:IB+N-1,JB:JB+NRHS-1) = 0
*
            CALL PDLASET( 'All', N-M, NRHS, ZERO, ZERO, B, IB+M, JB,
     $                    DESCB )
*
*           B(IB:IB+N-1,JB:JB+NRHS-1) := Q(1:N,:)' *
*                                        B(IB:IB+M-1,JB:JB+NRHS-1)
*
            CALL PDORMLQ( 'Left', 'Transpose', N, NRHS, M, A, IA, JA,
     $                    DESCA, WORK, B, IB, JB, DESCB, WORK( IPW ),
     $                    LWORK-LTAU, INFO )
*
*           workspace at least NRHS, optimally NRHS*NB
*
            SCLLEN = N
*
         ELSE
*
*           overdetermined system min || A' * X - B ||
*
*           B(IB:IB+N-1,JB:JB+NRHS-1) := Q * B(IB:IB+N-1,JB:JB+NRHS-1)
*
            CALL PDORMLQ( 'Left', 'No transpose', N, NRHS, M, A, IA, JA,
     $                    DESCA, WORK, B, IB, JB, DESCB, WORK( IPW ),
     $                    LWORK-LTAU, INFO )
*
*           workspace at least NRHS, optimally NRHS*NB
*
*           B(IB:IB+M-1,JB:JB+NRHS-1) := inv(L') *
*                                        B(IB:IB+M-1,JB:JB+NRHS-1)
*
            CALL PDTRSM( 'Left', 'Lower', 'Transpose', 'Non-unit', M,
     $                   NRHS, ONE, A, IA, JA, DESCA, B, IB, JB, DESCB )
*
            SCLLEN = M
*
         END IF
*
      END IF
*
*     Undo scaling
*
      IF( IASCL.EQ.1 ) THEN
         CALL PDLASCL( 'G', ANRM, SMLNUM, SCLLEN, NRHS, B, IB, JB,
     $                 DESCB, INFO )
      ELSE IF( IASCL.EQ.2 ) THEN
         CALL PDLASCL( 'G', ANRM, BIGNUM, SCLLEN, NRHS, B, IB, JB,
     $                 DESCB, INFO )
      END IF
      IF( IBSCL.EQ.1 ) THEN
         CALL PDLASCL( 'G', SMLNUM, BNRM, SCLLEN, NRHS, B, IB, JB,
     $                 DESCB, INFO )
      ELSE IF( IBSCL.EQ.2 ) THEN
         CALL PDLASCL( 'G', BIGNUM, BNRM, SCLLEN, NRHS, B, IB, JB,
     $                 DESCB, INFO )
      END IF
*
   10 CONTINUE
*
      WORK( 1 ) = DBLE( LWMIN )
*
      RETURN
*
*     End of PDGELS
*
      END
