!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     RPDORMQR:  Modified PDORMQR to be compatible with 
!     RPDGEQPF.
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      SUBROUTINE RPDORMQR( SIDE, TRANS, M, N, K, A, IA, JA, DESCA, TAU,
     $                    C, IC, JC, DESCC, WORK, LWORK, INFO )
*
*     Originally:
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 25, 2001
*
*     Modified date October 16, 2012
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            IA, IC, INFO, JA, JC, K, LWORK, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCC( * )
      DOUBLE PRECISION   A( * ), C( * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PDORMQR overwrites the general real M-by-N distributed matrix
*  sub( C ) = C(IC:IC+M-1,JC:JC+N-1) with
*
*                      SIDE = 'L'            SIDE = 'R'
*  TRANS = 'N':      Q * sub( C )          sub( C ) * Q
*  TRANS = 'T':      Q**T * sub( C )       sub( C ) * Q**T
*
*  where Q is a real orthogonal distributed matrix defined as the
*  product of k elementary reflectors
*
*        Q = H(1) H(2) . . . H(k)
*
*  as returned by PDGEQRF. Q is of order M if SIDE = 'L' and of order N
*  if SIDE = 'R'.
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
*  SIDE    (global input) CHARACTER
*          = 'L': apply Q or Q**T from the Left;
*          = 'R': apply Q or Q**T from the Right.
*
*  TRANS   (global input) CHARACTER
*          = 'N':  No transpose, apply Q;
*          = 'T':  Transpose, apply Q**T.
*
*  M       (global input) INTEGER
*          The number of rows to be operated on i.e the number of rows
*          of the distributed submatrix sub( C ). M >= 0.
*
*  N       (global input) INTEGER
*          The number of columns to be operated on i.e the number of
*          columns of the distributed submatrix sub( C ). N >= 0.
*
*  K       (global input) INTEGER
*          The number of elementary reflectors whose product defines the
*          matrix Q.  If SIDE = 'L', M >= K >= 0, if SIDE = 'R',
*          N >= K >= 0.
*
*  A       (local input) DOUBLE PRECISION pointer into the local memory
*          to an array of dimension (LLD_A,LOCc(JA+K-1)). On entry, the
*          j-th column must contain the vector which defines the elemen-
*          tary reflector H(j), JA <= j <= JA+K-1, as returned by
*          PDGEQRF in the K columns of its distributed matrix
*          argument A(IA:*,JA:JA+K-1). A(IA:*,JA:JA+K-1) is modified by
*          the routine but restored on exit.
*          If SIDE = 'L', LLD_A >= MAX( 1, LOCr(IA+M-1) );
*          if SIDE = 'R', LLD_A >= MAX( 1, LOCr(IA+N-1) ).
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
*  TAU     (local input) DOUBLE PRECISION array, dimension LOCc(JA+K-1).
*          This array contains the scalar factors TAU(j) of the
*          elementary reflectors H(j) as returned by PDGEQRF.
*          TAU is tied to the distributed matrix A.
*
*  C       (local input/local output) DOUBLE PRECISION pointer into the
*          local memory to an array of dimension (LLD_C,LOCc(JC+N-1)).
*          On entry, the local pieces of the distributed matrix sub(C).
*          On exit, sub( C ) is overwritten by Q*sub( C ) or Q'*sub( C )
*          or sub( C )*Q' or sub( C )*Q.
*
*  IC      (global input) INTEGER
*          The row index in the global array C indicating the first
*          row of sub( C ).
*
*  JC      (global input) INTEGER
*          The column index in the global array C indicating the
*          first column of sub( C ).
*
*  DESCC   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix C.
*
*  WORK    (local workspace/local output) DOUBLE PRECISION array,
*                                                     dimension (LWORK)
*          On exit, WORK(1) returns the minimal and optimal LWORK.
*
*  LWORK   (local or global input) INTEGER
*          The dimension of the array WORK.
*          LWORK is local input and must be at least
*          If SIDE = 'L',
*            LWORK >= MAX( (NB_A*(NB_A-1))/2, (NqC0 + MpC0)*NB_A ) +
*                     NB_A * NB_A
*          else if SIDE = 'R',
*            LWORK >= MAX( (NB_A*(NB_A-1))/2, ( NqC0 + MAX( NpA0 +
*                     NUMROC( NUMROC( N+ICOFFC, NB_A, 0, 0, NPCOL ),
*                             NB_A, 0, 0, LCMQ ), MpC0 ) )*NB_A ) +
*                     NB_A * NB_A
*          end if
*
*          where LCMQ = LCM / NPCOL with LCM = ICLM( NPROW, NPCOL ),
*
*          IROFFA = MOD( IA-1, MB_A ), ICOFFA = MOD( JA-1, NB_A ),
*          IAROW = INDXG2P( IA, MB_A, MYROW, RSRC_A, NPROW ),
*          NpA0 = NUMROC( N+IROFFA, MB_A, MYROW, IAROW, NPROW ),
*
*          IROFFC = MOD( IC-1, MB_C ), ICOFFC = MOD( JC-1, NB_C ),
*          ICROW = INDXG2P( IC, MB_C, MYROW, RSRC_C, NPROW ),
*          ICCOL = INDXG2P( JC, NB_C, MYCOL, CSRC_C, NPCOL ),
*          MpC0 = NUMROC( M+IROFFC, MB_C, MYROW, ICROW, NPROW ),
*          NqC0 = NUMROC( N+ICOFFC, NB_C, MYCOL, ICCOL, NPCOL ),
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
*
*  INFO    (global output) INTEGER
*          = 0:  successful exit
*          < 0:  If the i-th argument is an array and the j-entry had
*                an illegal value, then INFO = -(i*100+j), if the i-th
*                argument is a scalar and had an illegal value, then
*                INFO = -i.
*
*  Alignment requirements
*  ======================
*
*  The distributed submatrices A(IA:*, JA:*) and C(IC:IC+M-1,JC:JC+N-1)
*  must verify some alignment properties, namely the following
*  expressions should be true:
*
*  If SIDE = 'L',
*    ( MB_A.EQ.MB_C .AND. IROFFA.EQ.IROFFC .AND. IAROW.EQ.ICROW )
*  If SIDE = 'R',
*    ( MB_A.EQ.NB_C .AND. IROFFA.EQ.ICOFFC )
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LEFT, LQUERY, NOTRAN
      CHARACTER          COLBTOP, ROWBTOP
      INTEGER            IAROW, ICC, ICCOL, ICOFFC, ICROW, ICTXT, IINFO,
     $                   IPW, IROFFA, IROFFC, J, J1, J2, J3, JB, JCC,
     $                   LCM, LCMQ, LWMIN, MI, MPC0, MYCOL, MYROW, NI,
     $                   NPA0, NPCOL, NPROW, NQ, NQC0
*     ..
*     .. Local Arrays ..
      INTEGER            IDUM1( 4 ), IDUM2( 4 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, PCHK2MAT, PDLARFB,
     $                   PDLARFT, PDORM2R, PB_TOPGET, PB_TOPSET, PXERBLA
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ICEIL, ILCM, INDXG2P, NUMROC
      EXTERNAL           ICEIL, ILCM, INDXG2P, LSAME, NUMROC
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
         INFO = -(900+CTXT_)
      ELSE
         LEFT = LSAME( SIDE, 'L' )
         NOTRAN = LSAME( TRANS, 'N' )
*
*        NQ is the order of Q
*
         IF( LEFT ) THEN
            NQ = M
            CALL CHK1MAT( M, 3, K, 5, IA, JA, DESCA, 9, INFO )
         ELSE
            NQ = N
            CALL CHK1MAT( N, 4, K, 5, IA, JA, DESCA, 9, INFO )
         END IF
         CALL CHK1MAT( M, 3, N, 4, IC, JC, DESCC, 14, INFO )
         IF( INFO.EQ.0 ) THEN
            IROFFA = MOD( IA-1, DESCA( MB_ ) )
            IROFFC = MOD( IC-1, DESCC( MB_ ) )
            ICOFFC = MOD( JC-1, DESCC( NB_ ) )
            IAROW = INDXG2P( IA, DESCA( MB_ ), MYROW, DESCA( RSRC_ ),
     $                       NPROW )
            ICROW = INDXG2P( IC, DESCC( MB_ ), MYROW, DESCC( RSRC_ ),
     $                       NPROW )
            ICCOL = INDXG2P( JC, DESCC( NB_ ), MYCOL, DESCC( CSRC_ ),
     $                       NPCOL )
            MPC0 = NUMROC( M+IROFFC, DESCC( MB_ ), MYROW, ICROW, NPROW )
            NQC0 = NUMROC( N+ICOFFC, DESCC( NB_ ), MYCOL, ICCOL, NPCOL )
*
            IF( LEFT ) THEN
               LWMIN = MAX( ( DESCA( NB_ ) * ( DESCA( NB_ ) - 1 ) ) / 2,
     $                      ( MPC0 + NQC0 ) * DESCA( NB_ ) ) +
     $                 DESCA( NB_ ) * DESCA( NB_ )
            ELSE
               NPA0 = NUMROC( N+IROFFA, DESCA( MB_ ), MYROW, IAROW,
     $                        NPROW )
               LCM = ILCM( NPROW, NPCOL )
               LCMQ = LCM / NPCOL
               LWMIN =  MAX( ( DESCA( NB_ ) * ( DESCA( NB_ ) - 1 ) )
     $                  / 2, ( NQC0 + MAX( NPA0 + NUMROC( NUMROC(
     $                  N+ICOFFC, DESCA( NB_ ), 0, 0, NPCOL ),
     $                  DESCA( NB_ ), 0, 0, LCMQ ), MPC0 ) ) *
     $                  DESCA( NB_ ) ) + DESCA( NB_ ) * DESCA( NB_ )
            END IF
*
            WORK( 1 ) = DBLE( LWMIN )
            LQUERY = ( LWORK.EQ.-1 )
            IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
               INFO = -1
            ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) ) THEN
               INFO = -2
            ELSE IF( K.LT.0 .OR. K.GT.NQ ) THEN
               INFO = -5
            ELSE IF( .NOT.LEFT .AND. DESCA( MB_ ).NE.DESCC( NB_ ) ) THEN
               INFO = -(900+NB_)
            ELSE IF( LEFT .AND. IROFFA.NE.IROFFC ) THEN
               INFO = -12
            ELSE IF( LEFT .AND. IAROW.NE.ICROW ) THEN
               INFO = -12
            ELSE IF( .NOT.LEFT .AND. IROFFA.NE.ICOFFC ) THEN
               INFO = -13
            ELSE IF( LEFT .AND. DESCA( MB_ ).NE.DESCC( MB_ ) ) THEN
               INFO = -(1400+MB_)
            ELSE IF( ICTXT.NE.DESCC( CTXT_ ) ) THEN
               INFO = -(1400+CTXT_)
            ELSE IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
               INFO = -16
            END IF
         END IF
*
         IF( LEFT ) THEN
            IDUM1( 1 ) = ICHAR( 'L' )
         ELSE
            IDUM1( 1 ) = ICHAR( 'R' )
         END IF
         IDUM2( 1 ) = 1
         IF( NOTRAN ) THEN
            IDUM1( 2 ) = ICHAR( 'N' )
         ELSE
            IDUM1( 2 ) = ICHAR( 'T' )
         END IF
         IDUM2( 2 ) = 2
         IDUM1( 3 ) = K
         IDUM2( 3 ) = 5
         IF( LWORK.EQ.-1 ) THEN
            IDUM1( 4 ) = -1
         ELSE
            IDUM1( 4 ) = 1
         END IF
         IDUM2( 4 ) = 16
         IF( LEFT ) THEN
            CALL PCHK2MAT( M, 3, K, 5, IA, JA, DESCA, 9, M, 3, N, 4, IC,
     $                     JC, DESCC, 14, 4, IDUM1, IDUM2, INFO )
         ELSE
            CALL PCHK2MAT( N, 4, K, 5, IA, JA, DESCA, 9, M, 3, N, 4, IC,
     $                     JC, DESCC, 14, 4, IDUM1, IDUM2, INFO )
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PDORMQR', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 )
     $   RETURN
*
      CALL PB_TOPGET( ICTXT, 'Broadcast', 'Rowwise', ROWBTOP )
      CALL PB_TOPGET( ICTXT, 'Broadcast', 'Columnwise', COLBTOP )
*
      IF( ( LEFT .AND. .NOT.NOTRAN ) .OR.
     $    ( .NOT.LEFT .AND. NOTRAN ) ) THEN
         J1 = MIN( ICEIL( JA, DESCA( NB_ ) ) * DESCA( NB_ ), JA+K-1 )
     $                    + 1
         J2 = JA+K-1
         J3 = DESCA( NB_ )
      ELSE
         J1 = MAX( ( (JA+K-2) / DESCA( NB_ ) ) * DESCA( NB_ ) + 1, JA )
         J2 = MIN( ICEIL( JA, DESCA( NB_ ) ) * DESCA( NB_ ), JA+K-1 )
     $                    + 1
         J3 = -DESCA( NB_ )
      END IF
*
      IF( LEFT ) THEN
         NI  = N
         JCC = JC
         IF( NOTRAN ) THEN
            CALL PB_TOPSET( ICTXT, 'Broadcast', 'Rowwise', 'D-ring' )
         ELSE
            CALL PB_TOPSET( ICTXT, 'Broadcast', 'Rowwise', 'I-ring' )
         END IF
         CALL PB_TOPSET( ICTXT, 'Broadcast', 'Columnwise', ' ' )
      ELSE
         MI  = M
         ICC = IC
      END IF
*
*     Use unblocked code for the first block if necessary
*
      IF( ( LEFT .AND. .NOT.NOTRAN ) .OR. ( .NOT.LEFT .AND. NOTRAN ) )
     $   CALL PDORM2R( SIDE, TRANS, M, N, J1-JA, A, IA, JA, DESCA, TAU,
     $                 C, IC, JC, DESCC, WORK, LWORK, IINFO )
*
      IPW = DESCA( NB_ ) * DESCA( NB_ ) + 1
      DO 10 J = J1, J2, J3
         JB = MIN( DESCA( NB_ ), K-J+JA )
*
*        Form the triangular factor of the block reflector
*        H = H(j) H(j+1) . . . H(j+jb-1)
*
         CALL PDLARFT( 'Forward', 'Columnwise', NQ-J+JA, JB, A,
     $                 IA+J-JA, J, DESCA, TAU, WORK, WORK( IPW ) )
         IF( LEFT ) THEN
*
*           H or H' is applied to C(ic+j-ja:ic+m-1,jc:jc+n-1)
*
            MI  = M - J + JA
            ICC = IC + J - JA
         ELSE
*
*           H or H' is applied to C(ic:ic+m-1,jc+j-ja:jc+n-1)
*
            NI  = N - J + JA
            JCC = JC + J - JA
         END IF
*
*        Apply H or H'
*
         CALL PDLARFB( SIDE, TRANS, 'Forward', 'Columnwise', MI, NI,
     $                JB, A, IA+J-JA, J, DESCA, WORK, C, ICC, JCC,
     $                DESCC, WORK( IPW ) )
   10 CONTINUE
*
*     Use unblocked code for the last block if necessary
*
      IF( ( LEFT .AND. NOTRAN ) .OR. ( .NOT.LEFT .AND. .NOT.NOTRAN ) )
     $   CALL PDORM2R( SIDE, TRANS, M, N, J2-JA, A, IA, JA, DESCA, TAU,
     $                 C, IC, JC, DESCC, WORK, LWORK, IINFO )
*
      CALL PB_TOPSET( ICTXT, 'Broadcast', 'Rowwise', ROWBTOP )
      CALL PB_TOPSET( ICTXT, 'Broadcast', 'Columnwise', COLBTOP )
*
      WORK( 1 ) = DBLE( LWMIN )
*
      RETURN
*
*     End of PDORMQR
*
      END
!
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
     $                   IJ, TEMP3
      DOUBLE PRECISION   AJJ, ALPHA, TEMP, TEMP2
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
     $                   PXERBLA, PDELGET, DGSUM2D, PIELGET
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
!      TEMP3 = 0
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
!                TEMP3 = MAX0(TEMP3, IJ - J) ! for tracking "N - RANK"
                EXIT
              END IF
   11 CONTINUE
         ELSE
!            TEMP3 = TEMP3 + N-J
            PVT = J
!            CALL PDELGET( 'A', 'T', TEMP, WORK(IPN), 1, J, DESCN )
!            IF( TEMP.LT.TOL .AND. TEMP3.EQ.0  ) THEN
!              TEMP3 = 1
!            END IF
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
  
!      RANK = N - TEMP3
        ! Calculate numerical rank
        CALL DESCSET( DESCN, 1, DESCA( N_ ), 1, DESCA( NB_ ), 
     $              DESCA( RSRC_ ), DESCA( CSRC_ ), ICTXT, 1 )
        RANK = 0
        TEMP3 = 0
        DO 12 IJ = 1, N-1, 1
        CALL PIELGET( 'A', 'T', TEMP3, IPIV, 1, IJ, DESCN )
          IF ( IJ.LE.TEMP3 ) THEN
            RANK = RANK + 1
          END IF
   12 CONTINUE
!        CALL PIELGET( 'A', 'T', TEMP3, IPIV, 1, N, DESCN )
        CALL PDELGET( 'A', 'T', TEMP, WORK(IPN), 1, N, DESCN )
          IF( TEMP.LT.TOL ) THEN
            TEMP3 = 0
          ELSE
            TEMP3 = 1
          END IF
      RANK = RANK + TEMP3
      IF ( RANK.EQ.0 ) THEN
        RANK = 1
      END IF
      ! End of numerical rank calculation
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
!     RPDGELS:  Modified PDGELS to use custom RPDGEQPF
!     as well as many other modifications...
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      SUBROUTINE RPDGELS( TOL, TRANS, M, N, NRHS, 
     $                    A, IA, JA, DESCA, 
     $                    B, IB, JB, DESCB, 
     $                    FT, RSD,
     $                    TAU, WORK, LWORK, IPIV,
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
      DOUBLE PRECISION   A( * ), B( * ), WORK( * ),
     $                   FT ( * ), RSD( * ), TAU( * )
!
!     Removed the massive explanation because I got tired of looking
!     at it.
!
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
     $                   NRHSQB0, SCLLEN, ITMP
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
     $                   PDTZRZF, PDGEQPF, PDLACPY, RPDFTTD,
     $                   PDGEADD
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
        ! Copy B over to RSD for later residual calculation
         CALL PDLACPY('All', M, NRHS, B, IB, JB, DESCB,
     $                    RSD, IB, JB, DESCB)
!
!         CALL PDGEQRF( M, N, A, IA, JA, DESCA, WORK, WORK( IPW ),
!     $                 LWORK-LTAU, INFO )
         CALL RPDGEQPF( TOL, M, N, A, IA, JA, DESCA, IPIV, 
     $                  WORK, WORK( IPW ), LWORK-LTAU, RANK, INFO)
!         CALL PDTZRZF( M, N, A, IA, JA, DESCA, TAU, WORK, LWORK, 
!     $                 INFO )
!
!
!
!
!
         ITMP = N
         ! Adjust number of columns to fit numerical rank
         N = RANK
         DESCA(4) = N

    !      CALL PDLACPY('All', M, NRHS, B, IB, JB, DESCB,
    !     $              RSD, IB, JB, DESCB)
      
!          CALL RPDFTTD( M, N, NRHS, 
!     $              A, IA, JA, DESCA,
!     $              FT, RSD, IB, JB, DESCB,
!     $              WORK, WORK( IPW ), LWORK-LTAU, 
!     $              INFO )
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
!!!!!!!! DGEMM error ?!?!
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
!!!!!!!!!!!!!!!!!!!!!!!!!!
! Produce fitted.values = Ax = Q*(R*x)
            ! Copy over the first RANK elements of numerical soln X
            CALL PDLACPY('All', N, NRHS, B, IB, JB, DESCB,
     $                    FT, IB, JB, DESCB)
            ! Pretend A="QR" is the upper triangular R and compute R*FT
            CALL PDTRMM('L', 'U', 'N', 'N',
     $                   N, NRHS, 1.0D+0, 
     $                   A, IA, JA, DESCA, 
     $                   FT, IB, JB, DESCB )
            ! Compute Q*(R*FT)
            CALL PDORMQR( 'L', 'N', M, NRHS, N, A, IA, JA,
     $                     DESCA, WORK, FT, IB, JB, DESCB, WORK( IPW ),
     $                     LWORK-LTAU, INFO )
            ! Compute residual FT-b
            CALL PDGEADD( 'N', M, NRHS, -1.0D+0,
     $                     FT, IB, JB, DESCB, 1.0D+0, 
     $                     RSD, IB, JB, DESCB)
!(trans, m, n, alpha, a, ia, ja, desca, beta, c, ic, jc, descc)
!!!!!!!!!!!!!!!!!!!!!!!!!!
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
      
      ! pass back TAU
      TAU(1:IPW-1) = WORK(1:IPW-1)
      WORK( 1 ) = DBLE( LWMIN )
*
      RETURN
*
*     End of PDGELS
*
      END



!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     RPDFTTD:  Get fitted values and residuals from RPDGELS fit
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      SUBROUTINE RPDFTTD( M, N, NRHS,
     $                    A, IA, JA, DESCA, 
     $                    FT, RSD, IFT, JFT, DESCFT,
     $                    TAU, WORK, LWORK,
     $                    INFO )
!
      INTEGER            M, N, NRHS, IA, JA, IFT, JFT, 
     $                   LWORK, INFO
      INTEGER            DESCA( 9 ), DESCFT( 9 )
      DOUBLE PRECISION   A( * ), FT( * ), RSD( * ), TAU( * ),
     $                   WORK( LWORK )
      DOUBLE PRECISION   TMP
      EXTERNAL           PDTRMV, PDORMQR, PDTRMM
      
      ! Computes fitted values and residuals in the context of 
      ! minimizing Ax-b in l2 norm.  Let A=QR.
      
      ! Compute Rx
      call PDTRMM('L', 'U', 'N', 'N',
     $             N, NRHS, 1.0D0, 
     $             A, IA, JA, DESCA, 
     $             FT, IFT, JFT, DESCFT )
      
!      CALL PDTRMV('U', 'N', 'N', N, 
!     $             A, IA, JA, DESCA, 
!     $             FT, IFT, JFT, DESCFT, 1)
      WRITE (*,*) "adsfalsdfkjaslkdf"
      
      TMP = WORK(1)
      
      ! Copmute Q(Rx) = fitted values
      CALL PDORMQR( 'L', 'N', M, NRHS, N, 
     $               A, IA, JA, DESCA, TAU
     $               FT, IFT, JFT, DESCFT, 
     $               WORK, LWORK, INFO )
      WRITE (*,*) "adsfalsdfkjaslkdf"
      
      WORK(1) = TMP
      
      RETURN
      END







































      SUBROUTINE PDORGQR( M, N, K, A, IA, JA, DESCA, TAU, WORK, LWORK,
     $                    INFO )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 25, 2001
*
*     .. Scalar Arguments ..
      INTEGER            IA, INFO, JA, K, LWORK, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * )
      DOUBLE PRECISION   A( * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PDORGQR generates an M-by-N real distributed matrix Q denoting
*  A(IA:IA+M-1,JA:JA+N-1) with orthonormal columns, which is defined as
*  the first N columns of a product of K elementary reflectors of order
*  M
*
*        Q  =  H(1) H(2) . . . H(k)
*
*  as returned by PDGEQRF.
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
*          The number of rows to be operated on i.e the number of rows
*          of the distributed submatrix Q. M >= 0.
*
*  N       (global input) INTEGER
*          The number of columns to be operated on i.e the number of
*          columns of the distributed submatrix Q. M >= N >= 0.
*
*  K       (global input) INTEGER
*          The number of elementary reflectors whose product defines the
*          matrix Q. N >= K >= 0.
*
*  A       (local input/local output) DOUBLE PRECISION pointer into the
*          local memory to an array of dimension (LLD_A,LOCc(JA+N-1)).
*          On entry, the j-th column must contain the vector which
*          defines the elementary reflector H(j), JA <= j <= JA+K-1, as
*          returned by PDGEQRF in the K columns of its distributed
*          matrix argument A(IA:*,JA:JA+K-1). On exit, this array
*          contains the local pieces of the M-by-N distributed matrix Q.
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
*  TAU     (local input) DOUBLE PRECISION array, dimension LOCc(JA+K-1)
*          This array contains the scalar factors TAU(j) of the
*          elementary reflectors H(j) as returned by PDGEQRF.
*          TAU is tied to the distributed matrix A.
*
*  WORK    (local workspace/local output) DOUBLE PRECISION array,
*                                                   dimension (LWORK)
*          On exit, WORK(1) returns the minimal and optimal LWORK.
*
*  LWORK   (local or global input) INTEGER
*          The dimension of the array WORK.
*          LWORK is local input and must be at least
*          LWORK >= NB_A * ( NqA0 + MpA0 + NB_A ), where
*
*          IROFFA = MOD( IA-1, MB_A ), ICOFFA = MOD( JA-1, NB_A ),
*          IAROW = INDXG2P( IA, MB_A, MYROW, RSRC_A, NPROW ),
*          IACOL = INDXG2P( JA, NB_A, MYCOL, CSRC_A, NPCOL ),
*          MpA0 = NUMROC( M+IROFFA, MB_A, MYROW, IAROW, NPROW ),
*          NqA0 = NUMROC( N+ICOFFA, NB_A, MYCOL, IACOL, NPCOL ),
*
*          INDXG2P and NUMROC are ScaLAPACK tool functions;
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
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      CHARACTER          COLBTOP, ROWBTOP
      INTEGER            I, IACOL, IAROW, ICTXT, IINFO, IPW, J, JB, JL,
     $                   JN, LWMIN, MPA0, MYCOL, MYROW, NPCOL, NPROW,
     $                   NQA0
*     ..
*     .. Local Arrays ..
      INTEGER            IDUM1( 2 ), IDUM2( 2 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, PCHK1MAT, PDLARFB,
     $                   PDLARFT, PDLASET, PDORG2R, PB_TOPGET,
     $                   PB_TOPSET, PXERBLA
*     ..
*     .. External Functions ..
      INTEGER            ICEIL, INDXG2P, NUMROC
      EXTERNAL           ICEIL, INDXG2P, NUMROC
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX, MIN, MOD
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
         INFO = -(700+CTXT_)
      ELSE
         CALL CHK1MAT( M, 1, N, 2, IA, JA, DESCA, 7, INFO )
         IF( INFO.EQ.0 ) THEN
            IAROW = INDXG2P( IA, DESCA( MB_ ), MYROW, DESCA( RSRC_ ),
     $                       NPROW )
            IACOL = INDXG2P( JA, DESCA( NB_ ), MYCOL, DESCA( CSRC_ ),
     $                       NPCOL )
            MPA0 = NUMROC( M+MOD( IA-1, DESCA( MB_ ) ), DESCA( MB_ ),
     $                     MYROW, IAROW, NPROW )
            NQA0 = NUMROC( N+MOD( JA-1, DESCA( NB_ ) ), DESCA( NB_ ),
     $                     MYCOL, IACOL, NPCOL )
            LWMIN = DESCA( NB_ ) * ( MPA0 + NQA0 + DESCA( NB_ ) )
*
            WORK( 1 ) = DBLE( LWMIN )
            LQUERY = ( LWORK.EQ.-1 )
            IF( N.GT.M ) THEN
               INFO = -2
            ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
               INFO = -3
            ELSE IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
               INFO = -10
            END IF
         END IF
         IDUM1( 1 ) = K
         IDUM2( 1 ) = 3
         IF( LWORK.EQ.-1 ) THEN
            IDUM1( 2 ) = -1
         ELSE
            IDUM1( 2 ) = 1
         END IF
         IDUM2( 2 ) = 10
         CALL PCHK1MAT( M, 1, N, 2, IA, JA, DESCA, 7, 2, IDUM1, IDUM2,
     $                  INFO )
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PDORGQR', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.LE.0 )
     $   RETURN
*
      IPW = DESCA( NB_ )*DESCA( NB_ ) + 1
      JN = MIN( ICEIL( JA, DESCA( NB_ ) ) * DESCA( NB_ ), JA+K-1 )
      JL = MAX( ( (JA+K-2) / DESCA( NB_ ) ) * DESCA( NB_ ) + 1, JA )
      CALL PB_TOPGET( ICTXT, 'Broadcast', 'Rowwise', ROWBTOP )
      CALL PB_TOPGET( ICTXT, 'Broadcast', 'Columnwise', COLBTOP )
      CALL PB_TOPSET( ICTXT, 'Broadcast', 'Rowwise', 'D-ring' )
      CALL PB_TOPSET( ICTXT, 'Broadcast', 'Columnwise', ' ' )
*
      CALL PDLASET( 'All', JL-JA, JA+N-JL, ZERO, ZERO, A, IA, JL,
     $              DESCA )
*
*     Use unblocked code for the last or only block.
*
      CALL PDORG2R( M-JL+JA, JA+N-JL, JA+K-JL, A, IA+JL-JA, JL, DESCA,
     $              TAU, WORK, LWORK, IINFO )
*
*     Is there at least one block of columns to loop over ?
*
      IF( JL.GT.JN+1 ) THEN
*
*        Use blocked code
*
         DO 10 J = JL-DESCA( NB_ ), JN+1, -DESCA( NB_ )
            JB = MIN( DESCA( NB_ ), JA+N-J )
            I = IA + J - JA
*
            IF( J+JB.LE.JA+N-1 ) THEN
*
*              Form the triangular factor of the block reflector
*              H = H(j) H(j+1) . . . H(j+jb-1)
*
               CALL PDLARFT( 'Forward', 'Columnwise', M-I+IA, JB, A, I,
     $                       J, DESCA, TAU, WORK, WORK( IPW ) )
*
*              Apply H to A(i:ia+m-1,j+jb:ja+n-1) from the left
*
               CALL PDLARFB( 'Left', 'No transpose', 'Forward',
     $                       'Columnwise', M-I+IA, N-J-JB+JA, JB, A, I,
     $                       J, DESCA, WORK, A, I, J+JB, DESCA,
     $                       WORK( IPW ) )
            END IF
*
*           Apply H to rows i:ia+m-1 of current block
*
            CALL PDORG2R( M-I+IA, JB, JB, A, I, J, DESCA, TAU, WORK,
     $                    LWORK, IINFO )
*
*           Set rows ia:i-1 of current block to zero
*
            CALL PDLASET( 'All', I-IA, JB, ZERO, ZERO, A, IA, J, DESCA )
*
   10    CONTINUE
*
      END IF
*
*     Handle first block separately
*
      IF( JL.GT.JA ) THEN
*
         JB = JN - JA + 1
*
*        Form the triangular factor of the block reflector
*        H = H(j) H(j+1) . . . H(j+jb-1)
*
         CALL PDLARFT( 'Forward', 'Columnwise', M, JB, A, IA, JA, DESCA,
     $                 TAU, WORK, WORK( IPW ) )
*
*        Apply H to A(ia:ia+m-1,ja+jb:ja+n-1) from the left
*
         CALL PDLARFB( 'Left', 'No transpose', 'Forward', 'Columnwise',
     $                 M, N-JB, JB, A, IA, JA, DESCA, WORK, A, IA,
     $                 JA+JB, DESCA, WORK( IPW ) )
*
*        Apply H to rows ia:ia+m-1 of current block
*
         CALL PDORG2R( M, JB, JB, A, IA, JA, DESCA, TAU, WORK, LWORK,
     $                 IINFO )
*
      END IF
*
      CALL PB_TOPGET( ICTXT, 'Broadcast', 'Rowwise', ROWBTOP )
      CALL PB_TOPGET( ICTXT, 'Broadcast', 'Columnwise', COLBTOP )
*
      WORK( 1 ) = DBLE( LWMIN )
*
      RETURN
*
*     End of PDORGQR
*
      END









      SUBROUTINE PDTZRZF( M, N, A, IA, JA, DESCA, TAU, WORK, LWORK,
     $                    INFO )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 25, 2001
*
*     .. Scalar Arguments ..
      INTEGER            IA, INFO, JA, LWORK, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * )
      DOUBLE PRECISION   A( * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PDTZRZF reduces the M-by-N ( M<=N ) real upper trapezoidal matrix
*  sub( A ) = A(IA:IA+M-1,JA:JA+N-1) to upper triangular form by means
*  of orthogonal transformations.
*
*  The upper trapezoidal matrix sub( A ) is factored as
*
*     sub( A ) = ( R  0 ) * Z,
*
*  where Z is an N-by-N orthogonal matrix and R is an M-by-M upper
*  triangular matrix.
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
*          sub( A ) which is to be factored. On exit, the leading M-by-M
*          upper triangular part of sub( A ) contains the upper trian-
*          gular matrix R, and elements M+1 to N of the first M rows of
*          sub( A ), with the array TAU, represent the orthogonal matrix
*          Z as a product of M elementary reflectors.
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
*  TAU     (local output) DOUBLE PRECISION array, dimension LOCr(IA+M-1)
*          This array contains the scalar factors of the elementary
*          reflectors. TAU is tied to the distributed matrix A.
*
*  WORK    (local workspace/local output) DOUBLE PRECISION array,
*                                                    dimension (LWORK)
*          On exit, WORK(1) returns the minimal and optimal LWORK.
*
*  LWORK   (local or global input) INTEGER
*          The dimension of the array WORK.
*          LWORK is local input and must be at least
*          LWORK >= MB_A * ( Mp0 + Nq0 + MB_A ), where
*
*          IROFF = MOD( IA-1, MB_A ), ICOFF = MOD( JA-1, NB_A ),
*          IAROW = INDXG2P( IA, MB_A, MYROW, RSRC_A, NPROW ),
*          IACOL = INDXG2P( JA, NB_A, MYCOL, CSRC_A, NPCOL ),
*          Mp0   = NUMROC( M+IROFF, MB_A, MYROW, IAROW, NPROW ),
*          Nq0   = NUMROC( N+ICOFF, NB_A, MYCOL, IACOL, NPCOL ),
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
*  The  factorization is obtained by Householder's method.  The kth
*  transformation matrix, Z( k ), which is used to introduce zeros into
*  the (m - k + 1)th row of sub( A ), is given in the form
*
*     Z( k ) = ( I     0   ),
*              ( 0  T( k ) )
*
*  where
*
*     T( k ) = I - tau*u( k )*u( k )',   u( k ) = (   1    ),
*                                                 (   0    )
*                                                 ( z( k ) )
*
*  tau is a scalar and z( k ) is an ( n - m ) element vector.
*  tau and z( k ) are chosen to annihilate the elements of the kth row
*  of sub( A ).
*
*  The scalar tau is returned in the kth element of TAU and the vector
*  u( k ) in the kth row of sub( A ), such that the elements of z( k )
*  are in  a( k, m + 1 ), ..., a( k, n ). The elements of R are returned
*  in the upper triangular part of sub( A ).
*
*  Z is given by
*
*     Z =  Z( 1 ) * Z( 2 ) * ... * Z( m ).
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      CHARACTER          COLBTOP, ROWBTOP
      INTEGER            I, IACOL, IAROW, IB, ICTXT, IIA, IL, IN, IPW,
     $                   IROFFA, J, JM1, L, LWMIN, MP0, MYCOL, MYROW,
     $                   NPCOL, NPROW, NQ0
*     ..
*     .. Local Arrays ..
      INTEGER            IDUM1( 1 ), IDUM2( 1 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, INFOG1L, PCHK1MAT,
     $                   PDLATRZ, PDLARZB, PDLARZT, PB_TOPGET,
     $                   PB_TOPSET, PXERBLA
*     ..
*     .. External Functions ..
      INTEGER            ICEIL, INDXG2P, NUMROC
      EXTERNAL           ICEIL, INDXG2P, NUMROC
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX, MIN, MOD
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
            IROFFA = MOD( IA-1, DESCA( MB_ ) )
            IAROW = INDXG2P( IA, DESCA( MB_ ), MYROW, DESCA( RSRC_ ),
     $                       NPROW )
            IACOL = INDXG2P( JA, DESCA( NB_ ), MYCOL, DESCA( CSRC_ ),
     $                       NPCOL )
            MP0 = NUMROC( M+IROFFA, DESCA( MB_ ), MYROW, IAROW, NPROW )
            NQ0 = NUMROC( N+MOD( JA-1, DESCA( NB_ ) ), DESCA( NB_ ),
     $                    MYCOL, IACOL, NPCOL )
            LWMIN = DESCA( MB_ ) * ( MP0 + NQ0 + DESCA( MB_ ) )
*
            WORK( 1 ) = DBLE( LWMIN )
            LQUERY = ( LWORK.EQ.-1 )
            IF( N.LT.M ) THEN
               INFO = -2
            ELSE IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
               INFO = -9
            END IF
         END IF
         IF( LQUERY ) THEN
            IDUM1( 1 ) = -1
         ELSE
            IDUM1( 1 ) = 1
         END IF
         IDUM2( 1 ) = 9
         CALL PCHK1MAT( M, 1, N, 2, IA, JA, DESCA, 6, 1, IDUM1, IDUM2,
     $                  INFO )
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PDTZRZF', -INFO )
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
      IF( M.EQ.N ) THEN
*
         CALL INFOG1L( IA, DESCA( MB_ ), NPROW, MYROW, DESCA( RSRC_ ),
     $                 IIA, IAROW )
         IF( MYROW.EQ.IAROW )
     $      MP0 = MP0 - IROFFA
         DO 10 I = IIA, IIA+MP0-1
            TAU( I ) = ZERO
   10    CONTINUE
*
      ELSE
*
         L = N-M
         JM1 = JA + MIN( M+1, N ) - 1
         IPW = DESCA( MB_ ) * DESCA( MB_ ) + 1
         IN = MIN( ICEIL( IA, DESCA( MB_ ) ) * DESCA( MB_ ), IA+M-1 )
         IL = MAX( ( (IA+M-2) / DESCA( MB_ ) ) * DESCA( MB_ ) + 1, IA )
         CALL PB_TOPGET( ICTXT, 'Broadcast', 'Rowwise', ROWBTOP )
         CALL PB_TOPGET( ICTXT, 'Broadcast', 'Columnwise', COLBTOP )
         CALL PB_TOPSET( ICTXT, 'Broadcast', 'Rowwise', ' ' )
         CALL PB_TOPSET( ICTXT, 'Broadcast', 'Columnwise', 'D-ring' )
*
*        Use blocked code initially
*
         DO 20 I = IL, IN+1, -DESCA( MB_ )
            IB = MIN( IA+M-I, DESCA( MB_ ) )
            J = JA + I - IA
*
*           Compute the complete orthogonal factorization of the current
*           block A(i:i+ib-1,j:ja+n-1)
*
            CALL PDLATRZ( IB, JA+N-J, L, A, I, J, DESCA, TAU, WORK )
*
            IF( I.GT.IA ) THEN
*
*              Form the triangular factor of the block reflector
*              H = H(i+ib-1) . . . H(i+1) H(i)
*
               CALL PDLARZT( 'Backward', 'Rowwise', L, IB, A, I, JM1,
     $                       DESCA, TAU, WORK, WORK( IPW ) )
*
*              Apply H to A(ia:i-1,j:ja+n-1) from the right
*
               CALL PDLARZB( 'Right', 'No transpose', 'Backward',
     $                       'Rowwise', I-IA, JA+N-J, IB, L, A, I, JM1,
     $                       DESCA, WORK, A, IA, J, DESCA, WORK( IPW ) )
            END IF
*
   20    CONTINUE
*
*        Use unblocked code to factor the last or only block
*
         CALL PDLATRZ( IN-IA+1, N, N-M, A, IA, JA, DESCA, TAU, WORK )
*
         CALL PB_TOPSET( ICTXT, 'Broadcast', 'Rowwise', ROWBTOP )
         CALL PB_TOPSET( ICTXT, 'Broadcast', 'Columnwise', COLBTOP )
*
      END IF
*
      WORK( 1 ) = DBLE( LWMIN )
*
      RETURN
*
*     End of PDTZRZF
*
      END
