! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2012-2013, Schmidt


!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     RPDORMQR:  Modified PDORMQR to be compatible with 
!     RPDGEQPF.
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      SUBROUTINE RPDORMQR( SIDE, TRANS, M, N, K, A, IA, JA, DESCA, TAU,
     $                    C, IC, JC, DESCC, WORK, LWORK, INFO )
*
*     THIS IS A MODIFIED ROUTINE
*     See http://www.netlib.org/scalapack/ for the original version
*
*     Originally:
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
      END SUBROUTINE
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
*     THIS IS A MODIFIED ROUTINE
*     See http://www.netlib.org/scalapack/ for the original version
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
  
!      RANK = N - TEMP3
        ! Calculate numerical rank
      IF ( M.LT.N ) THEN ! .LT.?
        ISZ = M
      ELSE
        ISZ = N-1
      END IF
      CALL DESCSET( DESCN, 1, DESCA( N_ ), 1, DESCA( NB_ ), 
     $              DESCA( RSRC_ ), DESCA( CSRC_ ), ICTXT, 1 )
      RANK = 0
      TEMP3 = 0
      DO 12 IJ = 1, ISZ, 1
      CALL PIELGET( 'A', 'T', TEMP3, IPIV, 1, IJ, DESCN )
        IF ( IJ.LE.TEMP3 ) THEN
          RANK = RANK + 1
        END IF
   12 CONTINUE
      
      IF ( M.GE.N ) THEN
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
      END IF
      
      ! End of numerical rank calculation
*
      WORK( 1 ) = DBLE( LWMIN )
*
      RETURN
*
*     End of PDGEQPF
*
      END SUBROUTINE
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     RPDGELS:  Heavily modified PDGELS to use custom RPDGEQPF 
!     which uses R's 'limited pivoting strategy', as well as to 
!     return the things R wants in the return
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      SUBROUTINE RPDGELS( TOL, TRANS, M, N, NRHS, 
     $                    A, IA, JA, DESCA, 
     $                    B, IB, JB, DESCB, 
     $                    EFF, FT, RSD,
     $                    TAU, WORK, LWORK, IPIV,
     $                    RANK, INFO )
*
*     THIS IS A MODIFIED ROUTINE
*     See http://www.netlib.org/scalapack/ for the original version
*
*     Originally:
*  -- ScaLAPACK routine (version 1.7) --`
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
      DOUBLE PRECISION   TOL
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( 9 ), DESCB( 9 ), IPIV( * )
      DOUBLE PRECISION   A( * ), B( * ), WORK( * ), EFF(*),
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
!         ELSE
!            CALL CHK1MAT( N, 3, NRHS, 4, IB, JB, DESCB, 12, INFO )
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
!*
!*
!*        compute QR factorization of A
!*
!
        ! Copy B over to RSD for later residual calculation
         CALL PDLACPY('All', M, NRHS, B, IB, JB, DESCB,
     $                    RSD, IB, JB, DESCB)
!
         CALL RPDGEQPF( TOL, M, N, A, IA, JA, DESCA, IPIV, 
     $                  TAU, WORK( IPW ), LWORK-LTAU, RANK, INFO)
         ! Adjust number of columns to fit numerical rank
         ITMP = N ! original N
         N = RANK
         DESCA(4) = N
*
!*           Least-Squares Problem min || A * X - B ||
!*
!*           B(IB:IB+M-1,JB:JB+NRHS-1) := Q' * B(IB:IB+M-1,JB:JB+NRHS-1)
!*
            CALL PDORMQR( 'Left', 'Transpose', M, NRHS, N, A, IA, JA,
     $                    DESCA, TAU, B, IB, JB, DESCB, WORK( IPW ),
     $                    LWORK-LTAU, INFO )
!           "effects" variable
            CALL PDLACPY('All', M, NRHS, B, IB, JB, DESCB,
     $                    EFF, IB, JB, DESCB)
!*
!*           workspace at least NRHS, optimally NRHS*NB
!*
!*           B(IB:IB+N-1,JB:JB+NRHS-1) := inv(R) *
!*                                        B(IB:IB+N-1,JB:JB+NRHS-1)
!*
            CALL PDTRSM( 'Left', 'Upper', 'No transpose', 'Non-unit', N,
     $                   NRHS, ONE, A, IA, JA, DESCA, B, IB, JB, DESCB )
!!!!!!!!!!!!!!!!!!!!!!!!!!
! Produce fitted.values = Ax = Q*(R*x)
            ! Copy over the first RANK elements of numerical soln X
            CALL PDLACPY('All', N, NRHS, B, IB, JB, DESCB,
     $                    FT, IB, JB, DESCB)
            ! Pretend A="QR" is the upper triangular R and compute R*x
            CALL PDTRMM('L', 'U', 'N', 'N',
     $                   N, NRHS, 1.0D+0, 
     $                   A, IA, JA, DESCA, 
     $                   FT, IB, JB, DESCB )
            ! Compute fitted FT = Q*(R*x)
            CALL PDORMQR( 'L', 'N', M, NRHS, N, A, IA, JA,
     $                     DESCA, TAU, FT, IB, JB, DESCB, WORK( IPW ),
     $                     LWORK-LTAU, INFO )
            ! Compute residual RSD = FT-b
            CALL PDGEADD( 'N', M, NRHS, -1.0D+0,
     $                     FT, IB, JB, DESCB, 1.0D+0, 
     $                     RSD, IB, JB, DESCB)
!!!!!!!!!!!!!!!!!!!!!!!!!!
!*
            SCLLEN = N
!*
!*     Undo scaling
!*
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
!*
   10 CONTINUE

      WORK( 1 ) = DBLE( LWMIN )
!*
      RETURN
      END SUBROUTINE
