!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     PDLAPRNT:  Printing a distributed matrix
!     
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      SUBROUTINE RPDLAPRNT( M, N, A, IA, JA, DESCA, IRPRNT, 
     $                     ICPRNT, CMATNM, NOUT,
     $                     ICTXT, MYROW, MYCOL )
     
      INTEGER       M, N, IA, JA, DESCA( 9 ), IRPRNT, ICPRNT, 
     $              NOUT, ICTXT, MYROW, MYCOL

      DOUBLE PRECISION   A(DESCA(9), *),
     $                   WORK( DESCA(9) )

*WCC  CHARACTER     CMATNM
      CHARACTER*255 CMATNM

!     External Subroutines ..
*WCC  EXTERNAL      PDLAPRNT, BLACS_BARRIER
      EXTERNAL      BPRNT, BLACS_BARRIER

!     If I'm not in the process grid, go to the end of the program
      IF( MYROW.EQ.-1 )
     $   GO TO 10

!      IF ( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
!       WRITE (*,*) 'Row-major indexing'
!       WRITE (*,*) 'So X(x,y) means row x, column y of matrix X'
!      END IF

*WCC  CALL PDLAPRNT( M, N, A, IA, JA, DESCA, IRPRNT, ICPRNT,
      CALL BPRNT( M, N, A, IA, JA, DESCA, IRPRNT, ICPRNT,
     $               CMATNM, NOUT, WORK )

   10 CONTINUE

      CALL BLACS_BARRIER( ICTXT, 'All' )

      RETURN
      END
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     PDGEMR2D:  Data reshuffling
!     
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      SUBROUTINE RPDGEMR2D( M, N, X, DESCX, B, DESCB, CTXT,
     $                      LOCRX, LOCCX, LOCRB, LOCCB )
     
      INTEGER       M, N, DESCX( 9 ), DESCB( 9 ), CTXT,
     $              LOCRX, LOCCX, LOCRB, LOCCB,
     $              IX, JX, IB, JB
      PARAMETER   ( IX = 1, JX = 1, IB = 1, JB = 1 )

      DOUBLE PRECISION   X( LOCRX, LOCCX ),
     $                   B( LOCRB, LOCCB )

!     External Subroutines ..
      EXTERNAL      PDGEMR2D

      CALL PDGEMR2D( M, N,  
     $               X, IX, JX, DESCX,
     $               B, IB, JB, DESCB,
     $               CTXT )

      RETURN
      END
