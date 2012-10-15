!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     PDTRAN:  Matrix transpose
!     beta*C + alpha*t(A) --> C
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      SUBROUTINE RPDTRAN(A, C, ICTXT, MYROW, MYCOL,
     $                   DESCA, DESCC, M, N )

      INTEGER       IA, JA, IC, JC, M, N, 
     $              ICTXT, MYROW, MYCOL,
     $              DESCA( 9 ), DESCC( 9 )
      PARAMETER     ( IA = 1, JA = 1, IC = 1, JC = 1 )

      DOUBLE PRECISION   ONE, ZERO,
     $                   A( DESCA(9), * ),
     $                   C( DESCC(9), * )
      PARAMETER      ( ONE =  1.0D+0, ZERO =  0.0D+0 )

!     External Subroutines ..
      EXTERNAL      PDTRAN

!     If I'm not in the process grid, go to the end of the program
      IF( MYROW.EQ.-1 )
     $   GO TO 10

!     Compute transpose
      CALL PDTRAN( M, N, ONE, 
     $             A, IA, JA, DESCA, ZERO,
     $             C, IC, JC, DESCC )

   10 CONTINUE

      RETURN
      END
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     PDGEMM:  Matrix-Matrix multiplication
!     (alpha * op(A)) %*% op(B) + (beta * C) --> C
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      SUBROUTINE RPDGEMM( A, B, C, ICTXT, MYROW, MYCOL,
     $                   DESCA, DESCB, DESCC,
     $                   M, N, K )

      INTEGER       IA, JA, IB, JB, IC, JC, M, N, K, 
     $              ICTXT, MYROW, MYCOL,
     $              DESCA( 9 ), DESCB( 9 ), DESCC( 9 )
      PARAMETER     ( IA = 1, JA = 1, IB = 1, JB = 1 , 
     $                IC = 1, JC = 1 )

      DOUBLE PRECISION   ONE, ZERO,
     $                   A( DESCA(9), * ),
     $                   B( DESCB(9), * ),
     $                   C( DESCC(9), * )
      PARAMETER      ( ONE =  1.0D+0, ZERO =  0.0D+0 )

!     External Subroutines ..
      EXTERNAL      PDGEMM

!     If I'm not in the process grid, go to the end of the program
      IF( MYROW.EQ.-1 )
     $   GO TO 10

!     Compute product
        CALL PDGEMM( 'N', 'N', M, N, K, ONE,
     $              A, IA, JA, DESCA, 
     $              B, IB, JB, DESCB, ZERO, 
     $              C, IC, JC, DESCC )

   10 CONTINUE

      RETURN
      END
