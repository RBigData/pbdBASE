!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     PDGETRF:  Compute LU Factorization of matrix A
!     
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      SUBROUTINE RPDGETRF(A, ICTXT, MYROW, MYCOL,
     $                   DESCA, M, N,
     $                   LIPIV, INFO)

      INTEGER            IA, JA, M, N, LIPIV,
     $                   ICTXT, INFO, MYROW, MYCOL,
     $                   DESCA( 9 ), 
     $                   IPIV ( LIPIV )
      PARAMETER        ( IA = 1, JA = 1 )

      DOUBLE PRECISION   A( DESCA(9), * )

      EXTERNAL           PDGETRF

!     If I'm not in the process grid, go to the end of the program
      IF( MYROW.EQ.-1 )
     $   GO TO 10

!     Compute LU decomposition
      CALL PDGETRF( M, N, A, IA, JA, DESCA,
     $              IPIV, INFO )

   10 CONTINUE

      RETURN
      END
