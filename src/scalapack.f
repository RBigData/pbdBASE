!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     PDGESV:  Solve system AX=B
!     X --> B
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      SUBROUTINE RPDGESV(A, B, ICTXT, MYROW, MYCOL,
     $                   DESCA, DESCB, N, NRHS,
     $                   MXLDIMS, INFO)

      INTEGER            IA, JA, IB, JB, N, NRHS, 
     $                   ICTXT, INFO, MYROW, MYCOL,
     $                   MXRHSC, MXLDIMS(4),
     $                   DESCA( 9 ), DESCB( 9 ),
     $                   IPIV ( MXLDIMS(1) + DESCA(6) )
      PARAMETER        ( IA = 1, JA = 1, IB = 1, JB = 1 )

      DOUBLE PRECISION   A( DESCA(9), * ), 
     $                   B( DESCB(9), MXLDIMS(4) )

      EXTERNAL           PDGESV

!     If I'm not in the process grid, go to the end of the program
      IF( MYROW.EQ.-1 )
     $   GO TO 10

!     Solve the linear system A * X = B
      CALL PDGESV( N, NRHS, 
     $             A, IA, JA, DESCA, IPIV, 
     $             B, IB, JB, DESCB, INFO )

   10 CONTINUE

      RETURN
      END
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     PDGESVD:  Singular value decomposition
!     X --> B
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      ! First determine the size of LWORK.  Combining with the 
      ! subroutine below causes bizarre memory deallocation errors.
      SUBROUTINE RPDGESVDSZ(M, N, ASIZE, ICTXT, MYROW, MYCOL,
     $                    DESCA, DESCU, DESCVT, TEMP, INFO, 
     $                    JOBU, JOBVT)

      INTEGER          :: LWORK = -1
      INTEGER             M, N, ASIZE,
     $                    DESCA( 9 ), DESCU( 9 ), DESCVT( 9 ),
     $                    ICTXT, INFO, MYROW, MYCOL,
     $                    IA, JA, IU, JU, IVT, JVT
      PARAMETER         ( IA = 1, JA = 1, IU = 1, JU = 1,
     $                    IVT = 1, JVT = 1 )

      DOUBLE PRECISION    A( 1, 1 ), 
     $                    D( 1 ),
     $                    U( 1,1 ),
     $                    VT( 1,1 )
      DOUBLE PRECISION    TEMP( 1 )
      
      CHARACTER*1         JOBU, JOBVT

      EXTERNAL            PDGESVD

!     If I'm not in the process grid, go to the end of the program
      IF( MYROW.EQ.-1 )
     $   GO TO 10

      CALL PDGESVD(JOBU, JOBVT, M, N, 
     $            A, IA, JA, DESCA, 
     $            D, 
     $            U, IU, JU, DESCU, 
     $            VT, IVT, JVT, DESCVT, 
     $            TEMP, LWORK, INFO)

   10 CONTINUE

      RETURN
      END
      !
      ! Calculation of SVD
      !
      SUBROUTINE RPDGESVD(M, N, ASIZE, ICTXT, MYROW, MYCOL,
     $                    A, DESCA, D,
     $                    U, DESCU, VT, DESCVT,
     $                    INFO, LWORK, JOBU, JOBVT)

      INTEGER             M, N, ASIZE,
     $                    DESCA( 9 ), DESCU( 9 ), DESCVT( 9 ),
     $                    ICTXT, INFO, MYROW, MYCOL,
     $                    IA, JA, IU, JU, IVT, JVT
      PARAMETER         ( IA = 1, JA = 1, IU = 1, JU = 1,
     $                    IVT = 1, JVT = 1 )

      DOUBLE PRECISION    A( DESCA(9), * ), 
     $                    D( ASIZE ),
     $                    U( DESCU(9), * ),
     $                    VT( DESCVT(9), * )
      DOUBLE PRECISION, ALLOCATABLE :: WORK(:)

      CHARACTER*1         JOBU, JOBVT

      EXTERNAL            PDGESVD

!     If I'm not in the process grid, go to the end of the program
      IF( MYROW.EQ.-1 )
     $   GO TO 10

      ALLOCATE (WORK(LWORK))
      
      CALL PDGESVD(JOBU, JOBVT, M, N, 
     $            A, IA, JA, DESCA, 
     $            D, 
     $            U, IU, JU, DESCU, 
     $            VT, IVT, JVT, DESCVT, 
     $            WORK, LWORK, INFO)

   10 CONTINUE

      DEALLOCATE (WORK)

      RETURN
      END
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     PDGETRI:  Invert matrix A
!     inv(A) --> A
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      ! First calculate LWORK and LIWORK
      SUBROUTINE RPDGETRISZ( ICTXT, MYROW, MYCOL, DESCA, N, INFO,
     $                       TEMP, ITEMP )

      INTEGER              :: LWORK = -1
      INTEGER              :: LIWORK = -1
      INTEGER            IA, JA, N
     $                   ICTXT, INFO, MYROW, MYCOL,
     $                   DESCA( 9 ),
     $                   IPIV ( N + DESCA(6) ),
     $                   ITEMP( 1 )
      PARAMETER        ( IA = 1, JA = 1 )

      DOUBLE PRECISION   A( 1,1 ), 
     $                   TEMP( 1 )

      EXTERNAL           PDGETRI

!     If I'm not in the process grid, go to the end of the program
      IF( MYROW.EQ.-1 )
     $   GO TO 10

!     Determine LWORK and LIWORK
      CALL PDGETRI( N, A, IA, JA, DESCA,
     $              IPIV, TEMP, LWORK, ITEMP, LIWORK,
     $              INFO )

   10 CONTINUE

      RETURN
      END
      !
      ! Calculation of Inverse
      !
      SUBROUTINE RPDGETRI(A, ICTXT, MYROW, MYCOL,
     $                   DESCA, N,
     $                   INFO, LWORK, LIWORK)

      INTEGER, ALLOCATABLE :: IWORK(:)
      INTEGER            IA, JA, N, LWORK, LIWORK,
     $                   ICTXT, INFO, MYROW, MYCOL,
     $                   DESCA( 9 ),
     $                   IPIV ( N + DESCA(6) ),
     $                   LTEMP( 1 )
      PARAMETER        ( IA = 1, JA = 1 )

      DOUBLE PRECISION   A( DESCA(9), * ), 
     $                   TEMP( 1 )
      DOUBLE PRECISION, ALLOCATABLE :: WORK(:)

      EXTERNAL           PDGETRI, PDGETRF

!     If I'm not in the process grid, go to the end of the program
      IF( MYROW.EQ.-1 )
     $   GO TO 10

!     Compute LU decomposition
      CALL PDGETRF( N, N, A, IA, JA, DESCA,
     $              IPIV, INFO )

      ALLOCATE (WORK(LWORK))
      ALLOCATE (IWORK(LIWORK))

      CALL PDGETRI( N, A, IA, JA, DESCA,
     $              IPIV, WORK, LWORK, IWORK, LIWORK,
     $              INFO )

      DEALLOCATE(WORK)
      DEALLOCATE(IWORK)

   10 CONTINUE

      RETURN
      END
!
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
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     PDPOTRF:  Compute Cholesky factorization of symmetric, 
!             positive definite, distributed matrix A
!     
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      SUBROUTINE RPDPOTRF( A, ICTXT, MYROW, MYCOL,
     $                   DESCA, N,
     $                   UPLO, INFO)

      INTEGER            IA, JA, N, INFO,
     $                   ICTXT, MYROW, MYCOL,
     $                   DESCA( 9 )
      PARAMETER        ( IA = 1, JA = 1 )

      DOUBLE PRECISION   A( DESCA(9), * )

      CHARACTER*1         UPLO

      EXTERNAL           PDPOTRF

!     If I'm not in the process grid, go to the end of the program
      IF( MYROW.EQ.-1 )
     $   GO TO 10

!     Compute LU decomposition
      CALL PDPOTRF( UPLO, N,
     $            A, IA, JA, DESCA,
     $            INFO )

   10 CONTINUE

      RETURN
      END
!!
!!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!     PDGEQRF:  QR Factorization
!!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!      SUBROUTINE RPDGEQRF( A, ICTXT, MYROW, MYCOL,
!     $                   DESCA, N,
!     $                   UPLO, INFO)

!      INTEGER            IA, JA, N, INFO,
!     $                   ICTXT, MYROW, MYCOL,
!     $                   DESCA( 9 )
!      PARAMETER        ( IA = 1, JA = 1 )

!      DOUBLE PRECISION   A( DESCA(9), * )

!      CHARACTER*1         UPLO

!      EXTERNAL           PDGEQRF

!!     If I'm not in the process grid, go to the end of the program
!      IF( MYROW.EQ.-1 )
!     $   GO TO 10

!!     Compute QR factorization
!      CALL PDGEQRF( m, n, a, 
!                    ia, ja, desca, 
!                    tau, work, lwork, info)

!   10 CONTINUE

!      RETURN
!      END

