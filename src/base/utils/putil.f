! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2013, Schmidt




! Zero out the (global) triangle and/or diagonal of a distributed matrix
! INPUTS/OUTPUTS
  ! X = Array of data to be zeroed (in part or whole)
! INPUTS
  ! UPLO = Char specifying whether 'U'pper or 'L'ower triangle is to be zeroed,
    ! or whether 'B'oth should be zeroed, or whether 'N'either triangle should.
  ! DIAG = Char specifying whether 'Y'es, zero the diagonal too, or 'N'o, do not.
! EXAMPLE:  UPLO='B', DIAG='N' produces a diagonal matrix
!           UPLO='N', DIAG='Y' zeroes out only the diagonal
!           UPLO='L', DIAG='Y' zeros the diagonal and upper triangle
      SUBROUTINE PTRI2ZERO(UPLO, DIAG, X, DESCX)
      IMPLICIT NONE
      ! IN/OUT
      INTEGER             DESCX(9)
      DOUBLE PRECISION    X(DESCX(9), *)
      CHARACTER*1         UPLO, DIAG
      ! Local
      INTEGER              M, N, I, J, GI, GJ, 
     $                    LDM(2), BLACS(5)
      ! Parameter
      DOUBLE PRECISION    ZERO
      PARAMETER ( ZERO = 0.0D0 )
      ! External
      EXTERNAL            PDIMS, L2GPAIR, PDLASET
      
      
      ! Get local and proc grid info
      CALL PDIMS(DESCX, LDM, BLACS)
      
      M = LDM(1)
      N = LDM(2)
      
      
      ! Let ScaLAPACK do the work if possible:
      IF ( (UPLO.EQ.'U' .OR. UPLO.EQ.'L') .AND. DIAG.EQ.'Y') THEN
          CALL PDLASET(UPLO, DESCX(3), DESCX(4), ZERO, ZERO, X, 
     $                 1, 1, DESCX)
      ! Only do work if we own any local pieces
      ELSE IF (M.GT.0 .AND. N.GT.0) THEN
        
        ! Diagonal only (zeroing neither triangle
        IF (UPLO.EQ.'N') THEN
          ! Quick return if possible
          IF (DIAG.EQ.'N') THEN
            RETURN
          ELSE IF (DIAG.EQ.'Y') THEN
            DO J = 1, N
              DO I = 1, M
                CALL L2GPAIR(I, J, GI, GJ, DESCX, BLACS)
                IF (GI.EQ.GJ) THEN
                  X(I,J) = ZERO
                END IF
              END DO 
            END DO
          ELSE
            RETURN! "Invalid argument 'DIAG'"
          END IF
        END IF
        
        ! Zeroing both triangles
        IF (UPLO.EQ.'B') THEN
          ! Zero the diagonal as well
          IF (DIAG.EQ.'Y') THEN
            X(:,1:N) = ZERO
          ! Leave the diagonal alone
          ELSE IF (DIAG.EQ.'N') THEN
            DO J = 1, N
              DO I = 1, M
                CALL L2GPAIR(I, J, GI, GJ, DESCX, BLACS)
                IF (GI.NE.GJ) THEN
                  X(I,J) = ZERO
                END IF
              END DO 
            END DO
          ELSE
            RETURN! "Invalid argument 'DIAG'"
          END IF
        
        ! Zeroing the upper triangle
        ELSE IF (UPLO.EQ.'U') THEN
!          ! Zero the diagonal as well
!          IF (DIAG.EQ.'Y') THEN 
!            DO J = 1, N
!              DO I = 1, M
!                CALL L2GPAIR(I, J, GI, GJ, DESCX, BLACS)
!                IF (GI.LE.GJ) THEN
!                  X(I,J) = ZERO
!                END IF
!              END DO
!            END DO
          ! Leave the diagonal alone
          IF (DIAG.EQ.'N') THEN 
            DO J = 1, N
              DO I = 1, M
                CALL L2GPAIR(I, J, GI, GJ, DESCX, BLACS)
                IF (GI.LT.GJ) THEN
                  X(I,J) = ZERO
                END IF
              END DO
            END DO
          ELSE
            RETURN! "Invalid argument 'DIAG'"
          END IF
        
        ! Zeroing the lower triangle
        ELSE IF (UPLO.EQ.'L') THEN
!          ! Zero the diagonal as well
!          IF (DIAG.EQ.'Y') THEN 
!            DO J = 1, N
!              DO I = 1, M
!                CALL L2GPAIR(I, J, GI, GJ, DESCX, BLACS)
!                IF (GI.GE.GJ) THEN
!                  X(I,J) = ZERO
!                END IF
!              END DO
!            END DO
          ! Leave the diagonal alone
          IF (DIAG.EQ.'N') THEN 
            DO J = 1, N
              DO I = 1, M
                CALL L2GPAIR(I, J, GI, GJ, DESCX, BLACS)
                IF (GI.GT.GJ) THEN
                  X(I,J) = ZERO
                END IF
              END DO
            END DO
          ELSE
            RETURN! "Invalid argument 'DIAG'"
          END IF
        ELSE
          RETURN! "Invalid argument 'UPLO'"
        END IF
      
      END IF
      
      RETURN
      END SUBROUTINE


! Make the matrix symmetric via copying from one (global) triangle to the other.
! INPUTS/OUTPUTS
  ! X = Array of data to be zeroed (in part or whole)
! INPUTS
  ! UPLO = Copy FROM 'U'pper or copy FROM 'L'ower.  The modified data lives in 
    ! the opposite triangle from UPLO.
  ! IX/JX = 
  ! DESCX = ScaLAPACK descriptor array for X.
      SUBROUTINE PDMKSYM(UPLO, X, IX, JX, DESCX)
      IMPLICIT NONE
      ! IN/OUT
      INTEGER             IX, JX, DESCX(9)
      DOUBLE PRECISION    X(DESCX(9), *)
      CHARACTER*1         UPLO
      ! Local
      INTEGER             M, N, I, J, GI, GJ, LDM(2), BLACS(5)
      DOUBLE PRECISION, ALLOCATABLE :: CPX(:,:)
      ! Parameter
      DOUBLE PRECISION    ZERO, ONE, HALF
      PARAMETER ( ZERO = 0.0D0, ONE = 1.0D0, HALF = 0.5D0 )
      ! External
      EXTERNAL            PDIMS, L2GPAIR, PDGEADD, PTRI2ZERO
      
      
      ! Get local and proc grid info
      CALL PDIMS(DESCX, LDM, BLACS)
      
      M = LDM(1)
      N = LDM(2)
      
      ! Only do work if we own any local pieces
      IF (M.GT.0 .AND. N.GT.0) THEN
        
        ! Zero Lower
        IF (UPLO.EQ.'U') THEN
          CALL PTRI2ZERO('L', 'N', X, DESCX)
        ! Zero Upper
        ELSE IF (UPLO.EQ.'L') THEN
          CALL PTRI2ZERO('U', 'N', X, DESCX)
        ! Otherwise...
        ELSE IF (UPLO.NE.'B') THEN
          RETURN! "Invalid argument 'UPLO'"
        END IF
        
        ALLOCATE(CPX(M,N))
        CPX(1:M,1:N) = X(1:M,1:N)
        
        ! X = t(X) + X
        CALL PDGEADD('T', DESCX(3), DESCX(4), ONE, CPX, IX, JX, 
     $                DESCX, ONE, X, IX, JX, DESCX)
        
        DEALLOCATE(CPX)
        
        ! Correct diagonal
        DO J = 1, N
          DO I = 1, M
            CALL L2GPAIR(I, J, GI, GJ, DESCX, BLACS)
            IF (GI.EQ.GJ) THEN
              X(I,J) = X(I,J) * HALF
            END IF
          END DO 
        END DO
        
      END IF
      
      RETURN
      END SUBROUTINE


!!! Does what pdlacpy SHOULD do, namely copies the triangle of a global 
!!! distributed matrix onto the triangle of another matrix.
!!!
!!! If TRANS = 'Y', then A^T will be used.  This requires communication.
!!!
!!! For simplicity, we assume DESCA = DESCB.
!!! Not as efficient as it could be, but them's the breaks.
!!!
!!! PDMKSYM should be more efficient if A=B.
!!      SUBROUTINE PDGELACPY2(TRANS, UPLO, DIAG, A, IA, JA, DESCA, 
!!     $                      B, IB, JB, INFO)
!!      IMPLICIT NONE
!!      ! IN/OUT
!!      INTEGER             IA, JA, DESCA(9), IB, JB, INFO
!!      DOUBLE PRECISION    A(DESCA(9), *), B(DESCA(9), * )
!!      CHARACTER*1         TRANS, UPLO, DIAG
!!      ! Local
!!      INTEGER             M, N, I, J, GI, GJ, LDIM(2), BLACS(5)
!!      DOUBLE PRECISION, ALLOCATABLE :: CPX(:,:)
!!      ! Parameter
!!      DOUBLE PRECISION    ZERO, ONE
!!      PARAMETER ( ZERO = 0.0D0, ONE = 1.0D0 )
!!      ! External
!!      EXTERNAL            PDIMS, PDLACPY, PTRI2ZERO, PDGEADD
!!      
!!      
!!      INFO = 0
!!      
!!      ! Get local and proc grid info
!!      CALL PDIMS(DESCA, LDIM, BLACS)
!!      
!!      M = LDIM(1)
!!      N = LDIM(2)
!!      
!!      ! Upper part to be copied
!!      IF (TRANS.EQ.'N') THEN
!!        IF (UPLO.EQ.'U') THEN
!!          ! Copy over the diagonal if needed
!!          IF (DIAG.EQ.'Y') THEN
!!            DO J = 1, N
!!              DO I = 1, M
!!                CALL L2GPAIR(I, J, GI, GJ, DESCA, BLACS)
!!                IF (GI.GE.GJ) THEN
!!                  B(I,J) = A(I,J)
!!                END IF
!!              END DO 
!!            END DO
!!          ELSE IF (DIAG.EQ.'N') THEN
!!            DO J = 1, N
!!              DO I = 1, M
!!                CALL L2GPAIR(I, J, GI, GJ, DESCA, BLACS)
!!                IF (GI.GT.GJ) THEN
!!                  B(I,J) = A(I,J)
!!                END IF
!!              END DO 
!!            END DO
!!          ELSE
!!            INFO = -3 ! Invalid argument 'DIAG'
!!            RETURN
!!          END IF
!!      
!!      ELSE IF (TRANS.EQ.'Y') THEN
!!        IF (UPLO.EQ.'U' .OR. UPLO.EQ.'L') THEN
!!          ! Zero out to-be unused parts of B
!!          CALL PTRI2ZERO(UPLO, DIAG, B, DESCA) 
!!          
!!          ! Copy over via PDGEADD
!!          CALL PDGEADD('N', DESCA(3), DESCA(4), ONE, A, IA, JA, 
!!     $                  DESCA, ONE, A, IA, JA, DESCA)
!!        ELSE
!!          
!!        END IF
!!      ELSE
!!        INFO = -1 ! Invalid argument 'TRANS'
!!        RETURN 
!!      END IF
!!      
!!      RETURN
!!      END SUBROUTINE



! Grab diagonal of matrix
! INPUTS
  ! X = Submatrix whose global diagonal is to be taken
  ! IX/JX = 
  ! DESCX = Descriptor array for X
! OUTPUTS
  ! DIAG = Diagonal of the global matrix for which X is the submatrix.
      SUBROUTINE PDGDGTK(X, IX, JX, DESCX, DIAG, RDEST, CDEST)
      IMPLICIT NONE
      ! IN/OUT
      INTEGER             IX, JX, DESCX(9), RDEST, CDEST
      DOUBLE PRECISION    X(DESCX(9), *), DIAG( * )
!      CHARACTER*1         REDUCE
      ! Local
      INTEGER             K, M, N, I, J, GI, GJ, LDM(2), 
     $                    BLACS(5)
      ! Parameter
      DOUBLE PRECISION    ZERO
      PARAMETER ( ZERO = 0.0D0 )
      ! External
      EXTERNAL            PDIMS, DALLREDUCE
      
      
      ! Get local and proc grid info
      CALL PDIMS(DESCX, LDM, BLACS)
      
      M = LDM(1)
      N = LDM(2)
      
      K = MIN(DESCX(3), DESCX(4))
      
      DIAG(1:K) = ZERO
      
      DO J = 1, N
        DO I = 1, M
          CALL L2GPAIR(I, J, GI, GJ, DESCX, BLACS)
          IF (GI.EQ.GJ) THEN
            DIAG(GI) = X(I,J)
          END IF
        END DO
      END DO
      
      CALL DGSUM2D(DESCX(2), 'All', ' ', K, 1, DIAG, K, RDEST, CDEST)
      
      RETURN
      END SUBROUTINE


! Construct matrix containing diagonal DIAG
! INPUTS
  ! IX/JX = 
  ! DESCX = Descriptor array for X
  ! DIAG = Global vector which should be the diagonal of a distributed matrix
  ! LDIAG = Length of Diag.  Need not be equal to K = MIN(DESCX(3), DESCX(4)),
    ! though only the first K elements of DIAG will be used.
! OUTPUTS
  ! X = Submatrix of the global matrix with diagonal DIAG
      SUBROUTINE PDDIAGMK(X, IX, JX, DESCX, DIAG, LDIAG)
      IMPLICIT NONE
      ! IN/OUT
      INTEGER             IX, JX, DESCX(9), LDIAG
      DOUBLE PRECISION    X(DESCX(9), *), DIAG( * )
      ! Local
      INTEGER             M, N, I, J, GI, GJ, LDM(2), BLACS(5)
      ! Parameter
      DOUBLE PRECISION    ZERO
      PARAMETER ( ZERO = 0.0D0 )
      ! External
      EXTERNAL            PDIMS
      ! Function
      INTEGER             IND
      
      
      ! Get local and proc grid info
      CALL PDIMS(DESCX, LDM, BLACS)
      
      M = LDM(1)
      N = LDM(2)
      
      DO J = 1, N
        DO I = 1, M
          CALL L2GPAIR(I, J, GI, GJ, DESCX, BLACS)
          IF (GI.EQ.GJ) THEN
            X(I,J) = DIAG( IND(GI, LDIAG) )
          ELSE 
            X(I,J) = ZERO
          END IF
        END DO
      END DO
      
      RETURN
      END SUBROUTINE

