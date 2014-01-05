! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2013, Schmidt


! Wrapper for pdgemr2d
! INPUTS
  ! X = Input submatrix.
  ! IX/JX = 
  ! DESCX = Descriptor array for X.
  ! IY/JY = 
  ! DESCY = Descriptor array for Y.
  ! CMNCTXT = Common BLACS context for X and Y.
! OUTPUTS
  ! Y = 
!!!      SUBROUTINE REDIST(X, IX, JX, DESCX, Y, IY, JY, DESCY, CMNCTXT)
!!!      IMPLICIT NONE
!!!      ! IN/OUT
!!!      INTEGER             IX, JX, DESCX(9), IY, JY, DESCY(9), CMNCTXT
!!!      DOUBLE PRECISION    X( * ), Y( * )
!!!      ! Local
!!!      INTEGER             M, N, MXLDM, DESCA(9),
!!!     $                    LDMX(2), LDMY(2), BLACSX(4), BLACSY(4)
!!!      ! External
!!!      EXTERNAL            PDGEMR2D
!!!      
!!!      
!!!      ! Get local and proc grid info
!!!      CALL PDIMS(DESCX, LDMX, BLACSX)
!!!      CALL PDIMS(DESCX, LDMY, BLACSY)
!!!      
!!!      M = DESCX(3)
!!!      N = DESCX(4)
!!!      
!!!      ! Adjust LDA since PDGEMR2D crashes all the time when LDA=1
!!!      DESCA(3) = 1
!!!      DESCA(4) = 1
!!!      DESCA(9) = 1
!!!      
!!!      MXLDM = MAX(LDMX)
!!!      DESCA(2) = DESCX(2)
!!!      CALL IALLREDUCE(MXLDM, DESCA, 'MAX', 'All')
!!!      IF (DESCX(9).EQ.1 .AND. DESCX(3).GT.1) DESCX(9) = MXLDM
!!!      
!!!      MXLDM = MAX(LDMY)
!!!      DESCA(2) = DESCY(2)
!!!      CALL IALLREDUCE(MXLDM, DESCA, 'MAX', 'All')
!!!      IF (DESCY(9).EQ.1 .AND. DESCY(3).GT.1) DESCY(9) = MXLDM
!!!      
!!!      ! Redistribute
!!!      CALL PDGEMR2D(M, N, X, IX, JX, DESCX,
!!!     $              Y, IY, JY, DESCY, CMNCTXT)
!!!      
!!!      RETURN
!!!      END


! Construct local submatrix from global matrix
! INPUTS
  ! GBLX = Global, non-distributed matrix.  Owned by which processor(s) depends 
    ! on R/CSRC values
  ! DESCX = ScaLAPACK descriptor array for SUBX (not a typo).
  ! RSRC/CSRC = Row/Column process value corresponding to BLACS grid for the
    ! value in DESCX(2) (the ICTXT) on which the data is stored.  If RSRC = -1,
    ! then CSRC is ignored and total ownership is assumed, i.e., GBLX is owned 
    ! by all processors.
! OUTPUTS
  ! SUBX = Local submatrix.
      SUBROUTINE MKSUBMAT(GBLX, SUBX, DESCX)!, RSRC, CSRC)
      IMPLICIT NONE
      ! IN/OUT
      INTEGER             DESCX(9), RSRC, CSRC
      DOUBLE PRECISION    GBLX(DESCX(3), DESCX(4)), SUBX(DESCX(9), *)
      ! Local
      INTEGER             M, N, I, J, GI, GJ, RBL, CBL, TI, TJ,
     $                    LDM(2), BLACS(5)
      ! External
      EXTERNAL            PDIMS, L2GPAIR
      
      
      ! Get local and proc grid info
      CALL PDIMS(DESCX, LDM, BLACS)
      
      M = LDM(1)
      N = LDM(2)
      
      RBL = DESCX(5)
      CBL = DESCX(6)
      
      IF (M.GT.0 .AND. N.GT.0) THEN
        DO J = 1, N, CBL
          DO I = 1, M, RBL
            CALL L2GPAIR(I, J, GI, GJ, DESCX, BLACS)
            
            RBL = MIN(RBL, M-I+1)
            CBL = MIN(CBL, N-J+1)
            
            DO TJ = 0, CBL-1
              !$omp do simd
              DO TI = 0, RBL-1
                SUBX(I+TI, J+TJ) = GBLX(GI+TI, GJ+TJ)
              END DO
              !$omp end do simd
            END DO
          END DO 
        END DO
      END IF
      
      RETURN
      END


! Construct global matrix from local submatrix.
! INPUTS
  ! SUBX = Local submatrix.
  ! DESCX = ScaLAPACK descriptor array for SUBX.
  ! RDEST/CDEST = Row/Column process value corresponding to BLACS grid for the
    ! value in DESCX(2) (the ICTXT) on which the global matrix GBLX will be 
    ! stored.  If RDEST = -1, then CDEST is ignored and total ownership is 
    ! assumed, i.e., GBLX is given to all processors.
! OUTPUTS
  ! GBLX = Global, non-distributed matrix.
      SUBROUTINE MKGBLMAT(GBLX, SUBX, DESCX, RDEST, CDEST)
      IMPLICIT NONE
      ! IN/OUT
      INTEGER             DESCX(9), RDEST, CDEST, PROC
      DOUBLE PRECISION    GBLX(DESCX(3), DESCX(4)), SUBX(DESCX(9), *)
      ! Local
      INTEGER             M, N, I, J, GI, GJ, RBL, CBL, TI, TJ,
     $                    LDM(2), BLACS(5)
      ! Parameter
      DOUBLE PRECISION    ZERO
      PARAMETER ( ZERO = 0.0D0 )
      ! External
      EXTERNAL            PDIMS, L2GPAIR, DGSUM2D, DALLREDUCE
      
      
      ! Get local and proc grid info
      CALL PDIMS(DESCX, LDM, BLACS)
      
      M = LDM(1)
      N = LDM(2)
      
      GBLX = ZERO
      
      RBL = DESCX(5)
      CBL = DESCX(6)
      
      IF (M.GT.0 .AND. N.GT.0) THEN
        DO J = 1, N, CBL
          DO I = 1, M, RBL
            CALL L2GPAIR(I, J, GI, GJ, DESCX, BLACS)
            
            RBL = MIN(RBL, M-I+1)
            CBL = MIN(CBL, N-J+1)
            
            DO TJ = 0, CBL-1
              !$omp do simd
              DO TI = 0, RBL-1
                GBLX(GI+TI, GJ+TJ) = SUBX(I+TI, J+TJ)
              END DO
              !$omp end do simd
            END DO
          END DO 
        END DO
      END IF
      
      IF (RDEST.EQ.-1) THEN
        CALL DALLREDUCE(GBLX, DESCX, 'S', 'All')
      ELSE
        CALL DREDUCE(GBLX, DESCX, 'S', RDEST, CDEST, 'All')
      END IF
      
      RETURN
      END

