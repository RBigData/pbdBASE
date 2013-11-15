! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2013, Schmidt



! In the case of Matrix-Vector operations where the vector is global and not
! necessarily of "appropriate" length, this is used to adjust which element
! of the vector is used with the matrix.  Enables the equivalent of doing,
! for example, something like matrix(1, nrow=5, ncol=2) + 1:3 in R.
! INPUTS
  ! I = Index
  ! M = Modulus
      INTEGER FUNCTION IND(I, M)
      IMPLICIT NONE
      INTEGER             I, M
      
      
      IND = MOD(I, M)
      IF (IND.EQ.0) THEN
        IND = M
      END IF
      
      RETURN
      END


! Convert matrix indexing to vector indexing, or vice versa
!      SUBROUTINE INDMAT2VEC(STORAGE, M, N, I, J, K)
!      IMPLICIT NONE
!      CHARACTER*1         STORAGE
!      INTEGER             M, N, I, J, K
!      ! Function
!      INTEGER             IND
!      
!      
!      IF (STORAGE .EQ. 'C') THEN
!        
!      ELSE IF (STORAGE .EQ. 'R') THEN
!        
!      ELSE
!        IND = -1
!      END IF
!      
!      RETURN
!      END


!      SUBROUTINE INDVEC2MAT(STORAGE, M, N, I, J, K)
!      IMPLICIT NONE
!      CHARACTER*1         STORAGE
!      INTEGER             M, N, I, J, K
!      ! Function
!      INTEGER             IND
!      
!      
!      IF (STORAGE .EQ. 'C') THEN
!        I = IND(K, M)
!        J = K/M + 1
!      ELSE IF (STORAGE .EQ. 'R') THEN
!        I = IND(K, N)
!        J = K/N + 1
!      ELSE
!        I = -1
!        J = -1
!      END IF
!      
!      RETURN
!      END



! Takes a ScaLAPACK descriptor array and from it determines full (1) local
! dimension information (calling NUMROC) and (2) BLACS grid information.
! The local dimension LDIM is set to (/ 0, 0 /) if there is not actual 
! ownership of local data on that process.
! INPUTS
  ! DESC = ScaLAPACK descriptor array
! OUTPUTS
  ! LDM = [LDM1, LDM2] = Local dimension
  ! BLACS = [NPROCS, NPROW, NPCOL, MYPROW, MYPCOL]
      SUBROUTINE PDIMS(DESC, LDM, BLACS)
      IMPLICIT NONE
      ! IN/OUT
      INTEGER             DESC(9), LDM(2), BLACS(5)
      ! Functions
      INTEGER             NUMROC
      ! External
      EXTERNAL            BLACS_GRIDINFO
      
      
      CALL BLACS_GRIDINFO(DESC(2), BLACS(2), BLACS(3), 
     $                    BLACS(4), BLACS(5))
      
      IF (BLACS(2).EQ.(-1) .OR. BLACS(3).EQ.(-1)) THEN
        BLACS(1) = -1
      ELSE
        BLACS(1) = BLACS(2) * BLACS(3)
      END IF
      
      LDM(1) = NUMROC(DESC(3), DESC(5), BLACS(4), 
     $                DESC(7), BLACS(2))
      LDM(2) = NUMROC(DESC(4), DESC(6), BLACS(5), 
     $                DESC(8), BLACS(3))
      
      IF (LDM(1).LT.1 .OR. LDM(2).LT.1) THEN
        LDM(1) = 0
        LDM(2) = 0
      END IF
      
      RETURN
      END


! Local-to-global pair of indices; shorthand for calling INDXL2G twice.
! INPUTS
  ! I/J = Local coordinates.
  ! DESC = BLACS descriptor array.
  ! BLACS = BLACS process grid information, taken from PDIMS.
! OUTPUTS
  ! GI/GJ = Global coordinates.
      SUBROUTINE L2GPAIR(I, J, GI, GJ, DESC, BLACS)
      IMPLICIT NONE
      ! IN/OUT
      INTEGER             I, J, GI, GJ, DESC(9), BLACS(2)
      ! Functions
      INTEGER             INDXL2G
      
      
      GI = INDXL2G(I, DESC(5), BLACS(4), 0, BLACS(2))
      GJ = INDXL2G(J, DESC(6), BLACS(5), 0, BLACS(3))
      
      RETURN
      END




! Global-to-local pair of indices; shorthand for calling INDXG2L twice.
! INPUTS
  ! GI/GJ = Global coordinates.
  ! DESC = BLACS descriptor array.
  ! BLACS = BLACS process grid information, taken from PDIMS.
! OUTPUTS
  ! I/J = Local coordinates.
      SUBROUTINE G2LPAIR(I, J, GI, GJ, DESC, BLACS)
      IMPLICIT NONE
      ! IN/OUT
      INTEGER             I, J, GI, GJ, DESC(9), BLACS(2)
      ! Local
      INTEGER             DUM
      ! Functions
      INTEGER             INDXG2L
      
      
      I = INDXG2L(GI, DESC(5), DUM, DUM, BLACS(2))
      J = INDXG2L(GJ, DESC(6), DUM, DUM, BLACS(3))
      
      RETURN
      END
