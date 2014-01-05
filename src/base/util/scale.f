! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2013, Schmidt


! SWEEP array out of distributed matrix
! INPUTS/OUTPUTS
  ! X = Submatrix of data which should globally be "swept"
! INPUTS
  ! IX/JX = 
  ! DESCX = Descriptor array for X
  ! VEC = Vector to "sweep" through X
  ! LVEC = Length of VEC
  ! MARGIN = 1 for row sweeping, 2 for column sweeping
  ! FUN = Char with 4 possibilities, describing the type of sweep to perform:
    ! "+", "-", "*", "/"
      SUBROUTINE PDSWEEP(X, IX, JX, DESCX, VEC, LVEC, MARGIN, FUN)
      IMPLICIT NONE
      ! IN/OUT
      INTEGER             IX, JX, DESCX(9), MARGIN, LVEC
      DOUBLE PRECISION    X(DESCX(9), *), VEC(LVEC)
      CHARACTER*1         FUN
      ! Local
      INTEGER             K, M, N, POS, I, J, GI, GJ, LDM(2), BLACS(5)
      ! External
      EXTERNAL            PDIMS
      ! Function
      INTEGER             IND
      
      
      ! Get local and proc grid info
      CALL PDIMS(DESCX, LDM, BLACS)
      
      M = LDM(1)
      N = LDM(2)
      
      ! Only do work if we own any local pieces
      IF (M.GT.0 .AND. N.GT.0) THEN
        ! Addition
        IF (FUN.EQ."+") THEN
          IF (MARGIN.EQ.1) THEN
            K = DESCX(3)
            DO J = 1, N
              DO I = 1, M
                CALL L2GPAIR(I, J, GI, GJ, DESCX, BLACS)
                POS = IND(GI + K*(GJ-1), LVEC)
                X(I, J) = X(I, J) + VEC(POS)
              END DO
            END DO
          ELSE IF (MARGIN.EQ.2) THEN
            K = DESCX(4)
            DO J = 1, N
              DO I = 1, M
                CALL L2GPAIR(I, J, GI, GJ, DESCX, BLACS)
                POS = IND(GJ + K*(GI-1), LVEC)
                X(I, J) = X(I, J) + VEC(POS)
              END DO
            END DO
          END IF
        ! Subtraction
        ELSE IF (FUN.EQ."-") THEN
          IF (MARGIN.EQ.1) THEN
            K = DESCX(3)
            DO J = 1, N
              DO I = 1, M
                CALL L2GPAIR(I, J, GI, GJ, DESCX, BLACS)
                POS = IND(GI + K*(GJ-1), LVEC)
                X(I, J) = X(I, J) - VEC(POS)
              END DO
            END DO
          ELSE IF (MARGIN.EQ.2) THEN
            K = DESCX(4)
            DO J = 1, N
              DO I = 1, M
                CALL L2GPAIR(I, J, GI, GJ, DESCX, BLACS)
                POS = IND(GJ + K*(GI-1), LVEC)
                X(I, J) = X(I, J) - VEC(POS)
              END DO
            END DO
          END IF
        ! Multiplication
        ELSE IF (FUN.EQ."*") THEN
          IF (MARGIN.EQ.1) THEN
            K = DESCX(3)
            DO J = 1, N
              DO I = 1, M
                CALL L2GPAIR(I, J, GI, GJ, DESCX, BLACS)
                POS = IND(GI + K*(GJ-1), LVEC)
                X(I, J) = X(I, J) * VEC(POS)
              END DO
            END DO
          ELSE IF (MARGIN.EQ.2) THEN
            K = DESCX(4)
            DO J = 1, N
              DO I = 1, M
                CALL L2GPAIR(I, J, GI, GJ, DESCX, BLACS)
                POS = IND(GJ + K*(GI-1), LVEC)
                X(I, J) = X(I, J) * VEC(POS)
              END DO
            END DO
          END IF
        ! Division
        ELSE IF (FUN.EQ."/") THEN
          IF (MARGIN.EQ.1) THEN
            K = DESCX(3)
            DO J = 1, N
              DO I = 1, M
                CALL L2GPAIR(I, J, GI, GJ, DESCX, BLACS)
                POS = IND(GI + K*(GJ-1), LVEC)
                X(I, J) = X(I, J) / VEC(POS)
              END DO
            END DO
          ELSE IF (MARGIN.EQ.2) THEN
            K = DESCX(4)
            DO J = 1, N
              DO I = 1, M
                CALL L2GPAIR(I, J, GI, GJ, DESCX, BLACS)
                POS = IND(GJ + K*(GI-1), LVEC)
                X(I, J) = X(I, J) / VEC(POS)
              END DO
            END DO
          END IF
        END IF
      END IF
      
      RETURN
      END

