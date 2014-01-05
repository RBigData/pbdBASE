! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2013, Schmidt


! R level 2 BLAS
! INPUTS/OUTPUTS
  ! X = Submatrix of data which should globally be "swept"
! INPUTS
  ! IX/JX = 
  ! DESCX = Descriptor array for X
  ! VEC = Vector to "sweep" through X
  ! LVEC = Length of VEC
  ! FUN = Char with 4 possibilities, describing the type of sweep to perform:
    ! "+", "-", "*", "/"
      DOUBLE PRECISION FUNCTION FPMOD(A, B)
      DOUBLE PRECISION A, B
      
      IF (B.EQ.0) THEN
        FPMOD = NaN
      ELSE IF (B.GT.0) THEN
        IF (A.GE.0) THEN
          FPMOD = DMOD(A, B)
        ELSE
          FPMOD = B - DMOD(-A, B)
        END IF
      ELSE
        IF (A.EQ.0) THEN
          FPMOD = 0.0D0
        ELSE IF (A.GT.0) THEN
          FPMOD = B + DMOD(A, -B)
        ELSE
          FPMOD = -DMOD(-A, -B)
        END IF
      END IF
      
      RETURN
      END



      SUBROUTINE RL2BLAS(X, IX, JX, DESCX, VEC, LVEC, FUN)
      IMPLICIT NONE
      ! IN/OUT
      INTEGER             IX, JX, DESCX(9), LVEC, FUN
      DOUBLE PRECISION    X(DESCX(9), *), VEC(LVEC)
      ! Local
      INTEGER             K, M, N, POS, I, J, GI, GJ, LDM(2), BLACS(5)
      INTEGER             L
      ! Parameter
      DOUBLE PRECISION    ZERO, ONE
      PARAMETER ( ZERO = 0.0D0, ONE = 1.0D0 )
      ! External
      EXTERNAL            PDIMS, L2GPAIR
      ! Function
      INTEGER             IND
      DOUBLE PRECISION    FPMOD
      
      
      ! Get local and proc grid info
      CALL PDIMS(DESCX, LDM, BLACS)
      
      M = LDM(1)
      N = LDM(2)
      
      K = DESCX(3)
      
      
      ! Resorting to magic numbers because C strings are just too kludgy and terrible
      ! This is all very ad hoc anyway so I don't think I give a shit.
      
      ! Only do work if we own any local pieces
      IF (M.GT.0 .AND. N.GT.0) THEN
        ! Addition
        IF (FUN.EQ.0) THEN
          DO J = 1, N
            DO I = 1, M
              CALL L2GPAIR(I, J, GI, GJ, DESCX, BLACS)
              POS = IND(GI + K*(GJ-1), LVEC)
              X(I, J) = X(I, J) + VEC(POS)
            END DO
          END DO
        ! Subtraction
        ELSE IF (FUN.EQ.1) THEN
          DO J = 1, N
            DO I = 1, M
              CALL L2GPAIR(I, J, GI, GJ, DESCX, BLACS)
              POS = IND(GI + K*(GJ-1), LVEC)
              X(I, J) = X(I, J) - VEC(POS)
            END DO
          END DO
        ! Multiplication
        ELSE IF (FUN.EQ.2) THEN
          DO J = 1, N
            DO I = 1, M
              CALL L2GPAIR(I, J, GI, GJ, DESCX, BLACS)
              POS = IND(GI + K*(GJ-1), LVEC)
              X(I, J) = X(I, J) * VEC(POS)
            END DO
          END DO
        ! Division
        ELSE IF (FUN.EQ.3) THEN
          DO J = 1, N
            DO I = 1, M
              CALL L2GPAIR(I, J, GI, GJ, DESCX, BLACS)
              POS = IND(GI + K*(GJ-1), LVEC)
              X(I, J) = X(I, J) / VEC(POS)
            END DO
          END DO
        ! Power
        ELSE IF (FUN.EQ.4) THEN
          DO J = 1, N
            DO I = 1, M
              CALL L2GPAIR(I, J, GI, GJ, DESCX, BLACS)
              POS = IND(GI + K*(GJ-1), LVEC)
              X(I, J) = X(I, J) ** VEC(POS)
            END DO
          END DO
        ! %%
        ELSE IF (FUN.EQ.5) THEN
          DO J = 1, N
            DO I = 1, M
              CALL L2GPAIR(I, J, GI, GJ, DESCX, BLACS)
              POS = IND(GI + K*(GJ-1), LVEC)
              X(I, J) = FPMOD( X(I, J), VEC(POS) )
            END DO
          END DO 
        ! %/%
        ELSE IF (FUN.EQ.6) THEN
          DO J = 1, N
            DO I = 1, M
              CALL L2GPAIR(I, J, GI, GJ, DESCX, BLACS)
              POS = IND(GI + K*(GJ-1), LVEC)
              X(I, J) = FPMOD( VEC(POS), X(I, J) )
            END DO
          END DO
        ! <
        ELSE IF (FUN.EQ.7) THEN
          DO J = 1, N
            DO I = 1, M
              CALL L2GPAIR(I, J, GI, GJ, DESCX, BLACS)
              POS = IND(GI + K*(GJ-1), LVEC)
              IF (X(I, J) .LT. VEC(POS)) THEN
                X(I, J) = ONE
              ELSE
                X(I, J) = ZERO
              END IF
            END DO
          END DO
        ! >
        ELSE IF (FUN.EQ.8) THEN
          DO J = 1, N
            DO I = 1, M
              CALL L2GPAIR(I, J, GI, GJ, DESCX, BLACS)
              POS = IND(GI + K*(GJ-1), LVEC)
              IF (X(I, J) .GT. VEC(POS)) THEN
                X(I, J) = ONE
              ELSE
                X(I, J) = ZERO
              END IF
            END DO
          END DO
        ! <=
        ELSE IF (FUN.EQ.9) THEN
          DO J = 1, N
            DO I = 1, M
              CALL L2GPAIR(I, J, GI, GJ, DESCX, BLACS)
              POS = IND(GI + K*(GJ-1), LVEC)
              IF (X(I, J) .LE. VEC(POS)) THEN
                X(I, J) = ONE
              ELSE
                X(I, J) = ZERO
              END IF
            END DO
          END DO
        ! >=
        ELSE IF (FUN.EQ.10) THEN
          DO J = 1, N
            DO I = 1, M
              CALL L2GPAIR(I, J, GI, GJ, DESCX, BLACS)
              POS = IND(GI + K*(GJ-1), LVEC)
              IF (X(I, J) .GE. VEC(POS)) THEN
                X(I, J) = ONE
              ELSE
                X(I, J) = ZERO
              END IF
            END DO
          END DO
        ! <
        ELSE IF (FUN.EQ.11) THEN
          DO J = 1, N
            DO I = 1, M
              CALL L2GPAIR(I, J, GI, GJ, DESCX, BLACS)
              POS = IND(GI + K*(GJ-1), LVEC)
              IF (X(I, J) .EQ. VEC(POS)) THEN
                X(I, J) = ONE
              ELSE
                X(I, J) = ZERO
              END IF
            END DO
          END DO
        END IF
      END IF
      
      RETURN
      END



      LOGICAL FUNCTION CHECKPROC(I, J, DESC, BLACS)
      INTEGER I, J, DESC(9), BLACS(5)
      
      CHECKPROC = ( MOD( (I-1)/DESC(5), BLACS(2) ) .EQ. BLACS(4)
     $                      .AND.
     $              MOD( (J-1)/DESC(6), BLACS(3) ) .EQ. BLACS(5) )
      
      RETURN
      END


! R-style replacement.  (dist)Matrix-Vector and (dist)Matrix-(dist)Matrix (levels
! 2 and 3) are available.
! INPUTS/OUTPUTS
  ! X = Submatrix of data which should be replaced
! INPUTS
  ! IX/JX = 
  ! DESCX = Descriptor array for X
  ! LVEC = Length of VEC
      SUBROUTINE RL2INSERT(X, IX, JX, DESCX, VEC, LVEC, INDI, LINDI, 
     $                   INDJ, LINDJ)
      IMPLICIT NONE
      ! IN/OUT
      INTEGER             IX, JX, DESCX(9), LVEC, LINDI, LINDJ,
     $                    INDI(LINDI), INDJ(LINDJ)
      DOUBLE PRECISION    X(DESCX(9), *), VEC(LVEC)
      ! Local
      INTEGER             K, M, N, POS, I, J, TI, TJ, GI, GJ, 
     $                    LDM(2), BLACS(5)
      ! Parameter
      DOUBLE PRECISION    ZERO, ONE
      PARAMETER ( ZERO = 0.0D0, ONE = 1.0D0 )
      ! External
      EXTERNAL            PDIMS, G2LPAIR
      ! Function
      LOGICAL             CHECKPROC
      INTEGER             IND
      DOUBLE PRECISION    FPMOD
      
      
      ! Get local and proc grid info
      CALL PDIMS(DESCX, LDM, BLACS)
      
      M = LDM(1)
      N = LDM(2)
      K = DESCX(3)
      
      ! Only do work if we own any local pieces
      IF (M.GT.0 .AND. N.GT.0) THEN
        ! Insertion
        DO TJ = 1, LINDJ
          GJ = INDJ(TJ)
          DO TI = 1, LINDI
            GI = INDI(TI)
            CALL G2LPAIR(I, J, GI, GJ, DESCX, BLACS)
            IF (CHECKPROC(GI, GJ, DESCX, BLACS)) THEN
              POS = IND(GI + K*(GJ-1), LVEC)
              X(I, J) = VEC(POS)
            END IF
          END DO
        END DO
      END IF
      
      RETURN
      END



! Column YCOL of Y is copied onto column XCOL of X.  
! LCOLS = number of entries in XCOLS and YCOLS
      SUBROUTINE RCOLCPY(X, DESCX, XCOLS, Y, DESCY, YCOLS, LCOLS)
      IMPLICIT NONE
      ! IN/OUT
      INTEGER             DESCX(9), DESCY(9), LCOLS, XCOLS(LCOLS), 
     $                    YCOLS(LCOLS)
      DOUBLE PRECISION    X(DESCX(9), *), Y(DESCY(9), *)
      ! Local
      LOGICAL             IHAVE, INEED
      INTEGER             I, J, GI, GJ, MX, NX, MY, NY, GM, GN,
     $                    RBL, CBL, LDM(2), BLACS(5), 
     $                    LXCOL, LYCOL, XCOL, YCOL, COL,
     $                    RSRC, CSRC, RDEST, CDEST,
     $                    RXLEN, RYLEN, RXTOP, RYTOP
      ! External
      EXTERNAL            PDIMS, INDXG2L, INDXL2G
      ! Function
      INTEGER             INDXG2L, INDXL2G
      ! ################################################################
      
      
      GM = DESCX(3)
      GN = DESCX(4)
      
      IF (GM.NE.DESCY(3)) RETURN
      
      RBL = DESCX(5)
      CBL = DESCX(6)
      
      ! Get local and proc grid info
      CALL PDIMS(DESCY, LDM, BLACS)
      MY = LDM(1)
      NY = LDM(2)
      
      CALL PDIMS(DESCX, LDM, BLACS)
      
      MX = LDM(1)
      NX = LDM(2)
      
      DO COL = 1, LCOLS, 1
        XCOL = XCOLS(COL)
        YCOL = YCOLS(COL)
        
        LXCOL = INDXG2L(XCOL, DESCX(6), 1, 1, BLACS(3))
        LYCOL = INDXG2L(YCOL, DESCY(6), 1, 1, BLACS(3))
        DO GI = 1, GM, RBL
          ! Index juggling
          I = INDXG2L(GI, RBL, 1, 1, BLACS(2))
          
          RXTOP = MIN(I+RBL-1, MX)
          RXLEN = RXTOP - I + 1
          
          RYTOP = MIN(I+RBL-1, MY)
          RYLEN = RYTOP - I + 1
          
          ! Row and column (processor) source
          RSRC = MOD( (GI-1)/RBL, BLACS(2) )
          CSRC = MOD( (YCOL-1)/CBL, BLACS(3) )
          
          IHAVE = ( RSRC.EQ.BLACS(4) .AND. CSRC.EQ.BLACS(5) )
          
          ! Row and column (processor) destination
          RDEST = MOD( (GI-1)/RBL, BLACS(2) )
          CDEST = MOD( (XCOL-1)/CBL, BLACS(3) )
          
          INEED = ( RDEST.EQ.BLACS(4) .AND. CDEST.EQ.BLACS(5) )
          
          ! Copy
          IF (IHAVE) THEN ! Check if need to SEND
            IF (INEED) THEN ! Easy case
              X(I:RXTOP, LXCOL) = Y(I:RYTOP, LYCOL)
            ELSE ! Otherwise SEND
              ! Send
              CALL DGESD2D(DESCX(2), RYLEN, 1, Y(I, LYCOL), RYLEN, 
     $                     RDEST, CDEST)
            END IF
          ELSE IF (INEED) THEN ! Otherwise check if need to RECEIVE
            ! Receive
            CALL DGERV2D(DESCX(2), RXLEN, 1, X(I, LXCOL), RXLEN, 
     $                   RSRC, CSRC)
          END IF
          
          CALL BLACS_BARRIER(DESCX(2), 'A') 
          
        END DO ! GI
      END DO ! COL
      
      RETURN
      END



! Column YCOL of Y is copied onto column XCOL of X. 
! LXCOLS is allowed to be an integer multiple of LYCOLS
      SUBROUTINE RCOLCPY2(X, DESCX, XCOLS, LXCOLS, Y, DESCY, YCOLS, 
     $                    LYCOLS)
      IMPLICIT NONE
      ! IN/OUT
      INTEGER             DESCX(9), DESCY(9), LXCOLS, XCOLS(LXCOLS),
     $                    LYCOLS, YCOLS(LYCOLS)
      DOUBLE PRECISION    X(DESCX(9), *), Y(DESCY(9), *)
      ! Local
      INTEGER             COL
      ! External
      EXTERNAL            RCOLCPY
      ! ################################################################
      
      
      IF (MOD(LXCOLS, LYCOLS).NE.0) THEN
        RETURN
      END IF
      
      DO COL = 1, LXCOLS, LYCOLS
        CALL RCOLCPY(X, DESCX, XCOLS(COL), Y, DESCY, YCOLS, LYCOLS)
      END DO
      
      RETURN
      END


! Rows YROWS of Y are copied onto rows XROWS of X.  
! Obvious assumptions are made.
      SUBROUTINE RROWCPY(X, DESCX, XROWS, Y, DESCY, YROWS, LROWS)
      IMPLICIT NONE
      ! IN/OUT
      INTEGER             DESCX(9), DESCY(9), LROWS, XROWS(LROWS), 
     $                    YROWS(LROWS)
      DOUBLE PRECISION    X(DESCX(9), *), Y(DESCY(9), *)
      ! Local
      LOGICAL             IHAVE, INEED
      INTEGER             I, J, GI, GJ, MX, NX, MY, NY, GM, GN,
     $                    RBL, CBL, LDM(2), BLACS(5), 
     $                    LXROW, LYROW, XROW, YROW, ROW,
     $                    RSRC, CSRC, RDEST, CDEST,
     $                    CXLEN, CYLEN, CXTOP, CYTOP
      ! External
      EXTERNAL            PDIMS, INDXG2L, INDXL2G
      ! Function
      INTEGER             INDXG2L, INDXL2G
      ! ################################################################
      
      
      GM = DESCX(3)
      GN = DESCX(4)
      
      IF (GN.NE.DESCY(4)) RETURN
      
      RBL = DESCX(5)
      CBL = DESCX(6)
      
      ! Get local and proc grid info
      CALL PDIMS(DESCY, LDM, BLACS)
      MY = LDM(1)
      NY = LDM(2)
      
      CALL PDIMS(DESCX, LDM, BLACS)
      
      MX = LDM(1)
      NX = LDM(2)
      
      ! The work
      DO ROW = 1, LROWS, 1
        XROW = XROWS(ROW)
        YROW = YROWS(ROW)
        
        LXROW = INDXG2L(XROW, DESCX(6), 1, 1, BLACS(3))
        LYROW = INDXG2L(YROW, DESCY(6), 1, 1, BLACS(3))
        
!        J = INDXG2L(GJ, CBL, 1, 1, BLACS(3))
        
        DO GJ = 1, GN, CBL
          ! Index juggling
          J = INDXG2L(GJ, CBL, 1, 1, BLACS(2))
          
          CXTOP = MIN(J+CBL-1, NX)
          CXLEN = CXTOP - J + 1
          
          CYTOP = MIN(J+CBL-1, NY)
          CYLEN = CYTOP - J + 1
          
          ! Row and column (processor) source
          RSRC = MOD( (YROW-1)/RBL, BLACS(2) )
          CSRC = MOD( (GJ-1)/CBL, BLACS(3) )
          
          IHAVE = ( RSRC.EQ.BLACS(4) .AND. CSRC.EQ.BLACS(5) )
          
          ! Row and column (processor) destination
          RDEST = MOD( (XROW-1)/RBL, BLACS(2) )
          CDEST = MOD( (GJ-1)/CBL, BLACS(3) )
          
          INEED = ( RDEST.EQ.BLACS(4) .AND. CDEST.EQ.BLACS(5) )
          
          ! Copy
          IF (IHAVE) THEN ! Check if need to SEND
            IF (INEED) THEN ! Easy case
              X(LXROW, J:CXTOP) = Y(LYROW, J:CYTOP)
            ELSE ! Otherwise SEND
              ! Send
              CALL DGESD2D(DESCX(2), 1, CYLEN, Y(LYROW, J), CYLEN, 
     $                     RDEST, CDEST)
            END IF
          ELSE IF (INEED) THEN ! Otherwise check if need to RECEIVE
            ! Receive
            CALL DGERV2D(DESCX(2), 1, CXLEN, X(LXROW, J), CXLEN, 
     $                   RSRC, CSRC)
          END IF
          
          CALL BLACS_BARRIER(DESCX(2), 'A') 
          
        END DO ! GI
      END DO ! ROW
      
      RETURN
      END

! Rows YROWS of Y are copied onto rows XROWS of X.  
! LXROWS is allowed to be an integer multiple of LYROWS
      SUBROUTINE RROWCPY2(X, DESCX, XROWS, LXROWS, Y, DESCY, YROWS, 
     $                    LYROWS)
      IMPLICIT NONE
      ! IN/OUT
      INTEGER             DESCX(9), DESCY(9), LXROWS, XROWS(LXROWS),
     $                    LYROWS, YROWS(LYROWS)
      DOUBLE PRECISION    X(DESCX(9), *), Y(DESCY(9), *)
      ! Local
      INTEGER             ROW
      ! External
      EXTERNAL            RROWCPY
      ! ################################################################
      
      
      IF (MOD(LXROWS, LYROWS).NE.0) THEN
        RETURN
      END IF
      
      DO ROW = 1, LXROWS, LYROWS
        CALL RROWCPY(X, DESCX, XROWS(ROW), Y, DESCY, YROWS, LYROWS)
      END DO
      
      RETURN
      END










! Distributed matrix-distributed vector sum (like R)
! X is an M_1 x N distributed matrix, Y is a distributed M_2 x 1 vector
! X <-- X + Y
      SUBROUTINE PDMVSUM(X, DESCX, Y, DESCY)
      IMPLICIT NONE
      ! IN/OUT
      INTEGER             DESCX(9), DESCY(9)
      DOUBLE PRECISION    X(DESCX(9), *), Y(DESCY(9))
      ! Local
      INTEGER             LDM(2), BLACS(5),
     $                    I, J, GI, GJ, II, JJ, POS,
     $                    GM, GN, RBL, CBL, MY, NY, 
     $                    MX, NX, LVEC
      DOUBLE PRECISION    TMP
      
      
      ! Get local and proc grid info
      GM = DESCX(3)
      GN = DESCX(4)
      
      LVEC = DESCY(5)
      
      RBL = DESCX(5)
      CBL = DESCX(6)
      
      CALL PDIMS(DESCY, LDM, BLACS)
      MY = LDM(1)
      NY = LDM(2)
      
      CALL PDIMS(DESCX, LDM, BLACS)
      MX = LDM(1)
      NX = LDM(2)
      
      ! Quick return if possible
      IF (DESCY(4) .NE. 1) RETURN
      
      
      ! Easy case:  nrows(x) == length(y)
      IF (GM .EQ. DESCY(3)) THEN
        DO J = 1, NX
          DO I = 1, MX
            X(I, J) = X(I, J) + Y(I)
          END DO
        END DO
        
      ! General case
!!!      ELSE
!!!        
!!!        DO J = 1, NX!, CBL
!!!          DO I = 1, MX!, RBL
!!!            CALL L2GPAIR(I, J, GI, GJ, DESC, BLACS)
!!!            POS = IND(GI + GN*(GJ-1), LVEC)
!!!            
!!!  !          X(I, J) = X(I, J) + VEC(POS)
!!!            
!!!            
!!!            
!!!            I = INDXG2L(GI, RBL, 1, 1, BLACS(2))
!!!            
!!!            RXTOP = MIN(I+RBL-1, MX)
!!!            RXLEN = RXTOP - I + 1
!!!            
!!!            RYTOP = MIN(I+RBL-1, MY)
!!!            RYLEN = RYTOP - I + 1
!!!            
!!!            ! Row and column (processor) source
!!!            RSRC = MOD( (GI-1)/RBL, BLACS(2) )
!!!            CSRC = MOD( (YCOL-1)/CBL, BLACS(3) )
!!!            
!!!            IHAVE = ( RSRC.EQ.BLACS(4) .AND. CSRC.EQ.BLACS(5) )
!!!            
!!!            ! Row and column (processor) destination
!!!            RDEST = MOD( (GI-1)/RBL, BLACS(2) )
!!!            CDEST = MOD( (XCOL-1)/CBL, BLACS(3) )
!!!            
!!!            INEED = ( RDEST.EQ.BLACS(4) .AND. CDEST.EQ.BLACS(5) )
!!!            
!!!            ! Copy
!!!            IF (IHAVE) THEN ! Check if need to SEND
!!!              IF (INEED) THEN ! Easy case
!!!                X(I:RXTOP, LXCOL) = Y(I:RYTOP, LYCOL)
!!!              ELSE ! Otherwise SEND
!!!                ! Send
!!!                CALL DGESD2D(DESCX(2), RYLEN, 1, Y(I, LYCOL), RYLEN, 
!!!     $                       RDEST, CDEST)
!!!              END IF
!!!            ELSE IF (INEED) THEN ! Otherwise check if need to RECEIVE
!!!              ! Receive
!!!              CALL DGERV2D(DESCX(2), RXLEN, 1, X(I, LXCOL), RXLEN, 
!!!     $                     RSRC, CSRC)
!!!            END IF
!!!          
!!!            CALL BLACS_BARRIER(DESCX(2), 'A') 
!!!            
!!!            
!!!            
!!!            
!!!            
!!!            
!!!          END DO
!!!        END DO
!!!      
!!!      
!!!      
      
      END IF
      
      
      
      RETURN
      END






