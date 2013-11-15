! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2013, Schmidt


! "Optimal" process grid shape; tries to make the grid as close to square as 
! possible.  If a square grid is not possible, e.g. NPROCS=10, then the 
! larger number goes to NROWS, e.g. NROWS=5, NCOLS=2.
! INPUTS
  ! NPROCS = Number of processors to be used.
! OUTPUTS
  ! NROWS = "Optimal" number of process rows
  ! NCOLS = "Optimal" number of process columns
      SUBROUTINE OPTIMALGRID(NPROCS, NROWS, NCOLS)
      IMPLICIT NONE
      ! IN/OUT
      INTEGER             NPROCS, NROWS, NCOLS
      ! Local
      INTEGER             I, N
      
      
      N = INT(SQRT(REAL(NPROCS)))
      
      DO I = 0, N-1
        NCOLS = N - I
        NROWS = INT(MOD(NPROCS, NCOLS))
        IF (NROWS.EQ.0) EXIT
      END DO
      
      NROWS = NPROCS / NCOLS
      
      RETURN
      END




! Reductions using BLACS routines.  Bindings meant to mimic ALLREDUCE, which is
! much more friendly than the weird BLACS calls.
! INPUTS/OUTPUTS
  ! X = Array of data to be reduced
! INPUTS
  ! DESCX = Descriptor array of X
  ! OP = Operation: 'MIN', 'MAX', 'SUM'.  Default if you pass something else is 
    ! to sum.
  ! SCOPE = BLACS processor grid scope along which the reduction will take 
    ! place: 'All', 'Column', 'Row'
  ! RDEST/CDEST: Row/Column destination in the case of *REDUCE (rather than ALL)
      SUBROUTINE DALLREDUCE(X, DESCX, OP, SCOPE)
      IMPLICIT NONE
      ! IN/OUT
      INTEGER             DESCX(9)
      DOUBLE PRECISION    X( * )
      CHARACTER           OP, SCOPE
      ! Local
      INTEGER             M, N, LDA, ICTXT
      ! External
      EXTERNAL            DGSUM2D, DGAMX2D, DGMN2D
      
      
      M = DESCX(3)
      N = DESCX(4)
      LDA = DESCX(9)
      ICTXT = DESCX(2)
      
      IF (OP.EQ.'MIN') THEN
        CALL DGAMN2D(ICTXT,SCOPE,' ',M,N,X,LDA,-1,-1,-1,-1,-1)
      ELSE IF (OP.EQ.'MAX') THEN
        CALL DGAMX2D(ICTXT,SCOPE,' ',M,N,X,LDA,-1,-1,-1,-1,-1)
      ELSE ! default to sum
        CALL DGSUM2D(ICTXT, SCOPE, ' ', M, N, X, LDA, -1, -1)
      END IF
      
      RETURN 
      END


      SUBROUTINE DREDUCE(X, DESCX, OP, RDEST, CDEST, SCOPE)
      IMPLICIT NONE
      ! IN/OUT
      INTEGER             DESCX(9), RDEST, CDEST
      DOUBLE PRECISION    X( * )
      CHARACTER           OP, SCOPE
      ! Local
      INTEGER             M, N, LDA, ICTXT
      ! External
      EXTERNAL            DGSUM2D, DGAMX2D, DGMN2D
      
      
      M = DESCX(3)
      N = DESCX(4)
      LDA = DESCX(9)
      ICTXT = DESCX(2)
      
      IF (OP.EQ.'MIN') THEN
        CALL DGAMN2D(ICTXT, SCOPE, ' ', M, N, X, LDA, -1, -1, -1, 
     $               RDEST, CDEST)
      ELSE IF (OP.EQ.'MAX') THEN
        CALL DGAMX2D(ICTXT, SCOPE, ' ', M, N, X, LDA, -1, -1, -1,
     $               RDEST, CDEST)
      ELSE
        CALL DGSUM2D(ICTXT, SCOPE, ' ', M, N, X, LDA, RDEST, CDEST)
      END IF
      
      RETURN 
      END


      SUBROUTINE IALLREDUCE(X, DESCX, OP, SCOPE)
      IMPLICIT NONE
      ! IN/OUT
      INTEGER             DESCX(9), X( * )
      CHARACTER           OP, SCOPE
      ! Local
      INTEGER             M, N, LDA, ICTXT
      ! External
      EXTERNAL            IGSUM2D, IGAMX2D, IGMN2D
      
      
      M = DESCX(3)
      N = DESCX(4)
      LDA = DESCX(9)
      ICTXT = DESCX(2)
      
      IF (OP.EQ.'MIN') THEN
        CALL IGAMN2D(ICTXT, SCOPE, ' ', M, N, X, 1,-1,-1,-1,-1,-1)
      ELSE IF (OP.EQ.'MAX') THEN
        CALL IGAMX2D(ICTXT, SCOPE, ' ', M, N, X, 1,-1,-1,-1,-1,-1)
      ELSE
        CALL IGSUM2D(ICTXT, SCOPE, ' ', M, N, X, 1,-1,-1)
      END IF
      
      RETURN 
      END


      SUBROUTINE IREDUCE(X, DESCX, OP, RDEST, CDEST, SCOPE)
      IMPLICIT NONE
      ! IN/OUT
      INTEGER             DESCX(9), X( * ), RDEST, CDEST
      CHARACTER           OP, SCOPE
      ! Local
      INTEGER             M, N, LDA, ICTXT
      ! External
      EXTERNAL            IGSUM2D, IGAMX2D, IGMN2D
      
      
      M = DESCX(3)
      N = DESCX(4)
      LDA = DESCX(9)
      ICTXT = DESCX(2)
      
      IF (OP.EQ.'MIN') THEN
        CALL IGAMN2D(ICTXT, SCOPE, ' ', M, N, X, 1, -1, -1, -1, 
     $               RDEST, CDEST)
      ELSE IF (OP.EQ.'MAX') THEN
        CALL IGAMX2D(ICTXT, SCOPE, ' ', M, N, X, 1, -1, -1, -1,
     $               RDEST, CDEST)
      ELSE
        CALL IGSUM2D(ICTXT, SCOPE, ' ', M, N, X, 1, RDEST, CDEST)
      END IF
      
      RETURN 
      END

