!
! Not yet working; heavily under construction
!
!!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!     Parallel selection
!!     Based on: http://www.umiacs.umd.edu/research/EXPAR/papers/3494/node18.html#SECTION00051000000000000000
!!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!     INPUTS:
!!       k = desired rank
!!       n = length of vector
!!       VEC = vector from which rank k element is to be chosen
!!     OUTPUT:
!!       ret = rank k element of VEC

!      RECURSIVE SUBROUTINE PQKSLCT(ICTXT, NPROCS, MYPNUM, 
!     $                             k, VEC, l, r, n, ret)
!      ! Inputs
!      INTEGER ICTXT, NPROCS, MYPNUM, k, l, r, n, ct
!      DOUBLE PRECISION VEC(r)
!      ! Output
!      DOUBLE PRECISION ret
!      ! Local
!      INTEGER test, lbelow, lmeds, nmeds
!      DOUBLE PRECISION MED, mxvec, closest, allclosest
!      DOUBLE PRECISION MEDS(NPROCS)
!      ! Parameters
!      DOUBLE PRECISION inf, neginf, neginff
!      PARAMETER        ( inf = huge(1.),
!     $                   neginf = -huge(1.), 
!     $                   neginff = neginf / 10. )
!      
!      EXTERNAL PTTN, QKMED, qsort3, 
!     $         IGSUM2D, DGSUM2D, DGAMN2D, BLACS_BARRIER
!      
!      ret = neginf
!      
!      ! Local median
!      CALL QKMED(n, VEC, MED)

!      ! Create vector of medians, find median of medians
!      DO 10 i = 1, NPROCS
!        MEDS(i) = 0.
! 10   CONTINUE

!      MEDS(MYPNUM) = MED
!      CALL BLACS_BARRIER( ICTXT, 'A' )
!      CALL DGSUM2D( ICTXT, 'All', ' ', NPROCS, 1, MEDS, 
!     $              1, -1, -1 )

!      CALL PTTN(MEDS, 1, NPROCS, neginff, lmeds)

!      nmeds = NPROCS-lmeds
!      lmeds = lmeds+1

!      CALL QKMED(nmeds, MEDS(lmeds), MED)

!      ! Partition VEC by median of medians, sum up test values
!      CALL PTTN(VEC, l, r, MED, test)
!      lbelow = test

!      CALL BLACS_BARRIER( ICTXT, 'A' )
!      CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, test, 1, -1, -1 )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!      call system( 'sleep .5' )
!!      WRITE (*,*) "k=",k
!!      WRITE (*,*) "n=",n
!!      WRITE (*,*) "lbelow=",lbelow
!!      WRITE (*,*) "test=",test
!!      CALL BLACS_BARRIER( ICTXT, 'A' )
!!      WRITE (*,*) "------------------------------------"
!!      CALL BLACS_BARRIER( ICTXT, 'A' )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!      WRITE (*,*) "k=",k
!!      WRITE (*,*) "test=",test
!!      WRITE (*,*) "VEC=",VEC

!      ! tie-breaker
!!WCC  WRITE (*,*) "k=",k
!!WCC  WRITE (*,*) "test=",test
!!WCC  WRITE (*,*) "VEC=",VEC
!      
!      IF (test.EQ.0) THEN
!!        WRITE (*,*) "MED=",MED
!!        ct = 0
!!        DO 99 i = 1, SIZE(VEC)
!!          IF (VEC(1).EQ.VEC(i)) THEN
!!            ct = ct + 1
!!          END IF
!! 99     CONTINUE
!!        IF (ct.EQ.SIZE(VEC)) THEN
!!          ret = VEC(1)
!!          RETURN
!        ret = MED 
!        RETURN
!!        END IF
!      END IF
!      
!      ! otherwise...
!      IF (test.LT.k) THEN
!!        WRITE (*,*) "test < k"
!        
!        lbelow = lbelow+1
!        k = k-test
!        n = n-lbelow+1
!        CALL BLACS_BARRIER( ICTXT, 'A' )
!        CALL PQKSLCT(ICTXT, NPROCS, MYPNUM, k, VEC(lbelow), 
!     $               1, n, n, ret)
!      ELSE IF (test.GT.k) THEN
!!        WRITE (*,*) "test > k"
!        
!        n = lbelow
!        CALL BLACS_BARRIER( ICTXT, 'A' )
!        CALL PQKSLCT(ICTXT, NPROCS, MYPNUM, k, VEC, 
!     $               1, n, n, ret)
!      ELSE
!!        WRITE (*,*) "test == k"

!!WCC    WRITE (*,*) "----------------------------"
!        CALL BLACS_BARRIER( ICTXT, 'A' )

!        IF (SIZE(VEC).EQ.0) THEN
!          closest = neginf
!        ELSE
!          ret = neginf
!          DO 11 i = 1, SIZE(VEC)
!            IF (VEC(i).LT.MED) THEN
!              IF (VEC(i).GT.ret) THEN
!                ret = VEC(i)
!              END IF
!            END IF
! 11       CONTINUE
!          closest = MED - ret
!        END IF

!!        WRITE (*,*) "ret=",ret
!        
!        allclosest = closest
!        CALL DGAMN2D( ICTXT, 'A', ' ', 1, 1, allclosest, 1, 
!     $                -1, -1, -1, -1, -1 )
!        
!        IF (allclosest.NE.closest) THEN
!          ret = inf
!        END IF

!        CALL DGAMN2D( ICTXT, 'All', ' ', 1, 1, ret, 1,
!     $              -1, -1, -1, -1, -1 )

!!WCC    WRITE (*,*) "VEC=",VEC
!!WCC    WRITE (*,*) "test=",test
!      END IF

!      CALL BLACS_BARRIER( ICTXT, 'A' )
!      RETURN
!      END
