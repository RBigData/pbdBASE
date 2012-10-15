!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Selection Algorithm, taken from Numerical Recipes in Fortran 77.
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     INPUTS:
!       k = desired rank
!       n = length of vector
!       VEC = vector from which rank k element is to be chosen
!     OUTPUT:
!       RET = rank k element of VEC
      SUBROUTINE QKSLCT(k, n, VEC, RET)
      ! Inputs
      INTEGER k, n
      DOUBLE PRECISION VEC(n)
      ! Return
      DOUBLE PRECISION RET
      ! Local
      INTEGER i, ir, j, l, mid
      REAL A, TEMP

      l = 1
      ir = n
 2    IF (ir-l.LE.1) THEN
         IF (ir-1.EQ.1) THEN
            IF (VEC(ir).LT.VEC(l)) THEN
               TEMP = VEC(l)
               VEC(l) = VEC(ir)
               VEC(ir) = TEMP
            ENDIF
         ENDIF
         RET = VEC(k)
         RETURN
      ELSE
         mid = (l+ir)/2
         TEMP = VEC(mid)
         VEC(mid) = VEC(l+1)
         VEC(l+1) = TEMP
         IF (VEC(l).GT.VEC(ir)) THEN
            TEMP = VEC(l)
            VEC(l) = VEC(ir)
            VEC(ir) = TEMP
         ENDIF
         IF (VEC(l+1).GT.VEC(ir)) THEN
            TEMP = VEC(l+1)
            VEC(l+1) = VEC(ir)
            VEC(ir) = TEMP
         ENDIF
         IF (VEC(l).GT.VEC(l+1)) THEN
            TEMP = VEC(l)
            VEC(l) = VEC(l+1)
            VEC(l+1) = TEMP
         ENDIF
         i = l+1
         j = ir
         A = VEC(l+1)
 3       CONTINUE
         i = i+1
         IF (VEC(i).LT.A) GOTO 3
 4       CONTINUE
         j = j-1
         IF (VEC(j).GT.A) GOTO 4
         IF (j.LT.i) GOTO 5
         TEMP = VEC(i)
         VEC(i) = VEC(j)
         VEC(j) = TEMP
         GOTO 3
 5       VEC(l+1) = VEC(j)
         VEC(j) = A
         IF (j.GE.k) ir = j-1
         IF (j.LE.k) l = i
      ENDIF
      GOTO 2
      END

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Quick median calculator, using Quick select
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     INPUTS:
!       n = length of vector
!       VEC = vector from which median is to be chosen
!     OUTPUT:
!       MED = median of array VEC
      SUBROUTINE QKMED(n, VEC, MED)
      ! Input
      INTEGER n
      DOUBLE PRECISION VEC( n )
      ! Return
      DOUBLE PRECISION MED, TMP
      ! Local
      INTEGER odd, k1, k2
      DOUBLE PRECISION neginf
      
      EXTERNAL QKSLCT
      
      ! Need -Inf
      neginf = -huge(1.)
      
      MED = 0
      TMP = 0
      
      odd = MODULO(n, 2)
      k1 = n/2
      k2 = k1+1
      
      IF (n.LT.1) THEN
        MED = neginf
      ELSE IF (odd.EQ.1) THEN 
        CALL QKSLCT(k2, n, VEC, MED)
      ELSE
        CALL QKSLCT(k1, n, VEC, MED)
        CALL QKSLCT(k2, n, VEC, TMP)
        MED = (MED + TMP) / 2
      ENDIF
      
      RETURN
      END

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Partition function, based on Quicksort Algorithm
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     INPUTS:
!       VEC = vector from which median is to be chosen
!       l = leftmost index
!       r = rightmost index
!       pvt = pivot
!     OUTPUT:
!       ret = number of elements in VEC <= pvt
      SUBROUTINE PTTN(VEC, l, r, pvt, ret)
      ! Input
      INTEGER l, r
      DOUBLE PRECISION pvt, VEC( * )
      ! Return
      INTEGER ret
      ! Local
      DOUBLE PRECISION tmp
      INTEGER i
      
      ret = l
      DO 1 i = l, r
        IF (VEC(i).LT.pvt) THEN
          tmp = VEC(ret)
          VEC(ret) = VEC(i)
          VEC(i) = tmp
          ret = ret+1
        END IF
 1    CONTINUE
      
      ret = ret-1
      
      RETURN
      END
