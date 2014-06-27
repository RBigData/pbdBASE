! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2013, Schmidt


! compute inverse of a cholesky
subroutine pdchtri(uplo, x, ix, jx, descx, c, ic, jc, descc, info)
  implicit none
  ! in/out
  integer             ix, jx, descx(9), ic, jc, descc(9), info
  double precision    x(*), c(*)
  character*1         uplo
  ! local
  character*1         loup
  ! parameter
  double precision    one
  parameter ( one = 1.0d0 )
  ! external
  external            ptri2zero, pdtrtri, pdcrossprod
  
  
  if (uplo .eq. 'L') then
    loup = 'U'
  else if (uplo .eq. 'U') then
    loup = 'L'
  else 
    info = -1
    return
  end if
  
  ! zero triangle opposite uplo
  call ptri2zero(loup, 'N', x, descx)
  
  ! invert the uplo triangle
  call pdtrtri(uplo, 'N', descx(4), x, ix, jx, descx, info)
  
  ! 
  call pdcrossprod(uplo, 'T', one, x, ix, jx, descc, c, ic, jc, descc)
  
  return
end subroutine


! compute matrix inverse without having to understand scalapack peculiarities
! in place version (x is overwritten with x^-1)
subroutine pdinvip(x, ix, jx, descx, info)
  implicit none
  ! in/out
  integer             ix, jx, descx(9), info
  double precision    x(*)
  ! local
  integer             n, lwork, liwork, allocerr
  double precision    tmp
  integer, allocatable :: ipiv(:), iwork(:)
  double precision, allocatable :: work(:)
  ! external
  external           pdgetrf, pdgetri
  
  
  allocerr = 0
  n = descx(3)
  
  ! factor x=lu
  allocate(ipiv(n + descx(6)), stat=allocerr)
  if (allocerr.ne.0) return! "out of memory"
  
  call pdgetrf(n, n, x, ix, jx, descx, ipiv, info)
  if (info.ne.0) return
  
  ! invert x
  lwork = -1
  liwork = -1
  
  call pdgetri(n, x, ix, jx, descx, ipiv, tmp, lwork, liwork, liwork, info)
  if (info.ne.0) return
  
  lwork = int(tmp)
  allocate(work(lwork), stat=allocerr)
  if (allocerr.ne.0) return! "out of memory"
  
  allocate(iwork(liwork), stat=allocerr)
  if (allocerr.ne.0) return! "out of memory"
  
  call pdgetri(n, x, ix, jx, descx, ipiv, work, lwork, iwork, liwork, info)
  
  deallocate(ipiv)
  deallocate(work)
  deallocate(iwork)
  
  return
end subroutine


! non-in-place version of matrix inverse (on return, inv = x^-1)
subroutine pdinv(x, ix, jx, descx, inv, info)
  implicit none
  ! in/out
  integer             ix, jx, descx(9), info
  double precision    x(*), inv(*)
  ! external
  external           pdlacpy, pdinvip
  
  ! inv = x
  call pdlacpy('B', descx(3), descx(4), x, ix, jx, descx, inv, ix, jx, descx)
  
  call pdinvip(inv, ix, jx, descx, info)
  
  return
end subroutine


