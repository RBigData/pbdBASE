! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2013-2014, Schmidt


! x^t * x or x * x^t
! trans = 't' :  x^t*x, trans = 'n' : x*x^t
subroutine pdcrossprod(uplo, trans, alpha, x, ix, jx, descx, c, ic, jc, descc)
  implicit none
  ! in/out
  integer             ix, jx, descx(9), ic, jc, descc(9)
  double precision    x(*), c(*), alpha
  character*1         uplo, trans
  ! local
  integer             ldx, ldc
  character*1         nst
  ! external
  external            pdsyrk, pdmksym
  
  
  if (trans .eq. 'T') then
    nst = 'N'
    ldx = descx(4)
    ldc = descx(3)
  else 
    nst = 'T'
    ldx = descx(3)
    ldc = descx(4)
  end if
  
  ! compute 'uplo' triangle of x^t*x or x^t*x
  call pdsyrk(uplo, nst, ldc, ldx, alpha, x, ix, jx, descx, 0.0d0, c, ic, jc, descc)
  
  ! fill lower triangle (make symmetric)
  call pdmksym(uplo, c, ic, jc, descc)
  
  return
end subroutine

