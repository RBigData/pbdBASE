! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2013, 2016, Schmidt


! pdlange subroutine wrapper for use with f77_call in c
! also handles allocation of work vector
subroutine matnorm(value, norm, m, n, a, ia, ja, desca)
  implicit none
  ! in/out
  character           norm
  integer             m, n, ia, ja, desca( 9 )
  double precision    value, a( * )
  ! local
  integer             lwork, bl, ioffa, iamar, nprow, npcol, myprow, mypcol
  ! dynamic
  double precision, allocatable :: work(:)
  ! functions
  double precision    pdlange
  integer             indxg2p, numroc
  
  
  if (norm .eq. "M" .or. norm .eq. "F") then
    lwork = 0
  else 
    call blacs_gridinfo(desca(2), nprow, npcol, myprow, mypcol)
    
    if (norm .eq. "O" .or. norm .eq. "1") then
      bl = desca(5)
      ioffa = mod(ja-1, bl)
      
      iamar = indxg2p(ia, bl, mypcol, 0, npcol)
      
      lwork = numroc(n+ioffa, bl, mypcol, 0, npcol)
    else if (norm .eq. "I") then
      bl = desca(6)
      ioffa = mod(ia-1, bl)
      
      iamar = indxg2p(ia, bl, myprow, 0, nprow)
      
      lwork = numroc(m+ioffa, bl, myprow, 0, nprow)
    end if
  end if
  
  allocate(work(lwork))
  value = pdlange(norm, m, n, a, ia, ja, desca, work)
  deallocate(work)
  
  return
end subroutine



!     condition number estimator for general matrix
!       step 1:  get matrix norm of a
!       step 2:  factor a=lu
!       step 3:  call pdgecon
subroutine condnum(norm, m, n, a, ia, ja, desca, rcond, info)
  implicit none
  ! in/out
  character           norm
  integer             m, n, ldim1, ia, ja, desca( 9 ), info
  double precision    rcond, a( * )
  ! local
  integer             lwork, liwork, lipiv
  double precision    anorm, tmp
  ! dynamic
  double precision, allocatable :: work(:)
  integer, allocatable :: iwork(:), ipiv(:)
  ! functions
  external            matnorm, pdgetrf, pdgecon
  
  
  ! step 1:  get matrix norm of a
  call matnorm(anorm, norm, n, n, a, ia, ja, desca)
  
  ! step 2:  factor a=lu
  lipiv = desca(9) + desca(5) ! locr(m_a)+ mb_a
  allocate(ipiv(lipiv))
  
  call pdgetrf(m, n, a, ia, ja, desca, ipiv, info)
  
  if (info .ne. 0) then
    deallocate(ipiv)
    return
  end if
  
  ! step 3:  call pdgecon
  call pdgecon(norm, n, a, ia, ja, desca, anorm, rcond, tmp, -1, liwork, -1, info)
  
  lwork = int(tmp)
  
  allocate(work(lwork))
  allocate(iwork(liwork))
  
  call pdgecon(norm, n, a, ia, ja, desca, anorm, rcond, work, lwork, iwork, liwork, info)
  
  deallocate(ipiv)
  deallocate(iwork)
  deallocate(work)
  
  return
end subroutine
