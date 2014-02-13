! this source code form is subject to the terms of the mozilla public
! license, v. 2.0. if a copy of the mpl was not distributed with this
! file, you can obtain one at http://mozilla.org/mpl/2.0/.

! copyright 2013, schmidt


! subroutine wrapper for numroc()
subroutine numrocwrap(n, nb, iproc, nprocs, num)
  ! in/out
  integer             n, nb, iproc, nprocs, num
  ! functions
  integer             numroc
  
  
  num = numroc(n, nb, iproc, 0, nprocs)
  
  return
end subroutine



! in the case of matrix-vector operations where the vector is global and not
! necessarily of "appropriate" length, this is used to adjust which element
! of the vector is used with the matrix.  enables the equivalent of doing,
! for example, something like matrix(1, nrow=5, ncol=2) + 1:3 in r.
! inputs
  ! i = index
  ! m = modulus
integer function ind(i, m)
  implicit none
  integer             i, m
  
  
  ind = mod(i, m)
  if (ind.eq.0) then
    ind = m
  end if
  
  return
end function


! convert matrix indexing to vector indexing, or vice versa
!subroutine indmat2vec(storage, m, n, i, j, k)
!  implicit none
!  character*1         storage
!  integer             m, n, i, j, k
!  ! function
!  integer             ind
!  
!  
!  if (storage .eq. 'c') then
!    
!  else if (storage .eq. 'r') then
!    
!  else
!    ind = -1
!  end if
!  
!  return
!end subroutine


!subroutine indvec2mat(storage, m, n, i, j, k)
!  implicit none
!  character*1         storage
!  integer             m, n, i, j, k
!  ! function
!  integer             ind
!  
!  
!  if (storage .eq. 'c') then
!    i = ind(k, m)
!    j = k/m + 1
!  else if (storage .eq. 'r') then
!    i = ind(k, n)
!    j = k/n + 1
!  else
!    i = -1
!    j = -1
!  end if
!  
!  return
!end subroutine



! takes a scalapack descriptor array and from it determines full (1) local
! dimension information (calling numroc) and (2) blacs grid information.
! the local dimension ldim is set to (/ 0, 0 /) if there is not actual 
! ownership of local data on that process.
! inputs
  ! desc = scalapack descriptor array
! outputs
  ! ldm = [ldm1, ldm2] = local dimension
  ! blacs = [nprocs, nprow, npcol, myprow, mypcol]
subroutine pdims(desc, ldm, blacs)
  implicit none
  ! in/out
  integer             desc(9), ldm(2), blacs(5)
  ! functions
  integer             numroc
  ! external
  external            blacs_gridinfo
  
  
  call blacs_gridinfo(desc(2), blacs(2), blacs(3), blacs(4), blacs(5))
  
  if (blacs(2).eq.(-1) .or. blacs(3).eq.(-1)) then
    blacs(1) = -1
  else
    blacs(1) = blacs(2) * blacs(3)
  end if
  
  ldm(1) = numroc(desc(3), desc(5), blacs(4), desc(7), blacs(2))
  ldm(2) = numroc(desc(4), desc(6), blacs(5), desc(8), blacs(3))
  
  if (ldm(1).lt.1 .or. ldm(2).lt.1) then
    ldm(1) = 0
    ldm(2) = 0
  end if
  
  return
end subroutine


! local-to-global pair of indices; shorthand for calling indxl2g twice.
! inputs
  ! i/j = local coordinates.
  ! desc = blacs descriptor array.
  ! blacs = blacs process grid information, taken from pdims.
! outputs
  ! gi/gj = global coordinates.
subroutine l2gpair(i, j, gi, gj, desc, blacs)
  implicit none
  ! in/out
  integer             i, j, gi, gj, desc(9), blacs(5)
  ! functions
  integer             indxl2g
  
  
  gi = indxl2g(i, desc(5), blacs(4), 0, blacs(2))
  gj = indxl2g(j, desc(6), blacs(5), 0, blacs(3))
  
  return
end subroutine




! global-to-local pair of indices; shorthand for calling indxg2l twice.
! inputs
  ! gi/gj = global coordinates.
  ! desc = blacs descriptor array.
  ! blacs = blacs process grid information, taken from pdims.
! outputs
  ! i/j = local coordinates.
subroutine g2lpair(i, j, gi, gj, desc, blacs)
  implicit none
  ! in/out
  integer             i, j, gi, gj, desc(9), blacs(5)
  ! local
  integer             dum
  ! functions
  integer             indxg2l
  
  
  i = indxg2l(gi, desc(5), dum, dum, blacs(2))
  j = indxg2l(gj, desc(6), dum, dum, blacs(3))
  
  return
end subroutine

