! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2013, Schmidt


! r level 2 blas
! inputs/outputs
  ! x = submatrix of data which should globally be "swept"
! inputs
  ! ix/jx = 
  ! descx = descriptor array for x
  ! vec = vector to "sweep" through x
  ! lvec = length of vec
  ! fun = char with 4 possibilities, describing the type of sweep to perform:
    ! "+", "-", "*", "/"


subroutine rl2blas(x, ix, jx, descx, vec, lvec, fun)
  implicit none
  ! in/out
  integer :: ix, jx, descx(9), lvec, fun
  double precision :: x(descx(9), *), vec(lvec)
  ! local
  integer :: k, m, n, pos, i, j, gi, gj, ldm(2), blacs(5)
  integer :: l
  ! parameter
  double precision, parameter :: zero = 0.0d0, one = 1.0d0
  ! external
  external            pdims, l2gpair
  ! function
  integer :: ind
  double precision :: fpmod
  
  
  ! get local and proc grid info
  call pdims(descx, ldm, blacs)
  
  m = ldm(1)
  n = ldm(2)
  
  k = descx(3)
  
  
  ! resorting to magic numbers because c strings are just too kludgy and terrible
  ! this is all very ad hoc anyway so i don't think i give a shit.
  
  ! only do work if we own any local pieces
  if (m.gt.0 .and. n.gt.0) then
    ! addition
    if (fun.eq.0) then
      do j = 1, n
        do i = 1, m
          call l2gpair(i, j, gi, gj, descx, blacs)
          pos = ind(gi + k*(gj-1), lvec)
          x(i, j) = x(i, j) + vec(pos)
        end do
      end do
    ! subtraction
    else if (fun.eq.1) then
      do j = 1, n
        do i = 1, m
          call l2gpair(i, j, gi, gj, descx, blacs)
          pos = ind(gi + k*(gj-1), lvec)
          x(i, j) = x(i, j) - vec(pos)
        end do
      end do
    ! multiplication
    else if (fun.eq.2) then
      do j = 1, n
        do i = 1, m
          call l2gpair(i, j, gi, gj, descx, blacs)
          pos = ind(gi + k*(gj-1), lvec)
          x(i, j) = x(i, j) * vec(pos)
        end do
      end do
    ! division
    else if (fun.eq.3) then
      do j = 1, n
        do i = 1, m
          call l2gpair(i, j, gi, gj, descx, blacs)
          pos = ind(gi + k*(gj-1), lvec)
          x(i, j) = x(i, j) / vec(pos)
        end do
      end do
    ! power
    else if (fun.eq.4) then
      do j = 1, n
        do i = 1, m
          call l2gpair(i, j, gi, gj, descx, blacs)
          pos = ind(gi + k*(gj-1), lvec)
          x(i, j) = x(i, j) ** vec(pos)
        end do
      end do
    ! %%
    else if (fun.eq.5) then
      do j = 1, n
        do i = 1, m
          call l2gpair(i, j, gi, gj, descx, blacs)
          pos = ind(gi + k*(gj-1), lvec)
          x(i, j) = fpmod( x(i, j), vec(pos) )
        end do
      end do 
    ! %/%
    else if (fun.eq.6) then
      do j = 1, n
        do i = 1, m
          call l2gpair(i, j, gi, gj, descx, blacs)
          pos = ind(gi + k*(gj-1), lvec)
          x(i, j) = fpmod( vec(pos), x(i, j) )
        end do
      end do
    ! <
    else if (fun.eq.7) then
      do j = 1, n
        do i = 1, m
          call l2gpair(i, j, gi, gj, descx, blacs)
          pos = ind(gi + k*(gj-1), lvec)
          if (x(i, j) .lt. vec(pos)) then
            x(i, j) = one
          else
            x(i, j) = zero
          end if
        end do
      end do
    ! >
    else if (fun.eq.8) then
      do j = 1, n
        do i = 1, m
          call l2gpair(i, j, gi, gj, descx, blacs)
          pos = ind(gi + k*(gj-1), lvec)
          if (x(i, j) .gt. vec(pos)) then
            x(i, j) = one
          else
            x(i, j) = zero
          end if
        end do
      end do
    ! <=
    else if (fun.eq.9) then
      do j = 1, n
        do i = 1, m
          call l2gpair(i, j, gi, gj, descx, blacs)
          pos = ind(gi + k*(gj-1), lvec)
          if (x(i, j) .le. vec(pos)) then
            x(i, j) = one
          else
            x(i, j) = zero
          end if
        end do
      end do
    ! >=
    else if (fun.eq.10) then
      do j = 1, n
        do i = 1, m
          call l2gpair(i, j, gi, gj, descx, blacs)
          pos = ind(gi + k*(gj-1), lvec)
          if (x(i, j) .ge. vec(pos)) then
            x(i, j) = one
          else
            x(i, j) = zero
          end if
        end do
      end do
    ! <
    else if (fun.eq.11) then
      do j = 1, n
        do i = 1, m
          call l2gpair(i, j, gi, gj, descx, blacs)
          pos = ind(gi + k*(gj-1), lvec)
          if (x(i, j) .eq. vec(pos)) then
            x(i, j) = one
          else
            x(i, j) = zero
          end if
        end do
      end do
    end if
  end if
  
  return
end subroutine



logical function checkproc(i, j, desc, blacs)
  integer i, j, desc(9), blacs(5)
  
  checkproc = ( mod( (i-1)/desc(5), blacs(2) ) .eq. blacs(4) .and. &
                mod( (j-1)/desc(6), blacs(3) ) .eq. blacs(5) )
  
  return
end function


! r-style replacement.  (dist)matrix-vector and (dist)matrix-(dist)matrix (levels
! 2 and 3) are available.
! inputs/outputs
  ! x = submatrix of data which should be replaced
! inputs
  ! ix/jx = 
  ! descx = descriptor array for x
  ! lvec = length of vec
subroutine rl2insert(x, ix, jx, descx, vec, lvec, indi, lindi, indj, lindj)
  implicit none
  ! in/out
  integer :: ix, jx, lvec, lindi, lindj
  integer :: indi(lindi), indj(lindj)
  integer :: descx(9)
  double precision :: x(descx(9), *), vec(lvec)
  ! local
  integer :: k, m, n, pos, i, j, ti, tj, gi, gj
  integer :: ldm(2), blacs(5)
  ! parameter
  double precision :: zero, one
  parameter ( zero = 0.0d0, one = 1.0d0 )
  ! external
  external            pdims, g2lpair
  ! function
  logical :: checkproc
  integer :: ind
  double precision :: fpmod
  
  
  ! get local and proc grid info
  call pdims(descx, ldm, blacs)
  
  m = ldm(1)
  n = ldm(2)
  k = descx(3)
  
  ! only do work if we own any local pieces
  if (m.gt.0 .and. n.gt.0) then
    ! insertion
    do tj = 1, lindj
      gj = indj(tj)
      do ti = 1, lindi
        gi = indi(ti)
        call g2lpair(i, j, gi, gj, descx, blacs)
        if (checkproc(gi, gj, descx, blacs)) then
          pos = ind(gi + k*(gj-1), lvec)
          x(i, j) = vec(pos)
        end if
      end do
    end do
  end if
  
  return
end subroutine



! column ycol of y is copied onto column xcol of x.  
! lcols = number of entries in xcols and ycols
subroutine rcolcpy(x, descx, xcols, y, descy, ycols, lcols)
  implicit none
  ! in/out
  integer :: lcols
  integer :: descx(9), descy(9)
  integer :: xcols(lcols), ycols(lcols)
  double precision :: x(descx(9), *), y(descy(9), *)
  ! local
  logical :: ihave, ineed
  integer :: i, j, gi, gj, mx, nx, my, ny, gm, gn
  integer :: rbl, cbl
  integer :: ldm(2), blacs(5)
  integer :: lxcol, lycol, xcol, ycol, col
  integer :: rsrc, csrc, rdest, cdest
  integer :: rxlen, rylen, rxtop, rytop
  ! external
  external            pdims, indxg2l, indxl2g
  ! function
  integer :: indxg2l, indxl2g
  
  
  gm = descx(3)
  gn = descx(4)
  
  if (gm.ne.descy(3)) return
  
  rbl = descx(5)
  cbl = descx(6)
  
  ! get local and proc grid info
  call pdims(descy, ldm, blacs)
  my = ldm(1)
  ny = ldm(2)
  
  call pdims(descx, ldm, blacs)
  
  mx = ldm(1)
  nx = ldm(2)
  
  do col = 1, lcols, 1
    xcol = xcols(col)
    ycol = ycols(col)
    
    lxcol = indxg2l(xcol, descx(6), 1, 1, blacs(3))
    lycol = indxg2l(ycol, descy(6), 1, 1, blacs(3))
    do gi = 1, gm, rbl
      ! index juggling
      i = indxg2l(gi, rbl, 1, 1, blacs(2))
      
      rxtop = min(i+rbl-1, mx)
      rxlen = rxtop - i + 1
      
      rytop = min(i+rbl-1, my)
      rylen = rytop - i + 1
      
      ! row and column (processor) source
      rsrc = mod( (gi-1)/rbl, blacs(2) )
      csrc = mod( (ycol-1)/cbl, blacs(3) )
      
      ihave = ( rsrc.eq.blacs(4) .and. csrc.eq.blacs(5) )
      
      ! row and column (processor) destination
      rdest = mod( (gi-1)/rbl, blacs(2) )
      cdest = mod( (xcol-1)/cbl, blacs(3) )
      
      ineed = ( rdest.eq.blacs(4) .and. cdest.eq.blacs(5) )
      
      ! copy
      if (ihave) then ! check if need to send
        if (ineed) then ! easy case
          x(i:rxtop, lxcol) = y(i:rytop, lycol)
        else ! otherwise send
          ! send
          call dgesd2d(descx(2), rylen, 1, y(i, lycol), rylen, rdest, cdest)
        end if
      else if (ineed) then ! otherwise check if need to receive
        ! receive
        call dgerv2d(descx(2), rxlen, 1, x(i, lxcol), rxlen, rsrc, csrc)
      end if
      
      call blacs_barrier(descx(2), 'A') 
      
    end do ! gi
  end do ! col
  
  return
end subroutine



! column ycol of y is copied onto column xcol of x. 
! lxcols is allowed to be an integer multiple of lycols
subroutine rcolcpy2(x, descx, xcols, lxcols, y, descy, ycols, lycols)
  implicit none
  ! in/out
  integer :: lxcols, lycols
  integer :: descx(9), descy(9)
  integer :: xcols(lxcols), ycols(lycols)
  double precision :: x(descx(9), *), y(descy(9), *)
  ! local
  integer :: col
  ! external
  external            rcolcpy
  
  
  if (mod(lxcols, lycols) /= 0) then
    return
  end if
  
  do col = 1, lxcols, lycols
    call rcolcpy(x, descx, xcols(col), y, descy, ycols, lycols)
  end do
  
  return
end subroutine



! rows yrows of y are copied onto rows xrows of x.  
! obvious assumptions are made.
subroutine rrowcpy(x, descx, xrows, y, descy, yrows, lrows)
  implicit none
  ! in/out
  integer :: lrows
  integer :: descx(9), descy(9)
  integer :: xrows(lrows), yrows(lrows)
  double precision :: x(descx(9), *), y(descy(9), *)
  ! local
  logical :: ihave, ineed
  integer :: i, j, gi, gj, mx, nx, my, ny, gm, gn
  integer :: rbl, cbl
  integer :: ldm(2), blacs(5)
  integer :: lxrow, lyrow, xrow, yrow, row
  integer :: rsrc, csrc, rdest, cdest
  integer :: cxlen, cylen, cxtop, cytop
  ! external
  external            pdims, indxg2l, indxl2g
  ! function
  integer :: indxg2l, indxl2g
  
  
  gm = descx(3)
  gn = descx(4)
  
  if (gn.ne.descy(4)) return
  
  rbl = descx(5)
  cbl = descx(6)
  
  ! get local and proc grid info
  call pdims(descy, ldm, blacs)
  my = ldm(1)
  ny = ldm(2)
  
  call pdims(descx, ldm, blacs)
  
  mx = ldm(1)
  nx = ldm(2)
  
  ! the work
  do row = 1, lrows, 1
    xrow = xrows(row)
    yrow = yrows(row)
    
    lxrow = indxg2l(xrow, descx(6), 1, 1, blacs(3))
    lyrow = indxg2l(yrow, descy(6), 1, 1, blacs(3))
    
!        j = indxg2l(gj, cbl, 1, 1, blacs(3))
    
    do gj = 1, gn, cbl
      ! index juggling
      j = indxg2l(gj, cbl, 1, 1, blacs(2))
      
      cxtop = min(j+cbl-1, nx)
      cxlen = cxtop - j + 1
      
      cytop = min(j+cbl-1, ny)
      cylen = cytop - j + 1
      
      ! row and column (processor) source
      rsrc = mod( (yrow-1)/rbl, blacs(2) )
      csrc = mod( (gj-1)/cbl, blacs(3) )
      
      ihave = ( rsrc.eq.blacs(4) .and. csrc.eq.blacs(5) )
      
      ! row and column (processor) destination
      rdest = mod( (xrow-1)/rbl, blacs(2) )
      cdest = mod( (gj-1)/cbl, blacs(3) )
      
      ineed = ( rdest.eq.blacs(4) .and. cdest.eq.blacs(5) )
      
      ! copy
      if (ihave) then ! check if need to send
        if (ineed) then ! easy case
          x(lxrow, j:cxtop) = y(lyrow, j:cytop)
        else ! otherwise send
          ! send
          call dgesd2d(descx(2), 1, cylen, y(lyrow, j), cylen, rdest, cdest)
        end if
      else if (ineed) then ! otherwise check if need to receive
        ! receive
        call dgerv2d(descx(2), 1, cxlen, x(lxrow, j), cxlen, rsrc, csrc)
      end if
      
      call blacs_barrier(descx(2), 'A') 
      
    end do ! gi
  end do ! row
  
  return
end subroutine



! rows yrows of y are copied onto rows xrows of x.  
! lxrows is allowed to be an integer multiple of lyrows
subroutine rrowcpy2(x, descx, xrows, lxrows, y, descy, yrows, lyrows)
  implicit none
  ! in/out
  integer lxrows, lyrows
  integer :: descx(9), descy(9)
  integer :: xrows(lxrows), yrows(lyrows)
  double precision :: x(descx(9), *), y(descy(9), *)
  ! local
  integer :: row
  ! external
  external            rrowcpy
  
  
  if (mod(lxrows, lyrows).ne.0) then
    return
  end if
  
  do row = 1, lxrows, lyrows
    call rrowcpy(x, descx, xrows(row), y, descy, yrows, lyrows)
  end do
  
  return
end subroutine




! distributed matrix-distributed vector sum (like r)
! x is an m_1 x n distributed matrix, y is a distributed m_2 x 1 vector
! x <-- x + y
subroutine pdmvsum(x, descx, y, descy)
  implicit none
  ! in/out
  integer :: descx(9), descy(9)
  double precision :: x(descx(9), *), y(descy(9))
  ! local
  integer :: ldm(2), blacs(5)
  integer :: i, j, gi, gj, ii, jj, pos
  integer :: gm, gn, rbl, cbl
  integer :: my, ny, mx, nx
  integer :: lvec
  double precision :: tmp
  
  
  ! get local and proc grid info
  gm = descx(3)
  gn = descx(4)
  
  lvec = descy(5)
  
  rbl = descx(5)
  cbl = descx(6)
  
  call pdims(descy, ldm, blacs)
  my = ldm(1)
  ny = ldm(2)
  
  call pdims(descx, ldm, blacs)
  mx = ldm(1)
  nx = ldm(2)
  
  ! quick return if possible
  if (descy(4) .ne. 1) return
  
  
  ! easy case:  nrows(x) == length(y)
  if (gm .eq. descy(3)) then
    do j = 1, nx
      do i = 1, mx
        x(i, j) = x(i, j) + y(i)
      end do
    end do
    
  ! general case
  !!! todo
  
  end if
  
  
  
  return
end subroutine

