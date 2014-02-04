! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2013-2014, Schmidt


! sweep array out of distributed matrix
! inputs/outputs
  ! x = submatrix of data which should globally be "swept"
! inputs
  ! ix/jx = 
  ! descx = descriptor array for x
  ! vec = vector to "sweep" through x
  ! lvec = length of vec
  ! margin = 1 for row sweeping, 2 for column sweeping
  ! fun = char with 4 possibilities, describing the type of sweep to perform:
    ! "+", "-", "*", "/"
subroutine pdsweep(x, ix, jx, descx, vec, lvec, margin, fun)
  implicit none
  ! in/out
  integer             ix, jx, descx(9), margin, lvec
  double precision    x(descx(9), *), vec(lvec)
  character*1         fun
  ! local
  integer             k, m, n, pos, i, j, gi, gj, ldm(2), blacs(5)
  ! external
  external            pdims
  ! function
  integer             ind
  
  
  ! get local and proc grid info
  call pdims(descx, ldm, blacs)
  
  m = ldm(1)
  n = ldm(2)
  
  ! only do work if we own any local pieces
  if (m.gt.0 .and. n.gt.0) then
    ! addition
    if (fun.eq."+") then
      if (margin.eq.1) then
        k = descx(3)
        do j = 1, n
          do i = 1, m
            call l2gpair(i, j, gi, gj, descx, blacs)
            pos = ind(gi + k*(gj-1), lvec)
            x(i, j) = x(i, j) + vec(pos)
          end do
        end do
      else if (margin.eq.2) then
        k = descx(4)
        do j = 1, n
          do i = 1, m
            call l2gpair(i, j, gi, gj, descx, blacs)
            pos = ind(gj + k*(gi-1), lvec)
            x(i, j) = x(i, j) + vec(pos)
          end do
        end do
      end if
    ! subtraction
    else if (fun.eq."-") then
      if (margin.eq.1) then
        k = descx(3)
        do j = 1, n
          do i = 1, m
            call l2gpair(i, j, gi, gj, descx, blacs)
            pos = ind(gi + k*(gj-1), lvec)
            x(i, j) = x(i, j) - vec(pos)
          end do
        end do
      else if (margin.eq.2) then
        k = descx(4)
        do j = 1, n
          do i = 1, m
            call l2gpair(i, j, gi, gj, descx, blacs)
            pos = ind(gj + k*(gi-1), lvec)
            x(i, j) = x(i, j) - vec(pos)
          end do
        end do
      end if
    ! multiplication
    else if (fun.eq."*") then
      if (margin.eq.1) then
        k = descx(3)
        do j = 1, n
          do i = 1, m
            call l2gpair(i, j, gi, gj, descx, blacs)
            pos = ind(gi + k*(gj-1), lvec)
            x(i, j) = x(i, j) * vec(pos)
          end do
        end do
      else if (margin.eq.2) then
        k = descx(4)
        do j = 1, n
          do i = 1, m
            call l2gpair(i, j, gi, gj, descx, blacs)
            pos = ind(gj + k*(gi-1), lvec)
            x(i, j) = x(i, j) * vec(pos)
          end do
        end do
      end if
    ! division
    else if (fun.eq."/") then
      if (margin.eq.1) then
        k = descx(3)
        do j = 1, n
          do i = 1, m
            call l2gpair(i, j, gi, gj, descx, blacs)
            pos = ind(gi + k*(gj-1), lvec)
            x(i, j) = x(i, j) / vec(pos)
          end do
        end do
      else if (margin.eq.2) then
        k = descx(4)
        do j = 1, n
          do i = 1, m
            call l2gpair(i, j, gi, gj, descx, blacs)
            pos = ind(gj + k*(gi-1), lvec)
            x(i, j) = x(i, j) / vec(pos)
          end do
        end do
      end if
    end if
  end if
  
  return
end subroutine

