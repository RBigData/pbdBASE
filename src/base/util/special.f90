! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2013-2014, Schmidt


! construct hilbert matrices
subroutine dhilbmk(n, x)
  implicit none
  ! in/out
  integer             n
  double precision    x(n, n)
  ! local
  integer             i, j
  
  
  do j = 1, n
    do i = 1, n
      x(i, j) = 1.0d0/dble(i+j-1)
    end do
  end do
  
  return
end subroutine


subroutine pdhilbmk(x, descx)
  implicit none
  ! in/out
  integer             descx(9)
  double precision    x(descx(9), *)
  ! local
  integer             m, n, i, j, gi, gj, ldm(2), blacs(5)
  ! external
  external            pdims
  
  
  ! get local and proc grid info
  call pdims(descx, ldm, blacs)
  
  m = ldm(1)
  n = ldm(2)
  
  do j = 1, n
    do i = 1, m
      call l2gpair(i, j, gi, gj, descx, blacs)
      x(i,j) = 1.0d0/dble(gi+gj-1)
    end do
  end do
  
  return
end subroutine



! companion matrix constructor
! coef = (c_0, c_1, ..., c_{n-1}) where the c_i come from the polynomial:
! x^n + c_{n-1}x^{n-1} + ... + c_1 x + c_0
! coef is a global coefficients vector owned by all processors
subroutine pdmkcpn1(x, descx, coef)
  ! in/out
  integer             descx(9)
  double precision    x(descx(9), *), coef(*)
  ! local
  integer             i, j, m, n, gi, gj, ldm(2), blacs(5)
  ! parameters
  double precision    zero
  parameter( zero = 0.0d0 )
  
  
  ! get local and proc grid info
  call pdims(descx, ldm, blacs)
  
  m = ldm(1)
  n = ldm(2)
  
  ! fill first row
  do j = 1, n
    x(1, j) = zero
  end do
  
  call l2gpair(1, n, gi, gj, descx, blacs)
  if (gj.eq.descx(4)) then
    x(1, n) = -coef(1)
  end if
  
  ! fill identity submatrix
  do j = 1, n
    do i = 1, m
      call l2gpair(i, j, gi, gj, descx, blacs)
      if (gi.eq.gj+1) then
        x(i,j) = 1
      else 
        x(i,j) = zero
      end if
    end do
  end do
  
  ! fill last column
  do i = 1, m
    call l2gpair(i, n, gi, gj, descx, blacs)
    if (gj.eq.descx(4)) then
      x(i, n) = -coef(gi)
    end if
  end do
  
  
  return
end subroutine


