! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2013, Schmidt


! "optimal" process grid shape; tries to make the grid as close to square as 
! possible.  if a square grid is not possible, e.g. nprocs=10, then the 
! larger number goes to nrows, e.g. nrows=5, ncols=2.
! inputs
  ! nprocs = number of processors to be used.
! outputs
  ! nrows = "optimal" number of process rows
  ! ncols = "optimal" number of process columns
subroutine optimalgrid(nprocs, nrows, ncols)
  implicit none
  ! in/out
  integer             nprocs, nrows, ncols
  ! local
  integer             i, n
  
  
  n = int(sqrt(real(nprocs)))
  
  do i = 0, n-1
    ncols = n - i
    nrows = int(mod(nprocs, ncols))
    if (nrows==0) exit
  end do
  
  nrows = nprocs / ncols
  
  return
end subroutine




! reductions using blacs routines.  bindings meant to mimic allreduce, which is
! much more friendly than the weird blacs calls.
! inputs/outputs
  ! x = array of data to be reduced
! inputs
  ! descx = descriptor array of x
  ! op = operation: 'min', 'max', 'sum'.  default if you pass something else is 
    ! to sum.
  ! scope = blacs processor grid scope along which the reduction will take 
    ! place: 'all', 'column', 'row'
  ! rdest/cdest: row/column destination in the case of *reduce (rather than all)
subroutine dallreduce(x, descx, op, scope)
  implicit none
  ! in/out
  integer             descx(9)
  double precision    x( * )
  character           op, scope
  ! local
  integer             m, n, lda, ictxt
  ! external
  external            dgsum2d, dgamx2d, dgmn2d
  
  
  m = descx(3)
  n = descx(4)
  lda = descx(9)
  ictxt = descx(2)
  
  if (op == 'MIN') then
    call dgamn2d(ictxt, scope, ' ', m, n, x, lda, -1, -1, -1, -1, -1)
  else if (op == 'MAX') then
    call dgamx2d(ictxt, scope, ' ', m, n, x, lda, -1, -1, -1, -1, -1)
  else ! default to sum
    call dgsum2d(ictxt, scope, ' ', m, n, x, lda, -1, -1)
  end if
  
  return 
end subroutine


subroutine dreduce(x, descx, op, rdest, cdest, scope)
  implicit none
  ! in/out
  integer             descx(9), rdest, cdest
  double precision    x( * )
  character           op, scope
  ! local
  integer             m, n, lda, ictxt
  ! external
  external            dgsum2d, dgamx2d, dgmn2d
  
  
  m = descx(3)
  n = descx(4)
  lda = descx(9)
  ictxt = descx(2)
  
  if (op == 'MIN') then
    call dgamn2d(ictxt, scope, ' ', m, n, x, lda, -1, -1, -1, rdest, cdest)
  else if (op == 'MAX') then
    call dgamx2d(ictxt, scope, ' ', m, n, x, lda, -1, -1, -1, rdest, cdest)
  else
    call dgsum2d(ictxt, scope, ' ', m, n, x, lda, rdest, cdest)
  end if
  
  return 
end subroutine


subroutine iallreduce(x, descx, op, scope)
  implicit none
  ! in/out
  integer             descx(9), x( * )
  character           op, scope
  ! local
  integer             m, n, lda, ictxt
  ! external
  external            igsum2d, igamx2d, igmn2d
  
  
  m = descx(3)
  n = descx(4)
  lda = descx(9)
  ictxt = descx(2)
  
  if (op == 'MIN') then
    call igamn2d(ictxt, scope, ' ', m, n, x, 1, -1, -1, -1, -1, -1)
  else if (op == 'MAX') then
    call igamx2d(ictxt, scope, ' ', m, n, x, 1, -1, -1, -1, -1, -1)
  else
    call igsum2d(ictxt, scope, ' ', m, n, x, 1, -1, -1)
  end if
  
  return 
end subroutine


subroutine ireduce(x, descx, op, rdest, cdest, scope)
  implicit none
  ! in/out
  integer             descx(9), x( * ), rdest, cdest
  character           op, scope
  ! local
  integer             m, n, lda, ictxt
  ! external
  external            igsum2d, igamx2d, igmn2d
  
  
  m = descx(3)
  n = descx(4)
  lda = descx(9)
  ictxt = descx(2)
  
  if (op == 'MIN') then
    call igamn2d(ictxt, scope, ' ', m, n, x, 1, -1, -1, -1, rdest, cdest)
  else if (op == 'MAX') then
    call igamx2d(ictxt, scope, ' ', m, n, x, 1, -1, -1, -1, rdest, cdest)
  else
    call igsum2d(ictxt, scope, ' ', m, n, x, 1, rdest, cdest)
  end if
  
  return 
end subroutine

