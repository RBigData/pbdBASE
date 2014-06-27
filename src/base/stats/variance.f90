! this source code form is subject to the terms of the mozilla public
! license, v. 2.0. if a copy of the mpl was not distributed with this
! file, you can obtain one at http://mozilla.org/mpl/2.0/.

! copyright 2013, schmidt


! colmeans
subroutine pdclmn(x, descx, mn)
  implicit none
  ! in/out
  integer             descx( 9 )
  double precision    x(descx(9), *), mn(*)
  ! local
  integer             m, n, i, j, ldm(2), blacs(5), nn
  double precision    denom
  !external
  external            pdims, dgsum2d
  
  
  ! get local and proc grid info
  call pdims(descx, ldm, blacs)
  m = ldm(1)
  n = ldm(2)
  
  denom = dble(descx(3))
  
  do j = 1, n
    mn(j) = 0.0d0
    do i = 1, m
      mn(j) = mn(j) + x(i,j)
    end do
    mn(j) = mn(j) / denom
  end do
  
  nn = n
  call igamx2d(descx(2), 'col', ' ', 1, 1, nn,1,-1,-1,-1,-1,-1)
  
  call dgsum2d(descx(2), 'col', ' ', nn, 1, mn, nn, -1, -1)
  
  return
end subroutine


! colvar
subroutine pdclvar(x, descx, var)
  implicit none
  ! in/out
  integer             descx( 9 )
  double precision    x(descx(9), *), var(*)
  ! local
  integer             m, n, i, j, ldm(2), blacs(5), allocerr, nn
  double precision    scl, denom
  double precision, allocatable :: mn(:), work(:,:)
  ! parameter
  double precision    zero
  parameter ( zero = 0.0d0 )
  !external
  external            pdims, pdclmn, igamx2d, dgsum2d
  double precision    pdvar
  
  
  ! get local and proc grid info
  call pdims(descx, ldm, blacs)
  m = ldm(1)
  n = ldm(2)
  
  denom = dble(descx(3))
  
  nn = n
  call igamx2d(descx(2), 'col', ' ', 1, 1, nn,1,-1,-1,-1,-1,-1)
  
  allocerr = 0
  allocate(mn(nn), stat=allocerr)
  if (allocerr /= 0) return ! stop "out of memory"
  allocate(work(2, nn), stat=allocerr)
  if (allocerr /= 0) return ! stop "out of memory"
  
  mn(1:n) = zero
  call pdclmn(x, descx, mn)
  
  work(1:2, 1:nn) = zero
  
  do j = 1, n
    do i = 1, m
      work(1, j) = work(1, j) + (x(i,j) - mn(j))**2
      work(2, j) = work(2, j) + x(i,j) - mn(j)
    end do
  end do
  
  call dgsum2d(descx(2), 'col', ' ', 2*nn, 1, work, 2*nn, -1, -1)
  
  do j = 1, n
    var(j) = (work(1, j) - (work(2, j)**2)/denom) / (denom - 1)
  end do
  
  deallocate(mn)
  deallocate(work)
  
  return
end subroutine


