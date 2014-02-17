!module asdfmod
!  contains
!  
!  subroutine asdf(x)
!    integer, allocatable :: x(:,:)
!    
!    allocate(x(4,2))
!    
!    x = 1
!    
!    return
!  end subroutine
!  
!end module

!program main
!  use asdfmod
!  integer, allocatable :: x(:,:)
!  call asdf(x)
!  
!  print *, x
!  print *, ubound(x, dim=1), ubound(x, dim=2)
!end program


program test
  use subset
  integer :: my, ny
  integer :: i, j
  integer, allocatable :: rows(:)
  integer :: lrows
  integer, allocatable :: cols(:)
  integer :: lcols
  
  double precision, allocatable :: x(:,:)
  double precision , allocatable :: y(:,:)
  
  m = 3
  n = 3
  
  allocate(x(m, n))
  
  do j = 1, n
    do i = 1, m
      x(i, j) = i + m*(j-1)
    end do
  end do
  
  
  allocate(rows(1))
  lrows = -1
  
  allocate(cols(2))
  cols = (/ 1, 2 /)
  lcols = 2
  call subsetter(m, n, x, rows, lrows, cols, lcols, y, my, ny, info)
  print *, my, ny
  print *, y
  deallocate(cols)
  
  cols = (/ -3 /)
  lcols = 1
  call subsetter(m, n, x, rows, lrows, cols, lcols, y, my, ny, info)
  print *, my, ny
  print *, y
  deallocate(cols)
  
  
  deallocate(rows)
  allocate(cols(1))
  lcols = -1
  
  allocate(rows(2))
  rows = (/ 1, 2 /)
  lrows = 2
  call subsetter(m, n, x, rows, lrows, cols, lcols, y, my, ny, info)
  print *, my, ny
  print *, y
  deallocate(rows)
  
  rows = (/ -3 /)
  lrows = 1
  call subsetter(m, n, x, rows, lrows, cols, lcols, y, my, ny, info)
  print *, my, ny
  print *, y
  
end program



