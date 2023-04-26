program test
  use ex08module
  implicit none

complex(kind=8), dimension(:), allocatable :: x,y
complex(kind=8), dimension(:,:), allocatable :: mat
real(kind=8) :: re, im
integer :: N, ii

N=4
allocate(x(N),y(N))

do ii=1,N
  call random_number(re)
  call random_number(im)
  x(ii) = cmplx(re,im)
  write(*,*) x(ii)
end do
write(*,*)'--------------------------------'
do ii=1,N
  call random_number(re)
  call random_number(im)
  y(ii) = cmplx(re,im)
  write(*,*) y(ii)
end do


call dyadics(x,y,mat)
call write_a_complex_matrix('mat',mat)

end program test
