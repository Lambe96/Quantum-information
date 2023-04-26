program test
  implicit none

  !complex(kind=16), dimension(:,:), allocatable :: matrix
  real(kind=8), dimension(:), allocatable :: xx
  real(kind=8) :: Re_z, Im_z,x_max,step,freq,mass,h_bar
  integer(kind=8) ::  ii, jj, kk, N

  N=25
  allocate(xx(N))
  x_max = 7
  step = 2*x_max/real(N-1,8)

  do ii=1,N
    xx(ii) = -x_max+(ii-1)*step
  end do

  do ii=1,N
    write(*,*) xx(ii)
  end do



end program test
