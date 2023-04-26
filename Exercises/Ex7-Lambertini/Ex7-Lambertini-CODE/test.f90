program test
  use ex07module
  implicit none


  complex(kind=8), dimension(:,:), allocatable :: psi_evol
  complex(kind=8), dimension(100) :: psi
  integer(kind=8) :: t_resolution, ii, plan_tr, plan_atr
  real(kind=8) :: t_max, x_grid(10), t_grid(10), delta_t, omega, mass, h_bar, re, im

  t_max = 10
  t_resolution = 10

  delta_t = t_max/real(t_resolution-1,8)

  do ii=1,t_resolution
    t_grid(ii) = 0+(ii-1)*delta_t
    write(*,*) t_grid(ii)
  end do

  write(*,*) 'psi initialized'
  do ii=1,size(psi)
    call random_number(re)
    call random_number(im)
    psi(ii) = cmplx(100*re,100*im)
    write(*,*) psi(ii)
  end do

  call mod_squared_array(psi)

  do ii=1,size(psi)
    write(*,*) psi(ii)
  end do


  !call dfftw_plan_dft_1d(plan_tr,100,psi,psi,'FFTW_FORWARD','FFTW_ESTIMATE');
  !call dfftw_plan_dft_1d(plan_atr,100,psi,psi,'FFTW_BACKWARD','FFTW_ESTIMATE');

  !call dfftw_execute_dft(plan_tr,psi,psi)
  !write(*,*) 'psi trasformed:'
  !do ii=1,size(psi)
  !  write(*,*) psi(ii)
  !end do

  !call dfftw_execute_dft(plan_atr,psi,psi)
  !write(*,*) 'psi anti-trasformed:'
  !do ii=1,size(psi)
  !  write(*,*) psi(ii)
  !end do


end program test
