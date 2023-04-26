program qbits
  use ex08module
  implicit none

  type(pure_state) :: state
  integer(kind=8) :: dim, N
  logical :: separable
  complex(kind=8), dimension(:,:), allocatable :: rho, rho_a, rho_b
  real(kind=8) :: norm

  dim = 2
  N = 2
  separable = .true.
  norm = 2


  call init_pure_state(state,dim,N,separable)

  state%psi(1) = cmplx(1/sqrt(norm),0)
  state%psi(2) = 0
  state%psi(3) = 0
  state%psi(4) = cmplx(-1/sqrt(norm),0)

  call dyadics(state%psi,state%psi,rho)
  call partial_trace_a_he(rho, 2, 2, rho_b)
  call partial_trace_b_he(rho, 2, 2, rho_a)


  call write_a_complex_matrix('psi',rho)
  call write_a_complex_matrix('rho_reduced_over_a',rho_b)
  call write_a_complex_matrix('rho_reduced_over_b',rho_a)

end program qbits

!compile with 'gfortran -o ex08 ex08.f90 ex08_module.o -llapack -L/usr/local/lib -lfftw3'
