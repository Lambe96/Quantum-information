program ex10
  use ex10module
  implicit none

!declare the variable
complex(kind=8), dimension(:,:), allocatable :: identity, initial_system, &
H_a,h_b,H_l,H_r,H_ab,H_2N,H_2N_copy,P, H_l_prov,H_r_prov
real(kind=8), dimension(:), allocatable :: eigv,lambda_i,ground_states
real(kind=8) :: step,lambda_max
real(kind=8) :: norm_term
integer :: N, ii,jj,M,iter

!Create the lambda array
M = 100
lambda_max = 3
allocate(lambda_i(M),ground_states(M))

step = lambda_max/real(M-1,8)

do ii=1,M
  lambda_i(ii) = 0+(ii-1)*step
end do

!set parameters
N = 2
iter = 20
norm_term =  real(2,8)**iter

!initialize identities
allocate(identity(2**N,2**N))
do ii=1,2**N
  identity(ii,ii) = (1,0)
end do

!cycle over lambdas
do ii=1,M
!initialize the hamiltonian of the two subsystems
call H_ising_np_tr(N,lambda_i(ii),initial_system)

call sigma_x_i(N,N,H_l)
call sigma_x_i(N,1,H_r)
!Real Space Renormalization Group
do jj=1,iter

  call tensor_product(initial_system,identity,H_a)
  call tensor_product(identity,initial_system,H_b)
  call tensor_product(H_l,H_r,H_ab)

  H_2N = H_a + H_b + H_ab

  H_2N_copy = H_2N

  call diag_hermit_matrix(H_2N,eigv)

  P = H_2N(:,1:2**N)

  initial_system = MATMUL(TRANSPOSE(CONJG(P)),MATMUL(H_2N_copy,P))

  call tensor_product(H_l,identity,H_l_prov)
  call tensor_product(identity,H_r,H_r_prov)

  H_l = MATMUL( TRANSPOSE(CONJG(P)),MATMUL(H_l_prov,P))
  H_r = MATMUL( TRANSPOSE(CONJG(P)),MATMUL(H_r_prov,P))

  deallocate(H_a,H_b,H_ab,eigv,P,H_2N,H_l_prov,H_r_prov,H_2N_copy)

end do

ground_states(ii) = eigv(1)/norm_term/N

deallocate(initial_system)

end do

call write_arrays_on_file('ground_states.txt',lambda_i,ground_states)

end program ex10
