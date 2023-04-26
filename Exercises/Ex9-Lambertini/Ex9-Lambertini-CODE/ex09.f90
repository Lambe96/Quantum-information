program ex09
  use ex09module
  implicit none

complex(kind=8), dimension(:,:), allocatable :: sigma_z, sigma_x, identity, prov_prod, prod_id, ham
real(kind=8), dimension(:), allocatable :: eigv,lambda_i
real(kind=8), dimension(:,:), allocatable :: ground_states,four_eigv
real(kind=8) :: start,end,step,lambda_max
integer :: N, ii,jj,M,N_max
character(len=1024) :: filename

M = 100
N_max = 10
lambda_max = 3
allocate(lambda_i(M),ground_states(M,N_max))

step = lambda_max/real(M-1,8)

do ii=1,M
  lambda_i(ii) = 0+(ii-1)*step
end do


do N=2,N_max
  call cpu_time(start)
  allocate(four_eigv(M,4))
  do ii=1,M
    call H_ising_np_tr(N,lambda_i(ii),ham)
    call diag_hermit_matrix(ham,eigv)
    ground_states(ii,N-1) = eigv(1)/(N-1)
    four_eigv(ii,1) = eigv(1)/(N-1)
    four_eigv(ii,2) = eigv(2)/(N-1)
    four_eigv(ii,3) = eigv(3)/(N-1)
    four_eigv(ii,4) = eigv(4)/(N-1)
    deallocate(ham,eigv)
    write(*,*) 'N=',N,', M=',ii
  end do
  if ( N .lt. 10 ) then
    write (filename, "(A14,I1,A4)") "firstfoureigv_", N,'.txt'
  else
    write (filename, "(A14,I2,A4)") "firstfoureigv_", N,'.txt'
  end if
  call vec_and_col_on_file(filename,four_eigv,4,lambda_i)
  deallocate(four_eigv)
  call cpu_time(end)
  write(*,*) 'N=',N,' done, time needed:',end-start
end do


!call write_arrays_on_file('gs_N_3.txt',lambda_i,ground_states)
!call vec_and_col_on_file('gs_N.txt',ground_states,N_max-1,lambda_i)
!call write_a_complex_matrix('ham',ham)
!call write_a_vector('eigv',eigv)

end program ex09
