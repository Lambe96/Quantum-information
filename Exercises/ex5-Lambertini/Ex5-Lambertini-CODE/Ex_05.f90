program ex_05
  use ex05module
  implicit none

  !!!!!!!!!!!!!!!!!!!! declaration !!!!!!!!!!!!!!!!!!!!
  complex(kind=8), allocatable :: matrix (:,:), work(:)
  integer :: msize, i, j, N, lwork, info, N_bins
  real(kind=8), allocatable :: rwork(:), eigv(:), spacing_vec(:), bins_centers(:),probs(:),matrix1(:,:)
  real(kind=8) :: start, finish, av_delta,r_mid
  logical :: norm

  !!!!!!!!!!!!!!!!!!!! initialization !!!!!!!!!!!!!!!!!!!!
  msize = 5000
  N = 500
  norm = .true.
  N_bins = 71

  !call rand_init_hermit_mat(matrix,msize)
  call rand_init_diag_mat(matrix1,msize)
  allocate(eigv(msize))
  do i =1,msize
    eigv(i)=matrix1(i,i)
  end do
  call dlasrt('I', msize, eigv, info)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !!!!!!!!!!!!!!!!!!!!! solve diagonalization problem !!!!!!!!!!!!!!!!!!!!
  !call diag_hermit_matrix(matrix,eigv)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!! compute <r> !!!!!!!!!!!!!!!!!!!!!
  call compute_r_mid(eigv,r_mid)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!! Spacing !!!!!!!!!!!!!!!!!!!!
  !call norm_spacing(eigv,spacing_vec)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!! moving average Spacing!!!!!!!!!!!!!!!!!!!!
  call mov_av_norm_spacing(eigv,spacing_vec,N)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  call PDF_vec(spacing_vec, N_bins, bins_centers, probs, norm)

  open(unit=4,file='diag_N_frac_10.csv',Access = 'append')
    do i=1,N_bins
      write(4,'(F0.15,A,F0.15)')  bins_centers(i),',', probs(i)
    end do
  close(4)

  open(unit=5,file='R_mid.csv',Access = 'append')
      write(5,'(A,I0,A,F0.15)')  'matrix_kind = d, N =', N ,'<r> = ', r_mid
  close(5)


  !!!!!!!!!!!!!!!!!!!!! print statements !!!!!!!!!!!!!!!!!!!!
  !write(*,*) 'Time:',finish-start
  !write(*,*) 'av_delta', av_delta
  !call write_a_vector('spacing_vec', spacing_vec)
  call write_a_vector('bins_centers', bins_centers)
  call write_a_vector('probs', probs)
  !call write_a_vector('eig', eigv)
  write(*,*) 'r_mid =',r_mid
  !call write_a_complex_matrix('Original matrix', matrix)
  !call write_a_vector('Eigenvalues', eigv)
  !call write_a_complex_matrix('Eigenvectors', matrix)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end program ex_05
