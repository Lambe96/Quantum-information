program q_harm_osc
  use ex06module
  implicit none

  complex(kind=8), dimension(:,:), allocatable :: matrix
  real(kind=8), dimension(:), allocatable :: xx
  real(kind=8) :: Re_z, x_max,step,freq,mass,h_bar
  integer(kind=8) ::  ii, N, k_columns, k_entries
  double precision, dimension(:), allocatable :: eigv

  !set the parameters of the problem
  N=1000
  x_max = 10
  freq = 1
  mass = 0.5
  h_bar = 1
  step = 2*x_max/real(N-1,8)
  k_columns = 3
  k_entries = 100
  !initialize the matrix and store it in a file
  call init_tridiag_qao_matrix(matrix, x_max,xx, N, freq, mass, h_bar)
  call vec_and_col_on_file('original_matrix.txt',real(matrix), N, xx)
  !call write_a_vector('x',xx)

  !diagonalize and normalize the matrix
  call diag_hermit_matrix(matrix,eigv)
  matrix = matrix/sqrt(step)


  !store the first k eigenfunctions and the first k eigenvalues in a file
  call vec_and_col_on_file('x_and_columns.txt',real(matrix), k_columns, xx)
  call first_k_entries_on_file('first_k_eigenvalues.txt', eigv, k_entries)
end program q_harm_osc
