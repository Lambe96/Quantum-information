program matrix_mod_test
  use matrices
  implicit none

  type(dmatrix) :: mat_1, mat_1_t
  integer :: i,j
  complex(kind=8) :: ij
  real(kind=8) :: Re_z, Im_z, tr

  mat_1 = AB(4,4)   !initialize a 4X4 dmatrix

  do i = 1, mat_1%dim(1)
    do j = 1, mat_1%dim(2)
      CALL RANDOM_NUMBER(Re_z)
      CALL RANDOM_NUMBER(Im_z)      !fill the matrix with random numbers
      ij = cmplx(Re_z,Im_z)
      mat_1%elements(i,j) = ij
    end do
  end do



  mat_1%trace = TRA(mat_1) !calculate the trace
  mat_1%determinant = 0    !still not defined a function to calculate the determinant

  mat_1_t = ADJ(mat_1)    !calculate the transpose conjugate matrix

  call write_onf_dmat(mat_1,'mat_1.txt')    !write mat_1 into 'mat_1.txt'
  call write_onf_dmat(mat_1_t,'mat_1_t.txt')

  deallocate(mat_1%elements)    !deallocate memory
  deallocate(mat_1_t%elements)

end program matrix_mod_test

!compilation comments: compile this with
!gfortran -o [name executable] matrix_mod_test.f90 matrix_module_try1.o
