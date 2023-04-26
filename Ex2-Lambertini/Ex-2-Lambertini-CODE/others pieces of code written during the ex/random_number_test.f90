program test_rand

  implicit none

  integer,parameter :: seed = 86456
  integer(kind=4) :: a,b
  real:: Re_z, Im_z
  complex :: ij

  complex, dimension(:,:), allocatable :: elements

  a = 2
  b = 2
  allocate(elements(a,b))

  do i = 1, a
    do j = 1, b
      CALL RANDOM_NUMBER(Re_z)
      CALL RANDOM_NUMBER(Im_z)
      write(*,*) Re_z,Im_z
      ij = cmplx(Re_z,Im_z)
      elements(j,i) = ij
    end do
  end do
  write(*,*) elements
end program test_rand
