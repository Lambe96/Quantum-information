program test
  implicit none

  type dmatrix
    integer, dimension(2) :: dim
    complex(kind=8), dimension(:,:), allocatable :: elements
  end type dmatrix

    function Initialize_dmatrix(a,b)
      integer(kind=8) :: i,j
      integer(kind=4) :: a,b
      complex :: ij
      real(kind=8) :: Re_z,Im_z
      type(dmatrix) :: Initialize_dmatrix

      Initialize_dmatrix%dim(1) = a
      Initialize_dmatrix%dim(2) = b

      allocate(Initialize_dmatrix%elements(a,b))

      do i = 1, a
        do j = 1, b
          CALL RANDOM_NUMBER(Re_z)
          CALL RANDOM_NUMBER(Im_z)
          ij = cmplx(Re_z,Im_z)
          Initialize_dmatrix%elements(i,j) = ij
        end do
      end do

    end function Initialize_dmatrix
  end program test
