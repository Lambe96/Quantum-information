module matrices
  implicit none
  private                                        !define what is accesible
  public dmatrix, AB, TRA, ADJ, write_onf_dmat   !outside the module


  !derived type declaration
  type dmatrix
    integer, dimension(2) :: dim
    complex(kind=8), dimension(:,:), allocatable :: elements
    real(kind=8) :: trace, determinant
  end type dmatrix

  !Interfaces
  interface AB
    module procedure Initialize_dmatrix
  end interface

  interface TRA
    module procedure trace
  end interface

  interface ADJ
    module procedure Adjoint
  end interface


contains

  !function to initialize the new type. The commented parts are a first version initialization with random number.
  function Initialize_dmatrix(a,b)
    integer(kind=4) :: i,j
    integer(kind=4) :: a,b
    !complex :: ij
    !real(kind=8) :: Re_z,Im_z
    type(dmatrix) :: Initialize_dmatrix

    Initialize_dmatrix%dim(1) = a
    Initialize_dmatrix%dim(2) = b

    allocate(Initialize_dmatrix%elements(a,b))

    do i = 1, a
      do j = 1, b
        !CALL RANDOM_NUMBER(Re_z)
        !CALL RANDOM_NUMBER(Im_z)
        !ij = cmplx(Re_z,Im_z)
        !Initialize_dmatrix%elements(i,j) = ij
        Initialize_dmatrix%elements(i,j) = (0,0)
      end do
    end do

    Initialize_dmatrix%trace = 0
    Initialize_dmatrix%determinant = 0  !This will be always 0 untill we define
                                        !a funtion to calculate the determinant

    !if (abs(Initialize_dmatrix%dim(1) - Initialize_dmatrix%dim(2))>1e-1) then
    !  Initialize_dmatrix%trace = 999999999999999
    !  Initialize_dmatrix%determinant = 999999999999999
    !  write(*,*) 'This is not a square matrix so the trace and the determinant are not defined'
    !  write(*,*) 'The trace component is automatically set to', Initialize_dmatrix%trace
    !  write(*,*) 'The determinant component is automatically set to', Initialize_dmatrix%determinant
    !else
    !  do i = 1, Initialize_dmatrix%dim(1)
    !    Initialize_dmatrix%trace = Initialize_dmatrix%trace + Initialize_dmatrix%elements(i,i)
    !
    !    calculate the determinant

    !  end do

    !end if

  end function Initialize_dmatrix


  !function to compute the trace of a dmatrix.
  function trace(MM)
    type(dmatrix) :: MM
    integer(kind=4) :: i
    real(kind=8) :: trace

    trace = 0

    if (abs(MM%dim(1) - MM%dim(2))>1e-1) then
      ! write(*,*) 'This is not a square matrix'
      trace = 999999999999999
    else
      do i = 1, MM%dim(1)
        trace = trace + MM%elements(i,i)
      end do
    end if

  end function trace

  !Function to compute the conjugate transpose matrix of a dmatrix
  function Adjoint(MM)
    type(dmatrix) :: MM, Adjoint

    Adjoint%dim(1) = MM%dim(2)
    Adjoint%dim(2) = MM%dim(1)

    allocate(Adjoint%elements(Adjoint%dim(1),Adjoint%dim(2)))
    Adjoint%elements = conjg(transpose(MM%elements))

    adjoint%trace = MM%trace
    adjoint% determinant = 0  !still to be defined the function to compute the determinant

  end function Adjoint

  !subroutine useful for write on a file a dmatrix
  subroutine write_onf_dmat (MM,filename)
    type(dmatrix) :: MM
    character(len=*) :: filename
    integer :: i, j

    open (unit = 1, file = filename)
    Write(1,*) '# rows:', MM%dim(1)
    Write(1,*) '# columns:', MM%dim(2)
    write(1, *) ''

    write(1,*) 'Matrix:'
    write(1,*) ''
    do i = 1, MM%dim(1)
      do j = 1, MM%dim(2)
        write(1, '(ES0.0,"+i","(",ES0.0,")",3X)', advance='no') MM%elements(i,j)
      end do
      write(1, *) ''
    end do

    write(1,*) ''
    if (abs(MM%dim(1) - MM%dim(2))>1e-1) then
      write(1,*) 'Trace: undefined'
      write(1,*) 'Determinant: undefined'
    else
      write(1,*) 'Trace:', MM%trace
      write(1,*) 'Determinant:', MM%determinant
    end if
    close(1)
  end subroutine write_onf_dmat




    end module matrices

!compilation comments: compile this with
!gfortran -c matrix_module_try1.f90
