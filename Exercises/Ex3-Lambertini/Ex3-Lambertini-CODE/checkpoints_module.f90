!This module implement some simple check on matrix-like object.
!The first subroutine implement a logical control on the Execution
!of other routines and it is useful to check their proper execution and so
!keep track of errors and bug.
!The second and third subroutines are implementing the same dimensionality check
!on the objects to verify the possibility to multiply them.
!Author: Alessandro Lambertini --mat.1242885--



module debugging
  implicit none
  private
  public dim_check, checkpoints

  interface dim_check
    module procedure dimension_mul_check
    module procedure dimensionality_check
  end interface




contains

  !Subroutine useful for debugging that take a logical variable from another subroutines
  !and decide whether to continue the execution or print some info and stop the program.
  subroutine checkpoints (debug, sub_func_name, line)
    logical :: debug
    character(len=*) :: sub_func_name
    integer*4 :: line
    if (debug .eqv. .true.) then
      write(*,"(A,X,L,X,A,X,A,A,X,A,I0.0)") 'debug =', debug, &
      'for the func/subrt:', sub_func_name,'.','Line:',line
    else
      write(*,"(A,X,L,X,A,X,A,A,X,A,I0.0)") 'debug =', debug, &
      'for the func/subrt:', sub_func_name,'.','Line:',line
      stop
    end if

  end subroutine checkpoints

  !Subroutine to check if two matrices can be multiplied.
  !It takes as arguments the two matrices and an optional logical variable
  !useful to interact with checkpoints subroutine.
  subroutine dimension_mul_check(mm_1,mm_2, debug)
    logical, optional :: debug
    real, dimension(:,:) :: mm_1, mm_2
    integer,dimension(2) :: col_1, row_2
    col_1 = shape(mm_1)
    row_2 = shape(mm_2)
    if (col_1(2) .ne. row_2(1)) then
      write(*,*) "This two matrices can't be multiplied, the program will stop"
      if (present(debug)) then
        debug = .false.
      else
        stop
      end if
    else
      write(*,*) "This two matrices can be multiplied"
      if (present(debug)) then
        debug = .true.
      end if
    end if
  end subroutine dimension_mul_check

  !Subroutine to check if two matrices can be multiplied.
  !It takes as arguments two one-dim arrays with two entries,
  !the dimensions of the two matrices,
  !and an optional logical variable
  !useful to interact with checkpoints subroutine.
  subroutine dimensionality_check(dim1,dim2, debug)
    logical, optional :: debug
    integer,dimension(2) :: dim1, dim2

    debug = .TRUE.

    if (dim1(2) .ne. dim2(1)) then
      write(*,*) "This two matrices can't be multiplied, the program will stop!"
      if (present(debug)) then
        debug = .false.
      else
        stop
      end if
    else
      write(*,*) "This two matrices can be multiplied!"
      if (present(debug)) then
        debug = .true.
      end if
    end if
  end subroutine dimensionality_check


end module debugging
