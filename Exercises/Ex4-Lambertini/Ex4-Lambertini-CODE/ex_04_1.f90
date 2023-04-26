!Program useful to implement and test the matrix-matrix multiplication
!in Fortran. We implement two different loops to perform the multiplication
!and we compare the performance with the intrinsic Fortran function matmul().
!Moreover, we practice with the concept of debugging and error-handling through
!the use of debugging module.
!Author: Alessandro Lambertini --mat.1242885--


program matrix
  use debugging
  implicit none

  !define the variables
  integer, parameter :: dp = selected_real_kind(15, 307)
  real(kind=dp), dimension(:,:), allocatable :: m_1(:,:), m_2(:,:), m_fin(:,:),m_fin_2(:,:), m_int(:,:)
  real(kind=dp):: start, finish
  integer :: i, j, k
  integer, dimension(2) :: dim_m1, dim_m2
  logical :: debug
  !initialize matrices and variables
  debug = .true.

  open(unit=1, file='dim.txt')
  read(1,*) dim_m1(1),dim_m1(2),dim_m2(1),dim_m2(2)
  close(1)

  !dim_m1(1) = 4
  !dim_m1(2) = 4
  !dim_m2(1) = 4
  !dim_m2(2) = 3

  !Initialize the two matrices to be multiplied and check if the dimensions of
  !the two matrices are positive
  call rand_init_mat(m_1,dim_m1(1),dim_m1(2),debug)
  call checkpoints(debug, 'rand_init_mat', 23)

  call rand_init_mat(m_2,dim_m2(1),dim_m2(2),debug)
  call checkpoints(debug, 'rand_init_mat', 26)


  !check if the dimensions of the two matrices allows the multiplication
  !if so allocate the dimensions for the final matrices and print a friendly message
  !if not stop the program and print an error message

  call dim_check(dim_m1,dim_m2, debug)
  call checkpoints(debug, 'dim_check', 34)

  allocate(m_fin(dim_m1(1),dim_m2(2)),m_fin_2(dim_m1(1),dim_m2(2)), m_int(dim_m1(1),dim_m2(2)))


  !Write explicitely the matrix-matrix multiplication loop in two different orders.
  call cpu_time(start)

  call mat_mul_loop1(m_1,m_2,m_fin,debug)
  !Matrix multiplication with the first loop
  call checkpoints(debug, 'mat_mul_loop1',44)
  !check if the dimensions of the resulting matrix are the expected ones

  call cpu_time(finish)
  !write(*,*) 'Time loop 1', finish-start, 'seconds'
  open(unit=2,file='time_loop1.csv',Access = 'append')
  write(2,'(I0,A,F0.15)')  dim_m1(1),',', finish-start
  close(2)

  call cpu_time(start)

  call mat_mul_loop2(m_1,m_2,m_fin_2,debug)
  !Matrix multiplication with the second loop
  call checkpoints(debug, 'mat_mul_loop2', 54)
  !check if the dimensions of the resulting matrix are the expected ones

  call cpu_time(finish)
  !write(*,*) 'Time loop 2', finish-start, 'seconds'
  open(unit=3,file='time_loop2.csv',Access = 'append')
  write(3,'(I0,A,F0.15)')  dim_m1(1),',', finish-start
  close(3)

  !check that the dimensions are the ones expected after the multiplication
  !if not stop the program, print an error message and the expected dimensions vs the given ones

  !Use the Fortran intrinsic function.
  call cpu_time(start)
  m_int = matmul(m_1,m_2)
  call cpu_time(finish)
  !write(*,*) 'Intrinsic function', finish-start, 'seconds'
  open(unit=4,file='time_intrinsic.csv',Access = 'append')
  write(4,'(I0,A,F0.15)')  dim_m1(1),',', finish-start
  close(4)


  !print the results
  !write(*,*) 'm_1:'
  !call write_a_matrix(m_1)

  !write(*,*) 'm_2:'
  !call write_a_matrix(m_2)

  !write(*,*) 'm_fin:'
  !call write_a_matrix(m_fin)

  !write(*,*) 'm_fin_2:'
  !call write_a_matrix(m_fin_2)

  !write(*,*) 'm_int:'
  !call write_a_matrix(m_int)

  deallocate(m_1, m_2, m_fin, m_fin_2, m_int)

contains

  !subroutine that initialize the elements of a two-dimensional array with
  !real, double precision, random numbers between 1 and 100.
  !It takes as arguments the matrix name and the two dimensions
  !and an optional logical variable useful to interact with checkpoints
  !subroutine.
  subroutine rand_init_mat(mm,dim_1,dim_2,debug)
    logical, optional :: debug
    integer :: dim_1, dim_2, i, j
    real(kind=8), dimension(:,:), allocatable :: mm

    if (present(debug)) then
      if (dim_1 .ge. 1 .and. dim_2 .ge. 1) then
        debug = .TRUE.
      else
        debug = .FALSE.
      end if
    end if


    allocate(mm(dim_1,dim_2))

    do i=1,dim_2
      do j=1,dim_1
        call RANDOM_NUMBER(mm(j,i))
        mm(j,i) = mm(j,i)*100
      end do
    end do
  end subroutine rand_init_mat


  !subroutine that initialize the elements of a two-dimensional array with
  !real, double precision, random numbers between 1 and 100.
  !It takes as arguments the matrix name, the name of the file where are stored the
  !dimensions of the matrix, and optional logical variable useful to interact
  !with checkpoints subroutine.
  !subroutine rand_init_mat_file(mm,filename,debug)
  !  logical, optional :: debug
  !  integer :: dim_1, dim_2, i, j
  !  character(len=*) :: filename
  !  real(kind=8), dimension(:,:), allocatable :: mm

  !  open(unit=1, file=filename)
  !  read(1,*) dim_1,dim_2
  !  close(1)

  !  if (present(debug)) then
  !    if (dim_1 .ge. 1 .and. dim_2 .ge. 1) then
  !      debug = .TRUE.
  !    else
  !      debug = .FALSE.
  !    end if
  !  end if


  !  allocate(mm(dim_1,dim_2))

  !  do i=1,dim_2
  !    do j=1,dim_1
  !      call RANDOM_NUMBER(mm(j,i))
  !      mm(j,i) = mm(j,i)*100
  !    end do
  !  end do
  !end subroutine rand_init_mat_file



  !Subroutine that perform the mat-mat multiplication with the first loop
  !It takes as arguments the two matrices to be multiplied, the matrix
  !in which the result will be written, and an optional logical variable
  !useful to interact with checkpoints subroutine.
  subroutine mat_mul_loop1(mm_1,mm_2, m_fin, debug)
    logical, optional :: debug
    integer :: i,j,k
    real(kind=8), dimension(:,:) :: mm_1, mm_2, m_fin
    integer, dimension(2) :: dimm_1, dimm_2, dimm_fin

    dimm_1 = shape(mm_1)
    dimm_2 = shape(mm_2)
    dimm_fin = shape(m_fin)

    if (present(debug)) then
      if (dimm_fin(1) .eq. dimm_1(1) .and. dimm_fin(2) .eq. dimm_2(2)) then
        debug = .TRUE.
      else
        debug = .FALSE.
      end if
    end if

    do i=1,dimm_1(1)
      do j=1,dimm_2(2)
        do k=1,dimm_1(2)
          m_fin(i,j)=m_fin(i,j)+mm_1(i,k)*mm_2(k,j)
        end do
      end do
    end do
  end subroutine mat_mul_loop1

  !Subroutine that perform the mat-mat multiplication with the second loop
  !It takes as arguments the two matrices to be multiplied, the matrix
  !in which the result will be written, and an optional logical variable
  !useful to interact with checkpoints subroutine.
  subroutine mat_mul_loop2(mm_1,mm_2, m_fin, debug)
    logical, optional :: debug
    integer :: i,j,k
    real(kind=8), dimension(:,:) :: mm_1, mm_2, m_fin
    integer, dimension(2) :: dimm_1, dimm_2, dimm_fin

    dimm_1 = shape(mm_1)
    dimm_2 = shape(mm_2)
    dimm_fin = shape(m_fin)

    if (present(debug)) then
      if (dimm_fin(1) .eq. dimm_1(1) .and. dimm_fin(2) .eq. dimm_2(2)) then
        debug = .TRUE.
      else
        debug = .FALSE.
      end if
    end if

    do i=1,dimm_1(1)
      do k=1,dimm_1(2)
        do j=1,dimm_2(2)
          m_fin(i,j)=m_fin(i,j)+mm_1(i,k)*mm_2(k,j)
        end do
      end do
    end do
  end subroutine mat_mul_loop2


  !subroutine useful to write on screen a matrix in a readable way.
  subroutine write_a_matrix(mm)
    real(kind=8), dimension(:,:) :: mm
    integer :: i,j
    integer, dimension(2) :: dimm

    dimm = shape(mm)

    do i=1,dimm(1)
      do j=1,dimm(2)
          write(*,'(F0.15,3X)', advance='no') mm(i,j)
      end do
      write(*, *) ''
    end do
  end subroutine write_a_matrix


end program matrix
