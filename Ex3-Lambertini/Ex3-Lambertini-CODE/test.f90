!program test
!  implicit none
!
!  logical :: debug
!  real, dimension(2,2) :: mm_1
!  real, dimension(3,2) :: mm_2
!  integer,dimension(2) :: col_1, row_2
!
!
!  mm_1(1,1) = 1
!  mm_1(1,2) = 1
!  mm_1(2,1) = 1
!  mm_1(2,2) = 1
!
!  mm_2(1,1) = 1
!  mm_2(1,2) = 1
!  mm_2(2,1) = 1
!  mm_2(2,2) = 1
!  mm_2(3,1) = 1
!  mm_2(3,2) = 1
!
!
!  debug = .TRUE.
!  col_1 = shape(mm_1)
!  row_2 = shape(mm_2)
!
!  if (col_1(1) .ne. row_2(1)) then
!    write(*,*) "This two matrices can't be multiplied, the program will stop"
!    debug = .FALSE.
!    stop
!  else
!    write(*,*) "This two matrices can be multiplied"
!  end if
!
!end program test

program test_kind
  implicit none
  real(kind=4) :: m=345.45532
  complex :: a=(34,67)
  integer,parameter :: kc = kind('s')
  integer,parameter :: kl = kind(.true.)
  integer,parameter :: st = kind(m)
  integer,parameter :: cp = kind(a)
  integer :: cpp
  real(kind=8), dimension(:,:), allocatable :: m_1(:,:), m_2(:,:), m_fin(:,:),m_fin_2(:,:), m_int(:,:)
  real(kind=8):: start, finish
  integer :: i, j, k
  integer, dimension(2) :: dim_m1, dim_m2,dimm

  dim_m1(1) = 3
  dim_m1(2) = 5

  print *, "The default character kind is ", kc
  print *, "The default logical kind is ", kl
  print *, "The default double kind is ", st
  print *, "The default complex kind is ", cp

  call rand_init_mat(m_1, dim_m1(1), dim_m1(2))
  write(*,*) m_1
  dimm = shape(m_1)
  write(*,*) kind(dimm(1))

contains
  subroutine rand_init_mat(mm,dim_1,dim_2)
    integer :: dim_1, dim_2, i, j
    real(kind=8), dimension(:,:), allocatable :: mm

    allocate(mm(dim_1,dim_2))

    do i=1,dim_2
      do j=1,dim_1
        call RANDOM_NUMBER(mm(j,i))
        mm(j,i) = mm(j,i)*100
      end do
    end do
  end subroutine rand_init_mat

end program test_kind
