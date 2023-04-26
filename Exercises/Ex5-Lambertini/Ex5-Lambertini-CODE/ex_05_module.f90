module ex05module
  implicit none

contains

  !!!!!!!!!!!!!!!!!!! Initialize random hermitian matrix !!!!!!!!!!!!!!!!!!
  subroutine rand_init_hermit_mat(mm,dim_1,debug)
    logical, optional :: debug
    integer(kind=4) :: dim_1, i, j, k
    real(kind=8) :: Re_z, Im_z
    complex(kind=8) :: ij
    complex(kind=8), dimension(:,:), allocatable :: mm

    k=0

    if (present(debug)) then
      if (dim_1 .ge. 1) then
        debug = .TRUE.
      else
        debug = .FALSE.
      end if
    end if


    allocate(mm(dim_1,dim_1))

    do i=1,dim_1
      k = k+1
      do j=1,k
        if (i .ne. j) then
          call RANDOM_NUMBER(Re_z)
          Re_z = 2000*Re_z-1000
          call RANDOM_NUMBER(Im_z)
          Im_z = 2000*Im_z-1000
          mm(i,j) = cmplx(Re_z,Im_z)
          !mm(j,i) = conjg(mm(i,j))
        else if (i .eq. j) then
          call RANDOM_NUMBER(Re_z)
          Re_z = 2000*Re_z-1000
          Im_z = 0
          mm(i,j) = cmplx(Re_z,Im_z)
        end if
      end do
    end do
  end subroutine rand_init_hermit_mat
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!! Initialize random diagonal real matrix !!!!!!!!!!!!!!!!!!
  subroutine rand_init_diag_mat(mm,dim_1,debug)
    logical, optional :: debug
    integer(kind=4) :: dim_1, i
    real(kind=8), dimension(:,:), allocatable :: mm

    if (present(debug)) then
      if (dim_1 .ge. 1) then
        debug = .TRUE.
      else
        debug = .FALSE.
      end if
    end if


    allocate(mm(dim_1,dim_1))

    do i=1,dim_1
      call RANDOM_NUMBER(mm(i,i))
      mm(i,i) = mm(i,i)*2000-1000
    end do
  end subroutine rand_init_diag_mat
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!! solve diagonalization problem !!!!!!!!!!!!!!!!!!!!
  subroutine diag_hermit_matrix(matrix,eigv)
    complex(kind=8), allocatable :: work(:)
    complex(kind=8), dimension(:,:) :: matrix
    integer, dimension(2) :: dimm
    integer :: msize, lwork, info
    double precision, allocatable :: rwork(:), eigv(:)

    dimm = shape(matrix)
    msize = dimm(1)

    allocate(eigv(msize))
    allocate(work(msize*msize))
    allocate(rwork(3*msize-2))
    lwork = msize*msize

    call zheev('N','L',msize,matrix,msize,eigv,work,-1,rwork,info)
    lwork = work(1)

    call zheev('V','L',msize,matrix,msize,eigv,work,lwork,rwork,info)


    if (info.gt.0) then
      write(*,*)'The algorithm failed to compute eigenvalues.'
    end if
  end subroutine diag_hermit_matrix
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!! Spacing !!!!!!!!!!!!!!!!!!!!
  subroutine norm_spacing(eigv,spacing_vec)
    real(kind=8), allocatable :: spacing_vec(:)
    real(kind=8), dimension(:) :: eigv
    real(kind=8) :: av_delta
    integer :: i

    allocate(spacing_vec(size(eigv)-1))
    av_delta = 0
    do i=1,size(spacing_vec)
      spacing_vec(i) = eigv(i+1)-eigv(i)
      if (abs(spacing_vec(i)) < 1e-15) then
        spacing_vec(i) = 0
      end if
    end do
    av_delta = sum(spacing_vec)/size(spacing_vec)
    spacing_vec = spacing_vec/av_delta
  end subroutine norm_spacing
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!! moving average Spacing!!!!!!!!!!!!!!!!!!!!
  subroutine mov_av_norm_spacing(eigv,spacing_vec,N)
    real(kind=8), allocatable :: spacing_vec(:)
    real(kind=8), dimension(:) :: eigv
    real(kind=8) :: av_delta
    integer :: i, j, N

    allocate(spacing_vec(size(eigv)-1))
    av_delta = 0

    do i=1,size(spacing_vec)
      spacing_vec(i) = eigv(i+1)-eigv(i)
      if (abs(spacing_vec(i)) < 1e-15) then
        spacing_vec(i) = 0
      end if
    end do

    do i=1,size(spacing_vec)
      av_delta = 0
      if ((i > N/2) .and. (i<size(spacing_vec)-N/2)) then                    !be careful that int/int=int truncated to 0, then 1999/1000=1  (it should be i<size(spacing_vec)-(n/2+1))
        do j=(i-N/2),(i+N/2)
          av_delta=av_delta+spacing_vec(j)
        end do

        av_delta = av_delta/N

        spacing_vec(i) = spacing_vec(i)/av_delta


      else if (i.le.N/2) then
        av_delta = sum(spacing_vec(1:N))/N
        spacing_vec(i) = spacing_vec(i)/av_delta
      else if(i.ge.size(spacing_vec)-1-N/2) then
        av_delta = sum(spacing_vec(size(spacing_vec)-N:size(spacing_vec)))/N
        spacing_vec(i) = spacing_vec(i)/av_delta
      end if
    end do
  end subroutine mov_av_norm_spacing
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !!!!!!!!!!!!!!!!!!!!! compute the distribution of the elements of a real vector !!!!!!!!!!!!!!!!!!!!
  subroutine PDF_vec(spacing_vec, N_bins, bins_centers, probs, norm)
    real(kind=8), dimension(:) :: spacing_vec
    real(kind=8), dimension(:), allocatable :: probs, bins_centers, counts,bins_edges
    integer :: N_bins, i, j
    real(kind=8) :: area, ds
    logical :: norm

    allocate(probs(N_bins), bins_centers(N_bins), counts(N_bins),bins_edges(N_bins+1))

    ds = (maxval(spacing_vec)-minval(spacing_vec))/N_bins

    do i=1,N_bins
      bins_centers(i) = minval(spacing_vec) + (i - 0.5)*ds
    end do

    do i=1,N_bins+1
      bins_edges(i) = minval(spacing_vec) + (i-1)*ds
    end do

    counts = 0

    do i=1,size(spacing_vec)
      do j=1,N_bins
        if ((spacing_vec(i).ge.bins_edges(j)) .and. (spacing_vec(i).le.bins_edges(j+1))) then
          counts(j)=counts(j)+1
        end if
      end do
    end do

    if (norm .eqv. .true.) then
      area= sum(counts*ds)
      probs = counts/area
    else if (norm .eqv. .false.) then
      probs = counts
    end if
  end subroutine PDF_vec
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!! Compute <r> !!!!!!!!!!!!!!!!!!!!
  subroutine compute_r_mid(eigv,r_mid)
    real(kind=8), allocatable :: spacing_vec(:), rr(:)
    real(kind=8), dimension(:) :: eigv
    real(kind=8) :: r_mid
    integer :: i

    allocate(spacing_vec(size(eigv)-1), rr(size(eigv)-2))

    do i=1,size(spacing_vec)
      spacing_vec(i) = eigv(i+1)-eigv(i)
      if (abs(spacing_vec(i)) < 1e-15) then
        spacing_vec(i) = 0
      end if
    end do

    do i=1,size(rr)
      rr(i) = min(spacing_vec(i),spacing_vec(i+1))/max(spacing_vec(i),spacing_vec(i+1))
    end do

    r_mid = sum(rr)/size(rr)
  end subroutine compute_r_mid
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  !!!!!!!!!!!!!!!!!!!!! write a complex matrix !!!!!!!!!!!!!!!!!!!!
  subroutine write_a_complex_matrix(title,mm)
    complex(kind=8), dimension(:,:) :: mm
    integer(kind=8) :: i,j
    integer(kind=8), dimension(2) :: dimm
    character(len=*) :: title

    dimm = shape(mm)

    write(*,*) title

    do i=1,dimm(1)
      do j=1,dimm(2)
          write(*,'(ES0.5,"+i","(",ES0.5,")",3X)', advance='no') mm(i,j)
      end do
      write(*,*) ''
    end do
  end subroutine write_a_complex_matrix
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!! write a real vector !!!!!!!!!!!!!!!!!!!!
  subroutine write_a_vector(title,mm)
    real(kind=8), dimension(:) :: mm
    integer(kind=8) :: i
    integer(kind=8) :: dimm
    character(len=*) :: title

    dimm = size(mm)

    write(*,*) title

    do i=1,dimm
      write(*,'(ES0.5)') mm(i)
    end do
      write(*,*) ''
  end subroutine write_a_vector
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module ex05module
