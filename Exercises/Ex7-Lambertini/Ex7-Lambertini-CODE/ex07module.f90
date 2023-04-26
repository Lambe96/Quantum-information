module ex07module
  implicit none

  contains

    !!!!!!!!!!!!!!!!!!! Initialize a tridiagonal matrix !!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!that describes properly the q.a.o.!!!!!!!!!!!!!!!!!!!!
    subroutine init_tridiag_qao_matrix(matrix, x_max,xx, N, freq, mass, h_bar)
      complex(kind=8), dimension(:,:), allocatable :: matrix
      real(kind=8), dimension(:), allocatable :: xx
      real(kind=8) :: Re_z, x_max,step,freq,mass,h_bar
      integer(kind=8) ::  ii, N

      allocate(xx(N))
      allocate(matrix(N,N))

      step = 2*x_max/real(N-1,8)

      do ii=1,N
        xx(ii) = -x_max+(ii-1)*step
      end do

      !initialize the diagonal
      do ii=1,N
        Re_z = ((2*h_bar**2)/(2*mass) + (mass/2)*freq**2*xx(ii)**2*step**2)/step**2
        matrix(ii,ii) = cmplx(Re_z,0)
      end do

      Re_z = (-h_bar**2/(2*mass)/step**2)

      !initialize the lower elements
      do ii =1,N-1
        matrix(ii,ii+1) = cmplx(Re_z,0)
      end do

      !initialize the upper elements
      do ii = 2,N
        matrix(ii,ii-1) = cmplx(Re_z,0)
      end do

    end subroutine init_tridiag_qao_matrix
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

    !!!!!!!!!!!!!!!!!!!!! Temporal evolution !!!!!!!!!!!!!!!!!!!!
    subroutine temporal_evolution(psi, t_max, t_resolution, x_grid, omega, mass, h_bar, psi_evol)
      !declaration
      complex(kind=8), dimension(:,:), allocatable :: psi_evol
      complex(kind=8) :: psi(:)
      integer(kind=8) :: t_resolution, plan_forward, plan_backward, ii, jj
      real(kind=8) :: t_max, x_grid(:), t_grid(t_resolution), delta_t, delta_x, omega, mass, h_bar, area

      allocate(psi_evol(size(x_grid),t_resolution+1))

      !time grid initialization
      delta_t = t_max/real(t_resolution-1,8)
      delta_x = 2*maxval(x_grid)/real(size(x_grid),8)

      psi_evol(:,1) = psi
      do ii=1,t_resolution
        t_grid(ii) = 0+(ii-1)*delta_t
      end do

      !prepare the system for the fourier transform
      call dfftw_plan_dft_1d(plan_forward,size(x_grid),psi,psi,-1,64);
      call dfftw_plan_dft_1d(plan_backward,size(x_grid),psi,psi,1,64);


      psi_evol(:,1) = psi
      do ii = 1,t_resolution

        !propagate the first half of the potential
        do jj = 1, size(x_grid)
          psi(jj) = psi(jj) * exp(cmplx(0,-delta_t/2*0.5*mass*omega**2/h_bar)*(x_grid(jj)-t_grid(ii)/t_max)**2)
        end do

        !momentum space
        call dfftw_execute_dft(plan_forward,psi,psi)
        !propagate the kinetic term

        do jj=0, size(x_grid)/2-1
          psi(jj+1)=psi(jj+1)*exp(cmplx(0, -0.5*delta_t*4*acos(-1.0)**2/&
          (h_bar*mass*(delta_x*size(x_grid))**2))*jj**2)
        end do

        do jj=size(x_grid)/2, size(x_grid)-1
          psi(jj+1)=psi(jj+1)*exp(cmplx(0, -0.5*delta_t*4*acos(-1.0)**2/&
          (h_bar*mass*(delta_x*size(x_grid))**2))*(size(x_grid)-jj)**2)
        end do

        !back to the position space
        call dfftw_execute_dft(plan_backward,psi,psi)

        !propagate the second half of the potential
        do jj = 1, size(x_grid)
          psi(jj) = psi(jj) * exp(cmplx(0,-delta_t/2*0.5*mass*omega**2/h_bar)*(x_grid(jj)-t_grid(ii)/t_max)**2)
        end do

        !normalize the evolved state
        psi = psi/real(size(x_grid),8)

        !store the evolved state in a matrix
        psi_evol(:,ii+1) = psi
      end do
      call dfftw_destroy_plan(plan_forward)
      call dfftw_destroy_plan(plan_backward)
    end subroutine temporal_evolution
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!! compute the squared modulus of the  !!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!! elements of a matrix  !!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mod_squared_matrix(mm)
      complex(kind=8), dimension(:,:) :: mm
      integer(kind=8), dimension(2) :: dimm
      integer(kind=8) :: ii,jj

      dimm = shape(mm)

      do ii=1,dimm(1)
        do jj=1,dimm(2)
          mm(ii,jj) = mm(ii,jj)*conjg(mm(ii,jj))
        end do
      end do
    end subroutine mod_squared_matrix
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!! compute the squared modulus of the  !!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!! elements of a matrix  !!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mod_squared_array(vv)
      complex(kind=8), dimension(:) :: vv
      integer(kind=8) :: ii

      do ii=1,size(vv)
        vv(ii) = vv(ii)*conjg(vv(ii))
      end do
    end subroutine mod_squared_array
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


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
        write(*,'(F10.3)') mm(i)
      end do
        write(*,*) ''
    end subroutine write_a_vector
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!! write up to 3 real arrays in a file !!!!!!!!!!!!!!!!!!!!!
    subroutine write_arrays_on_file(title,x,y,z)
      character(len=*) title
      real(kind=8), dimension(:) :: x
      real(kind=8), optional, dimension(:) :: y,z
      integer :: i

      if (present(y) .and. present(z)) then

        if ((size(x).eq.size(y)) .and. (size(x).eq.size(z))) then

          open(1,file=title,Access = 'append')
          do i=1,size(x)
            write(1,'(F0.15,A,F0.15,A,F0.15)') x(i),',', y(i), ',', z(i)
          end do
          close(1)

        else

          write(*,*) 'arrays have not the same size'

        end if

      elseif (present(y)) then

        if (size(x).eq.size(y)) then

          open(1,file=title,Access = 'append')
          do i=1,size(x)
            write(1,'(F0.15,A,F0.15)') x(i),',', y(i)
          end do
          close(1)

        else

          write(*,*) 'arrays have not the same size'

        end if

      else

        open(1,file=title,Access = 'append')
        do i=1,size(x)
          write(1,'(F0.15)') x(i)
        end do
        close(1)

      end if
    end subroutine write_arrays_on_file
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !!!!!!!!!!! write the first k entries of a real array in a file !!!!!!!!!!!
    subroutine first_k_entries_on_file(title,vec,kk)
      character(len=*) title
      real(kind=8), dimension(:) :: vec
      integer(kind=8) :: ii, kk

      if ( kk .le. size(vec) ) then

        open(1,file=title,Access = 'append')
        do ii=1,kk
          write(1,'(I0,A,F0.15)') (ii-1),',',vec(ii)
        end do
        close(1)
      else
        !error handling
        write(*,*) 'size(vec) > k'

      end if
    end subroutine first_k_entries_on_file
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!! write up to 2 real arrays in a file and the first k !!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!! columns of a matrix in a file !!!!!!!!!!!!!!!!!!!!!!
    subroutine vec_and_col_on_file(title,mm,kk,vec_1,vec_2)
      character(len=*) title
      real(kind=8), optional, dimension(:) :: vec_1,vec_2
      real(kind=8), dimension(:,:) :: mm
      integer(kind=8) :: ii, jj, kk

      !both vec_1 and vec_2 are present and written as 1st and 2nd columns
      if ( present(vec_2) .and. present(vec_1) ) then

        if ( (size(vec_2) .eq. size(vec_1) ) .and. (size(vec_2) .eq. size(mm, dim=1) )) then

          open(1,file=title,Access = 'append')
          do ii=1,size(vec_1)
            write(1,*) vec_1(ii), ',', vec_2(ii), ',', (mm(ii,jj),',',jj = 1, kk)
          end do
          close(1)


        else
          !error handling
          write(*,*) 'size mismatch'

        end if

      !only vec_1 is present and written as 1st column
      elseif( present(vec_1) ) then
        if ( size(vec_1) .eq. size(mm, dim=1)) then

          open(1,file=title,Access = 'append')
          do ii=1,size(vec_1)
            write(1,*) vec_1(ii),',', (mm(ii,jj),',',jj = 1, kk)
          end do
          close(1)

        else
          !error handling
          write(*,*) 'size mismatch'
        end if

      !only the matrix is present and the first k columns are written
      else

        open(1,file=title,Access = 'append')
        do ii=1,size(vec_1)
          write(1,*) (mm(ii,jj),',',jj = 1, kk)
        end do
        close(1)
      end if

    end subroutine vec_and_col_on_file






  end  module ex07module
