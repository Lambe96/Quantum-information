program testing

implicit none

integer :: i, j, k

k = 1

do i=1,4
   do j=1,10
      write(*, '(I2,X)', advance='no') k
      k = k + 1
   end do
   write(*, *) ''  ! this gives you the line break
end do

end program testing
