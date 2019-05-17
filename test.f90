program main
	implicit none
	integer,parameter :: pi_num = 1000
	real(8) :: x,y,r,pi,tpi(pi_num),sigma
	integer :: i,j,c,N=1000000
	

	call random_seed()
	do j = 1,pi_num
		c = 0
		do i = 1,N
		   call random_number(x)
		   call random_number(y)
		   r = sqrt(x**2+y**2)
		   if (r<1) then
			   c = c+1
		   endif
		end do
		pi = dble(c)/dble(N)*4
!		write(*,*) pi
		tpi(j) = pi
	end do
	
	pi = 0
	do i = 1,pi_num
		pi = pi + tpi(pi_num)
	end do
	pi = pi/dble(pi_num)
	write(*,*) pi

	sigma = 0
	do i = 1,pi_num
		sigma = sigma + (tpi(i)-pi)**2
	end do
	sigma = sigma/dble(pi_num-1)
!	write(*,*) sigma
	
!	write(*,*) tpi

	

end program
