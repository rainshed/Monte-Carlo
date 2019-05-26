program main
	implicit none
	real(8) :: x,prob,A,xx
	real(8) :: last,next
	integer :: n,i,clock,seed(1)


	xx= 0
	n = 0
	
	last = 0 

	call random_seed()


	do i = 1, 10000000
		call random_number(x)
		next = (x-0.5d0)*200
		if (abs(next)<abs(last)) then
			A = 1
		else
			A = exp(-0.5d0*(next**2-last**2))
		endif
		call random_number(x)
		if (x < A) then
				last = next
			if (i > 1000) then
				xx = xx + last**2
				n = n+1
			endif
		else
			if (i > 1000) then
				xx = xx + last**2
				n = n+1
			endif
		end if
	end do
	
	xx = xx/dble(n)
	write(*,*) xx

end program
