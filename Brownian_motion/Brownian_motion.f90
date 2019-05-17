module fun
	implicit none
	integer,parameter :: Nx = 101, Ny = 101
contains

	function initial() result(xy)
		integer :: xy(2)
		real(8) :: x
		call random_number(x)
		xy(1) = int(x*Nx+1)
		call random_number(x)
		xy(2) = int(x*Ny+1)
	end function

	function go(input) result(output)
		integer :: input(2),output(2)
		real(8) :: x,a
		call random_number(x)
		output(1) = Mod(input(1) + int(x*3) -1+Nx, Nx)
		call random_number(x)
		output(2) = Mod(input(2) + int(x*3) -1+Ny, Ny)
	end function


end module fun

program main
	use fun
	implicit none
	integer :: pos(Nx,Ny),init(2),next(2)
	integer :: i,j


	pos = 0
	pos((Nx+1)/2,(Ny+1)/2) = 1
	call random_seed()
	do i = 0,300
	init = initial()
	if (pos(init(1),init(2)) == 0) then
		do while (.true.)
			next = go(init)
!			write(*,*) next
			if (pos(next(1),next(2)) == 1) then
				pos(init(1),init(2)) = 1
				write(*,*)init
				exit
			endif
			init = next
		end do
	endif
	end do

	open(unit = 10,file='data.dat',status='unknown')
	do i = 1,Nx
	  do j = 1,Ny
	    write(10,*) i,j,pos(i,j)
	  end do
	end do


end program	
