program main
	
	use mpi
	implicit none
	integer :: i,j,k,clock,n,rc,ntasks,id
	real :: x
	integer,dimension(:),allocatable :: SEED,status
!--------------------------------------------------------------

	integer,parameter :: D = 64
	integer :: state(D,D), r1, r2, neighbor_sum
	integer :: warm_up = 100000, total = 10000000
	real(8) ::  eps = -1, E1, E2, p,T=1
	integer :: fact


	call MPI_INIT(rc)
	call MPI_COMM_SIZE(MPI_COMM_WORLD,ntasks,rc)
	call MPI_COMM_RANK(MPI_COMM_WORLD,id,rc)
	allocate(Status(MPI_STATUS_SIZE))

	if (id == 0) then
		call SYSTEM_CLOCK(clock)
		call RANDOM_SEED(size=n)
		allocate(SEED(n))
		do i = 1,n
		  SEED(i) = clock + 37*i
		end do
		call RANDOM_SEED(PUT=SEED)
		deallocate(SEED)
		call RANDOM_NUMBER(x)

		do i = 1, ntasks - 1
		  call RANDOM_NUMBER(x)
		  clock = clock + int(x*1000000)
		  call MPI_SEND(clock, 1 , MPI_LONG, i , i, MPI_COMM_WORLD,rc)
		end do

!-----------------------------------------------------------------------------------------------------
	call fact(2,i)
	write(*,*) i	


!------------------------------------------------------------------------------------------------------------

	else
		!-----------------------------------
				


			



!-----------------------------------------------------------------------------------------------

	end if

!	call MPI_REDUCE(pi,tpi,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,rc)

	call MPI_FINALIZE(rc)

end program

	recursive subroutine fact(n,ans) 
		integer :: n,a
		integer :: ans
		if (n == 0) then
			ans = 1
		else
			call fact(n-1,a)	
			ans = n*a
		endif
	end
