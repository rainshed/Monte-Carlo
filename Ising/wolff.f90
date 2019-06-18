program main
	
	use mpi
	implicit none
	integer :: i,j,k,clock,n,rc,ntasks,id
	real :: x
	integer,dimension(:),allocatable :: SEED,status
!--------------------------------------------------------------

	integer,parameter :: D = 32
	integer :: state(D,D), r1, r2, neighbor_sum
	integer :: warm_up = 20, total = 100
	real(8) ::  eps = -1, E1, E2, p,T=1


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
		


!------------------------------------------------------------------------------------------------------------

	else
		call MPI_RECV(clock, 1, MPI_LONG, 0, id, MPI_COMM_WORLD, status, rc)
		call RANDOM_SEED(size=n)
		allocate(SEED(n))
		do i = 1,n
		  SEED(i) = clock + 37*i
		end do
		call RANDOM_SEED(PUT=SEED)
		deallocate(SEED)
		call RANDOM_NUMBER(x)

!------------------------------------------------------------------------------------------------
		do k = 1,500
		T = dble(k)/dble(100)

		!-----initial----------------------
!		do i = 1,D
!			do j = 1,D
!				call random_number(x)
!				if (x<0.5) then
!					state(i,j) = 1
!				else 
!					state(i,j) = -1
!				endif
!			end do
!		end do
		state = 1
		!----------------------------------


		!-----dynamic----------------------
		do i = 1,warm_up+total
			call random_number(x)
			r1 = int(x*D) + 1
			call random_number(x)
			r2 = int(x*D) + 1
			call nextfilp(state,r1,r2,T,D)

		end do
		open(unit=11,file='test.dat',status='unknown')
		
		write(11,*) T,abs(sum(state))
		end do 


		!-----------------------------------
				


			



!-----------------------------------------------------------------------------------------------

	end if

!	call MPI_REDUCE(pi,tpi,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,rc)

	call MPI_FINALIZE(rc)

end program

recursive subroutine nextfilp(state,i,j,T,D)
	integer :: D
	integer :: state(D,D),i,j
	real(8) :: p,x,T
	p = 1-exp(-2/T)
	if (state(i,j) == 1) then
		state(i,j) = -1
	else
		state(i,j) = 1
	endif
	if (state(i,j) /= state(mod(i,D)+1,j)) then
!		write(*,*) state(i,j),state(mod(i,D)+1,j)
		call random_number(x)
		if (x<p) then
		!	state(mod(i,D)+1,j) = state(i,j)
			call nextfilp(state,mod(i,D)+1,j,T,D)
		endif
	endif
	if (state(i,j) /= state(mod(i+D-1,D),j)) then
		call random_number(x)
		if (x<p) then
		!	state(mod(i+D-1,D),j) = state(i,j)
			call nextfilp(state,mod(i+D-1,D),j,T,D)
		endif
	endif
	if (state(i,j) /= state(i,mod(j,D)+1)) then
		call random_number(x)
		if (x<p) then
		!	state(i,mod(j,D)+1) = state(i,j)
			call nextfilp(state,i,mod(j,D)+1,T,D)
		endif
	endif
	if (state(i,j) /= state(i,mod(j+D-1,D))) then
		call random_number(x)
		if (x<p) then
		!	state(i,mod(j+D-1,D)) = state(i,j)
			call nextfilp(state,i,mod(j+D-1,D),T,D)
		endif
	endif

end subroutine

!	recursive function nextfilp(state,i,j,T,D)
!		integer :: D
!		integer :: state(D,D),i,j
!		real(8) :: p,x,T
!		p = 1-exp(-2/T)
!		if (state(i,j) == 1) then
!			state(i,j) = -1
!		else
!			state(i,j) = 1
!		endif
!		if (state(i,j) /= state(mod(i,D)+1,j)) then
!			call random_number(x)
!			if (x<p) then
!				state(mod(i,D)+1,j) = state(i,j)
!				nextfilp(state,mod(i,D)+1,j,T,D)
!			endif
!		endif
!		if (state(i,j) /= state(mod(i+D-1,D),j)) then
!			call random_number(x)
!			if (x<p) then
!				state(mod(i+D-1,D),j) = state(i,j)
!				nextfilp(state,mod(i+D-1,D),j,T,D)
!			endif
!		endif
!		if (state(i,j) /= state(i,mod(j,D)+1)) then
!			call random_number(x)
!			if (x<p) then
!				state(i,mod(j,D)+1) = state(i,j)
!				nextfilp(state,i,mod(j,D)+1,T,D)
!			endif
!		endif
!		if (state(i,j) /= state(i,mod(j+D-1,D))) then
!			call random_number(x)
!			if (x<p) then
!				state(i,mod(j+D-1,D)) = state(i,j)
!				nextfilp(state,i,mod(j+D-1,D),T,D)
!			endif
!		endif
!	
!	end function
