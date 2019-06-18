program main
	
	use mpi
	implicit none
	integer :: i,j,k,clock,n,rc,ntasks,id
	real :: x
	integer,dimension(:),allocatable :: SEED,status
!--------------------------------------------------------------

	integer,parameter :: D = 32
	integer :: state(D,D), r1, r2, neighbor_sum
	integer :: warm_up = 100, total = 1000
	real(8) ::  eps = -1, E1, E2, p,T=1
	real(8) :: m,m1,m2,chi,c,energy,energy2,energy_t
	character(2) :: num


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
	endif

!------------------------------------------------------------------------------------------------
		state = 1



		do k = id,500,ntasks
		T = dble(k)/dble(100)
		if (id==1) then
			write(*,*) k/5,'%'
		endif

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
		!----------------------------------

		m = 0
		m1 = 0
		m2 = 0
		energy_t = 0
		energy2 = 0
		c = 0

		!-----dynamic---------------------- 
		do i = 1,warm_up+total
			call random_number(x)
			r1 = int(x*D) + 1
			call random_number(x)
			r2 = int(x*D) + 1
			call nextfilp(state,r1,r2,T,D)
			if (i>warm_up) then
				m = m + dble(abs(sum(state)))/D/D
				m2 = m2 + dble(abs(sum(state)))**2/D/D/D/D
				energy = 0
				do j = 1,D
				do n = 1,D
				neighbor_sum = dble(state(mod(j+D-2,D)+1,n)+state(mod(j,D)+1,n)+state(j,mod(n-2+D,D)+1)+state(j,mod(n,D)+1))
!				write(*,*) mod(j-2+D,D)+1,n,state(mod(j-2+D,D)+1,n)
				energy = energy- dble(state(j,n))*neighbor_sum
				end do
				end do
!				energy = energy/D/D
				energy2 = energy2 + energy**2 
				energy_t = energy_t + energy
			endif
		end do
!		if (T>2.6) then 
!		write(*,*) energy_t,energy2
!		stop
!		endif
		m = m/total
		m2 = m2/total
		chi = m2-m**2
		energy_t = energy_t/total
		energy2 = energy2/total
		c = energy2 - energy_t**2


		write(num,'(i2)') id
		open(unit=id,file='wolff/'//'sp'//trim(adjustl(num))//'M_T_s32.dat',status='unknown')
		open(unit=id+30,file='wolff/'//'sp'//trim(adjustl(num))//'chi_T_s32.dat',status='unknown')
		open(unit=id+100,file='wolff/'//'sp'//trim(adjustl(num))//'C_T_s32.dat',status='unknown')
		write(id,*) T,m
		write(id+30,*) T,chi
		write(id+100,*) T,c

	
		write(11,*) T,m
		write(12,*) T,chi
		write(13,*) T,c


		end do 

!		open(unit=10,file='Ising.dat',status='unknown')
!		do i = 1,D
!			do j = 1,D
!			write(10,*) i,j,state(i,j)
!			end do
!		end do
		!-----------------------------------
				


			



!-----------------------------------------------------------------------------------------------


!	call MPI_REDUCE(pi,tpi,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,rc)

	call MPI_FINALIZE(rc)

	call execute_command_line("cat wolff/sp*C*s32.dat > data/C_T_s32.dat")
	call execute_command_line("cat wolff/sp*chi*s32.dat > data/chi_T_s32.dat")
	call execute_command_line("cat wolff/sp*M*s32.dat > data/M_T_s32.dat")
	call execute_command_line("rm wolff/sp*")

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
		call random_number(x)
		if (x<p) then
		!	state(mod(i,D)+1,j) = state(i,j)
			call nextfilp(state,mod(i,D)+1,j,T,D)
		endif
	endif
	if (state(i,j) /= state(mod(i+D-2,D)+1,j)) then
		call random_number(x)
		if (x<p) then
		!	state(mod(i+D-1,D),j) = state(i,j)
			call nextfilp(state,mod(i+D-2,D)+1,j,T,D)
		endif
	endif
	if (state(i,j) /= state(i,mod(j,D)+1)) then
		call random_number(x)
		if (x<p) then
		!	state(i,mod(j,D)+1) = state(i,j)
			call nextfilp(state,i,mod(j,D)+1,T,D)
		endif
	endif
	if (state(i,j) /= state(i,mod(j+D-2,D)+1)) then
		call random_number(x)
		if (x<p) then
		!	state(i,mod(j+D-1,D)) = state(i,j)
			call nextfilp(state,i,mod(j+D-2,D)+1,T,D)
		endif
	endif

end subroutine
