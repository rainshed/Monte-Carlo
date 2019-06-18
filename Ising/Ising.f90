program main
	
	use mpi
	implicit none
	integer :: i,j,k,clock,n,rc,ntasks,id
	real :: x
	integer,dimension(:),allocatable :: SEED,status
!--------------------------------------------------------------

	integer,parameter :: D = 32
	integer :: state(D,D), r1, r2, neighbor_sum
	integer :: warm_up = 100000, total = 100000
	real(8) ::  eps = -1, E1, E2, p,T=1
	real(8) :: m,m1,m2,chi,c,energy,energy2,energy_t


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
			neighbor_sum = state(mod(r1-2+D,D)+1,r2)+state(mod(r1,D)+1,r2)+state(r1,mod(r2-2+D,D)+1)+state(r1,mod(r2,D)+1)
			E1 =dble(eps*state(r1,r2)*neighbor_sum)
			E2 = -E1
			if (E2 < E1) then
				state(r1,r2) = -state(r1,r2)
			else
				p = exp(-(E2-E1)/T)
			!	write(*,*) p
			!	stop
				call random_number(x)
				if (x<p) then
					state(r1,r2) = -state(r1,r2)
				endif
			endif
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
		open(unit=11,file='data/M_T_s32.dat',status='unknown')
		open(unit=12,file='data/chi_T_s32.dat',status='unknown')
		open(unit=13,file='data/C_T_s32.dat',status='unknown')
		
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

end program
