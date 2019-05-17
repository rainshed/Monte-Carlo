!program main
!	implicit none
!	real(8) :: x,y,r,pi
!	integer :: i,c,N=10000
!	
!
!	call random_seed()
!	c = 0
!	do i = 1,N
!	   call random_number(x)
!	   call random_number(y)
!	   r = sqrt(x**2+y**2)
!	   if (r<1) then
!		   c = c+1
!	   endif
!	end do
!	pi = dble(c)/dble(N)
!	write(*,*) pi
!
!	
!
!end program

program main
	
	use mpi
	implicit none
	integer :: i,clock,n,rc,ntasks,id
	real :: x
	integer,dimension(:),allocatable :: SEED,status

	integer :: j,c,seq_num 
	integer,parameter :: total=100000
	real(8) :: y,r,pi
	real(8),allocatable :: pi_seq(:)
	real(8) , allocatable :: pi_all(:,:)

	call MPI_INIT(rc)
	call MPI_COMM_SIZE(MPI_COMM_WORLD,ntasks,rc)
	call MPI_COMM_RANK(MPI_COMM_WORLD,id,rc)
	allocate(Status(MPI_STATUS_SIZE))
	seq_num = ntasks -1
	allocate(pi_seq(ntasks-1))

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
		write(*,*) id,x

		do i = 1, ntasks - 1
		  call RANDOM_NUMBER(x)
		  clock = clock + int(x*1000000)
		  call MPI_SEND(clock, 1 , MPI_LONG, i , i, MPI_COMM_WORLD,rc)
		end do

		allocate(pi_all(ntasks-1,seq_num))
		do i = 1, ntasks - 1
		  call MPI_RECV(pi, 1, MPI_DOUBLE_PRECISION, i, i+ntasks, MPI_COMM_WORLD, status, rc)
		  write(*,*) pi
		  pi_seq(i) = pi
		end do
		
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
!		do j = 1,seq_num
			c = 0
	   		do i = 1,total
	   		   call random_number(x)
	   		   call random_number(y)
	   		   r = sqrt(x**2+y**2)
	   		   if (r<1) then
	   			   c = c+1
	   		   endif
	   		end do
	   		pi = dble(c)/dble(total)*4
!	   		pi_seq(j) = pi
!		end do
!		write(*,*) pi
		call MPI_SEND(pi, 1 , MPI_DOUBLE_PRECISION, 0 , id+ntasks, MPI_COMM_WORLD,rc)
	end if

!	call MPI_REDUCE(pi,tpi,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,rc)

	call MPI_FINALIZE(rc)

end program
