
program main
	
	use mpi
	implicit none
	integer :: i,clock,n,rc,ntasks,id
	real :: x
	integer,dimension(:),allocatable :: SEED,status

	integer :: j,c
	integer,parameter :: total=10000,seq_num = 1000
	real(8) :: y,r,pi,sigma
	real(8) :: pi_seq(seq_num)
	real(8) , allocatable :: pi_all(:,:)

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

		allocate(pi_all(ntasks-1,seq_num))
		do i = 1, ntasks - 1
		  call MPI_RECV(pi_seq, seq_num, MPI_DOUBLE_PRECISION, i, i, MPI_COMM_WORLD, status, rc)
		  pi_all(i,:) = pi_seq
		end do

		pi = 0
		do i = 1,ntasks -1
		  do j = 1,seq_num
		    pi = pi + pi_all(i,j)
		  end do
		end do
		pi = pi/dble(seq_num*(ntasks-1))

		write(*,*) 'pi average',pi

		sigma = 0
		do i = 1,ntasks - 1
		  do j = 1,seq_num
		    sigma = sigma + (pi_all(i,j)-pi)**2
		  end do
		end do
		sigma = sigma/dble(seq_num*(ntasks-1))
		write(*,*) 'sigma',sigma

!		write(*,*) pi_all(1,:)
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
		do j = 1,seq_num
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
	   		pi_seq(j) = pi
		end do
!		write(*,*) 'seq',pi_seq
		call MPI_SEND(pi_seq, seq_num , MPI_DOUBLE_PRECISION, 0 , id, MPI_COMM_WORLD,rc)
	end if

!	call MPI_REDUCE(pi,tpi,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,rc)

	call MPI_FINALIZE(rc)

end program
