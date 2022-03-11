program assignment1
	!Assignment 1
	!SID : 500633928
	!Name : Yihua Hu
	implicit none
	include 'mpif.h'
	
	!Parameter define
	integer :: i,pid,nprocs,ierr,tag1,tag2, left, right, n_local,j, n, left_column, right_column, index_x, index_y, iter_num
	integer, parameter :: alpha = 50
	integer, parameter :: beta = 4
	integer, parameter :: gamma = 100
	integer, parameter :: Lx = 1 ! Length in x dimention
	integer, parameter :: Ly = 1 ! Length in y dimention
	integer, parameter :: total_size = 100 ! matrix size
	! For the memory allocation in line 18,19, should be changed if the number of processor is changed, local_size = n_local and should be an integer
	integer, parameter :: local_size = 25 
	integer, parameter :: pi = 3.14159265358979323846! pi
	integer :: status(MPI_STATUS_SIZE,4) 
	integer :: request(4),num_req
	real, parameter :: criterion = 1e-4
	
	!Memory allocation
	real(kind = 8) :: Tem(total_size, local_size + 2)! Old temperature, +2 is for the ghost node on left and right
	real(kind = 8) :: T(total_size, local_size + 2) ! New Temperature, +2 is for the ghost node on left and right
	real(kind = 8) :: R, R0, R_all ! Residual
	real(kind = 8) :: T_final(total_size, total_size) ! Final global temperature field
	real(kind = 8) :: U(total_size, total_size) ! U table
	real(kind = 8) :: V(total_size, total_size) ! V table
	real(kind = 8) :: delta_x, delta_y ! Grid size
	
	! Calculate delta_x
	delta_x = Lx/real((total_size-1))
	! Let delta_y = delta_x
	delta_y = delta_x 
	!Communication tag
	tag1 = 1
	tag2 = 2
	!Number of request
	num_req = 4
	!Overall residual
	R_all = criterion + 1
	iter_num = 0
	!Define table of U and V in all the temperature point
	!Because the index of row in Fortran is from top to the bottom, so the relationship between index and y should be transferred by certain algorithm
	do j=1,total_size	!x index
		do i=1,total_size !y index
			U(i, j) = COS(beta*2*pi*((j-1)*delta_x)) * SIN(beta*2*pi*(total_size-i)*delta_y)
			V(i, j) =-SIN(beta*2*pi*((j-1)*delta_x)) * COS(beta*2*pi*(total_size-i)*delta_y)
		enddo
	enddo
		
	! MPI initialize
	call MPI_INIT(ierr) 
    call MPI_COMM_RANK(MPI_COMM_WORLD,pid,ierr) 
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr) 
	
	! Judge n_local is an integer
	if (mod(total_size,nprocs).ne.0) then
		write ( *, * ) 'total_size/nprocs should be zero'
		stop
	endif
	
	! Size of column in each process
	n_local = total_size/nprocs
	
	!Initialize the temperature with all zero 
	do j=1,n_local+2 
        do i=1,total_size 
            Tem(i,j)=0.0 
        end do 
    end do 

	!Define the communication, if the process is the process 0, it will not send data to the left, and the right process will not send
	!data to the right
	if(pid > 0) then 
        left=pid-1 
    else 
        left=MPI_PROC_NULL 
    end if 
    if(pid < nprocs - 1) then 
        right=pid+1 
    else 
        right=MPI_PROC_NULL 
    end if 
	
	!Iteration, if residual is bigger than criterion, loop will not end unitil R_all<criterion
	do while(R_all>criterion)
		iter_num = iter_num + 1
		!Using MPI non-blocking communication to swap ghost node
		!Data in the send/recv buffer is the pointer towards the start point of the array
		!e.g. Tem(1,2) is the data in first row and second column, but the length is equal to total size, so it will send second column
		call MPI_ISEND(Tem(1,2),total_size,MPI_REAL8,left,tag1,MPI_COMM_WORLD,request(1),ierr)
		call MPI_IRECV(Tem(1,1),total_size,MPI_REAL8,left,tag2,MPI_COMM_WORLD,request(2),ierr)
		call MPI_IRECV(Tem(1,n_local+2),total_size,MPI_REAL8,right,tag1,MPI_COMM_WORLD,request(3),ierr)
		call MPI_ISEND(Tem(1,n_local+1),total_size,MPI_REAL8,right,tag2,MPI_COMM_WORLD,request(4),ierr)
		!Synchronization for all the request
		call MPI_WAITALL(num_req, request(:), status, ierr)
		! Define the calculation area, if processor 0, calculation will be start in column 3,because row 1 is for the ghost node, row 2 is the boundary
		! with zero temperature
		if (pid == 0)then
			left_column = 3
		else 
			left_column = 2
		endif
		! If processor at the end of right, calculation will be start in column 2, end up with column(length of temperature in each processor)
		! because the column(len(temperature))+1 is the boudary condition, and column(len(temperature))+2 is the ghost node.
		if (pid == nprocs - 1) then 
			right_column = n_local
		else
			right_column = n_local + 1
		endif
		
		do j = left_column, right_column
			do i = 2, total_size-1
				! Index in each processor
				index_x = pid * n_local + j
				index_y = i
				! Formula in the report by using jacobi solver
				T(i, j) = (2*gamma*delta_x**2+(2-alpha*U(index_y, index_x)*delta_x)*Tem(i,j+1)+&
											  (2+alpha*U(index_y, index_x)*delta_x)*Tem(i,j-1)+&
											  (2-alpha*V(index_y, index_x)*delta_x)*Tem(i+1,j)+&
											  (2+alpha*V(index_y, index_x)*delta_x)*Tem(i-1,j))/8
			end do
		end do
		
		!Boundary condition with all zero in 4 edges
		T(1, :) = 0
		T(total_size, :) = 0
		if (pid == 0) then
			T(:, left_column - 1) = 0
		elseif (pid == nprocs-1) then
			T(:, right_column + 1) = 0
		else
		endif
		
		!Calculate all the residual in each processor
		R = sum(abs(T(:,left_column:right_column)-Tem(:,left_column:right_column)))
		
		!Send each residual to the root(0 processor) and sum them up
		call MPI_REDUCE(R, R0, 1,MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
		
		!Brodcast all residual to each processor and do the calculation for the next step
		call MPI_BCAST(R0,1,MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

		!Overall residual should be divided by the processor number
		R_all = R0/nprocs
		
		!Swap data for the next iteration
		Tem(:,:) = T(:,:)

	end do


	!Gather all temperature to the final temperature matrix in root processor
	CALL MPI_Gather(Tem(1,2), n_local*total_size, MPI_REAL8, T_final(:,:), n_local*total_size, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
	
	open(unit=2,file="Temp.dat")
		write(2,*) 'title="Tecplot"'
		write(2,*) 'variables="x" "y" "Tem"'
		write(2,*) 'zone t="Fortran"  j=',total_size,' i=',total_size,' f=point'
    do j=1,total_size
       do i=1,total_size
			write(2,*) j,i,T_final(i,j)
       end do
    end do
	close(2)

	!End MPI
	call MPI_FINALIZE(ierr)
	
end program assignment1
