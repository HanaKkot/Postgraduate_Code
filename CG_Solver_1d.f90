! Major Project
! Programmer: Hana Hu, Narcy Liu
PROGRAM CGSolver3D
	implicit none
	include 'mpif.h'
	integer :: kx,ky,kz,il,ih,jl,jh,ml,mh
	integer :: i,niter_precon,j,m, front_surface, back_surface
	integer :: pid, nprocs, ierr, tag1, tag2, front, back, n_local,index_z
	integer :: domain_size
	integer :: request(4),num_req
	integer :: request2(4),num_req2
	integer :: status(MPI_STATUS_SIZE,4)

	logical :: precon
	real(KIND=8), parameter :: pi = 3.14159265358979323846 
	real(KIND=8) :: time1,time2,timetotal,r0,r_all
	real(KIND=8) :: meshlx,meshly,meshlz,gam,rcurrent,rmax,dx,dy,dz
	real(KIND=8), allocatable :: aw(:,:,:),ae(:,:,:),ap(:,:,:),an(:,:,:),as(:,:,:),t(:,:,:),b(:,:,:),x(:),y(:),z(:)
	real(KIND=8), allocatable :: af(:,:,:),ab(:,:,:),t_final(:,:,:)
	! Data IO
    integer :: UNIT 
    character(len = 80) :: filename 
	num_req = 4
	num_req2 = 4
	!-----------------------------------------
	! USER INPUT------------------------------
	tag1 = 1
	tag2 = 2

	kx = 33; ky = 33 ;kz = 33	 					! 
	meshlx = 1.0; meshly = 1.0; meshlz = 1.0				! DOMAIN SIZE
	rmax = 1E-7	 							! MAX RESIDUAL 	  
  	gam = 100.0 							! HEAT FORCING TERM
	niter_precon = 10 						! NUMBER OF PRECON ITERATIONS
	precon = .TRUE.  						! USE A PRECONDITIONER
	
	! MPI initialize
	call MPI_INIT(ierr) 
    call MPI_COMM_RANK(MPI_COMM_WORLD,pid,ierr) 
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr) 
	write(*,*) "MPI initialize finished"
	! Allocate memory

	! Size of domain in z direction
	n_local = kz/nprocs
	CALL AllocateMemory()
	
	! Define the communication operator
	! If the processor is the processor 0, local data will not swap to the front surface. 
	! If the processor is the last processor, local data will not swap to the back surface. 
	if(pid > 0) then 
        front = pid-1 
    else  
        front = MPI_PROC_NULL 
    end if 
    if(pid < nprocs - 1) then 
        back = pid+1 
    else 
        back = MPI_PROC_NULL 
    end if 
	
	! Define the front index and back index in z direction
	
	! If processor is processor 0, calculation will be start in column 3, end up with column n_local + 1 
	if (pid == 0)then
		front_surface = 3
	else 
		front_surface = 2
	endif
	! If processor is the last processor, calculation will be start in column 2, end up with column n_local 
	if (pid == nprocs - 1) then 
		back_surface = n_local
	else	
		back_surface = n_local + 1
	endif
	! Porcessors in the middle starts from 2, end up with n_local + 1
	
	
	! Initialize Solution Parameters
	WRITE(*,*) "Initializing for 1D Parallel CG Solver for 3D Poission equation "	
	CALL SolInit()

	! Start Timer
	CALL CPU_TIME(time1)
	
	! Solve Using CG method
	CALL CGSolve(aw,ae,an,as,ap,af,ab,b,t,rmax)
	
	! Finnish Timer
  	CALL CPU_TIME(time2)
  	WRITE(*,'(a)') "# Simulation Finnished "
	WRITE(*,'(a15,f14.10)') "# Total WTime = ",  time2 - time1
	
	! Write out full Temperature Field to disk for visualization
    filename = 'cg_3d.dat'
    UNIT = 1001 

	! Output the solution to show the contour plot on Tecplot 360
	if(pid==0) then
		open(unit=2,file= filename)
			write(2,*) 'title="Tecplot"'
			write(2,*) 'variables="x" "y" "z" "Tem"'
			write(2,*) 'zone t="Fortran"  i=' ,kx,' j=',ky,' k=', kz ,' f=point'
		
			do j=1,kx
				do i=1,ky
					do m=1,kz
						write(2,*) j,i,m,real(t_final(i,j,m))
					end do
				end do
			end do
		close(2)
	endif
Contains

!-------------------------------------------------
! Subroutine to allocate memory
SUBROUTINE AllocateMemory()
	IMPLICIT NONE
	
	ALLOCATE(aw(kx,ky,n_local+2))
	ALLOCATE(ae(kx,ky,n_local+2))
	ALLOCATE(an(kx,ky,n_local+2))
  	ALLOCATE(as(kx,ky,n_local+2))
  	ALLOCATE(ap(kx,ky,n_local+2))
	ALLOCATE(af(kx,ky,n_local+2))
	ALLOCATE(ab(kx,ky,n_local+2))
	ALLOCATE(t_final(kx,ky,kz))! 100*100*100  kx=ky=kz
	ALLOCATE(t(kx,ky,n_local+2))
	ALLOCATE(b(kx,ky,n_local+2))
	ALLOCATE(x(kx))
	ALLOCATE(y(ky))
	ALLOCATE(z(kz))

END SUBROUTINE AllocateMemory

!-------------------------------------------------
! Subroutine to Initialize Solution and establish Variables
SUBROUTINE SolInit()
	IMPLICIT NONE
	INTEGER :: i,j,m
	
	! Interior Indices
	il = 2; jl = 2; ml = 2
	ih = kx-1; jh = ky-1; mh = kz-1
	

	! Grid Size
	dx = meshlx/REAL(kx-1)
	dy = meshly/REAL(ky-1)
	dz = meshlz/REAL(kz-1)

  ! x coordinate
	DO i = 1,kx
		x(i) = (i-1)*dx											
	ENDDO
	
	! y coordinate
	DO j = 1,ky
		y(j) = (j-1)*dy											
	ENDDO

	! z coordinate
	DO m = 1,kz
		z(m) = (m-1)*dz											
	ENDDO
	! Initialize solution matrices
	t(:,:,:) = 0.0
	b = 0.0

	! Initialize A matrix and b
	! b is different in each domain, so it should be relevant to the processor id by line 163
	DO j = 1,ky
		DO i = 1,kx
			do m = 2,n_local+1
				index_z = pid*(n_local-1) + m ! b in each domain is different
				b(i,j,m) = gam*x(i)*y(j)*z(index_z)*sin(x(i)*pi/meshlx)*sin(y(j)*pi/meshly)*sin(z(index_z)*pi/meshlz)  
				aw(i,j,m) = -1/dx**2										
				ae(i,j,m) = -1/dx**2		
				an(i,j,m) = -1/dy**2										
				as(i,j,m) = -1/dy**2								
				ap(i,j,m) =  2/dx**2 + 2/dy**2 + 2/dz**2
				af(i,j,m) = -1/dz**2
				ab(i,j,m) = -1/dz**2 
			enddo	
		ENDDO
	ENDDO
	

END SUBROUTINE Solinit


! Conjugent Gradient Solver with Jacobi Preconditioner
SUBROUTINE CGSolve(aw,ae,an,as,ap,af,ab,b,t,rmax)
  IMPLICIT NONE
  REAL(KIND = 8), DIMENSION(kx,ky,n_local+2), INTENT(IN) :: aw,ae,an,as,ap,af,ab,b
  REAL(KIND = 8), DIMENSION(kx,ky,n_local+2), INTENT(INOUT) :: t
  REAL(KIND = 8) :: rmax,rcurrent
  
  REAL(KIND = 8), DIMENSION(kx,ky,n_local+2) :: d,q,s,r,k  
  REAL(KIND = 8) :: delta,delta_o,alpha,beta,gamma, gamma_0, delta_root, delta_all, gamma_all
  
  INTEGER :: i,j,iter,jj,m
  
  ! Initialize the solution constants and vectors
  delta = 1.0 ; delta_o = 1.0
  alpha = 1.0; beta = 1.0; gamma = 1.0
  d = 0.0; q(:,:,:) = 0.0; s(:,:,:) = 0.0
  
  ! Iteration Counter
  iter = 0
  
  ! Initialiaze 1/ap Matrix
  do j = 1,ky
	do i = 1,kx
		do m = 1,n_local+2
			k(i,j,m) = 1/ap(i,j,m)
		enddo
	enddo
  enddo
  
  ! Initialize Precondtioned search vector

  if (precon) then
  	do j = jl,jh
		do i = il,ih
			do m = front_surface,back_surface
				d(i,j,m) = r(i,j,m)*k(i,j,m)
			enddo
		enddo
	enddo
  ENDIF
  
  ! Swap ghost node before caiculating residual
  ! Ghost node of Temperature is used in GetResidual subroutine
  call MPI_ISEND(t(:,:,2),kx*ky,MPI_REAL8,front,tag1,MPI_COMM_WORLD,request(1),ierr)
  call MPI_IRECV(t(:,:,1),kx*ky,MPI_REAL8,front,tag2,MPI_COMM_WORLD,request(2),ierr)
  call MPI_ISEND(t(:,:,n_local+1),kx*ky,MPI_REAL8,back,tag2,MPI_COMM_WORLD,request(4),ierr)
  call MPI_IRECV(t(:,:,n_local+2),kx*ky,MPI_REAL8,back,tag1,MPI_COMM_WORLD,request(3),ierr)

  call MPI_WAITALL(num_req, request(:), status, ierr) 
 
  CALL GetResidual(aw,ae,an,as,ap,af,ab,b,t,r)
  
  ! Calculate residual in each processor
  rcurrent = SUM(SUM(SUM(ABS(r(il:ih,jl:jh,ml:n_local+1)),1),1),1) / ((kx-2)*(ky-2))*(n_local-2)
  
  ! Sum up for each processors and store to r0
  call MPI_REDUCE(rcurrent, r0, 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

  !Brodcast all residual to each processor and do the calculation for the next step
  call MPI_BCAST(r0, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

  r_all = r0/nprocs
  write(*,*) "r_all =",r_all

	IF (precon) then
  	! Get Preconditoned Matrix 'd' 
  		CALL JacobiPrecon(aw,ae,an,as,ap,af,ab,r,d,k)
  ELSE
    ! Or else search vector is the residual
  	d = r
  ENDIF

  ! Get delta and delta_o  
  
  CALL DotProduct(d,r,delta)
 
  ! Reduce local delta to the global
  call MPI_REDUCE(delta, delta_root, 1,MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(delta_root,1,MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
  delta_all = delta_root/nprocs
  delta_o = delta_all
  write(*,*) "delta_all =",delta_all
  
  ! Iteration
  DO WHILE (r_all > rmax)		 
    
  	! Swap ghost node of d
	! Ghost node of d array is used in GetResidual subroutine
	call MPI_ISEND(d(:,:,2),kx*ky,MPI_REAL8,front,tag1,MPI_COMM_WORLD,request(1),ierr)
	call MPI_IRECV(d(:,:,1),kx*ky,MPI_REAL8,front,tag2,MPI_COMM_WORLD,request(2),ierr)
	call MPI_ISEND(d(:,:,n_local+1),kx*ky,MPI_REAL8,back,tag2,MPI_COMM_WORLD,request(4),ierr)
	call MPI_IRECV(d(:,:,n_local+2),kx*ky,MPI_REAL8,back,tag1,MPI_COMM_WORLD,request(3),ierr)

	call MPI_WAITALL(num_req, request(:), status, ierr)
	
	call MatrixMultiply(aw,ae,an,as,ap,af,ab,d,q)
	
  	! Get 'gamma' = d^T*q
	! Reduce local delta to the global
	call DotProduct(q,d,gamma)
	call MPI_REDUCE(gamma, gamma_0, 1,MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
	! Send gamma to each processors
	call MPI_BCAST(gamma_0,1,MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
	! Update the gamma to the global
	gamma_all = gamma_0/nprocs

  	! Get 'alpha'
  	alpha = delta_all/gamma_all
  	! Update 't' 
  		DO j = jl,jh
			DO i = il,ih
				do m = front_surface,back_surface
					t(i,j,m) = t(i,j,m) + alpha*d(i,j,m)
				enddo
			ENDDO
		ENDDO	
		
  	! Update Residual
  	IF (MOD(iter,50).eq.0) THEN
		
		call MPI_ISEND(t(:,:,2),kx*ky,MPI_REAL8,front,tag1,MPI_COMM_WORLD,request(1),ierr)
		call MPI_IRECV(t(:,:,1),kx*ky,MPI_REAL8,front,tag2,MPI_COMM_WORLD,request(2),ierr)
		call MPI_ISEND(t(:,:,n_local+1),kx*ky,MPI_REAL8,back,tag2,MPI_COMM_WORLD,request(4),ierr)
		call MPI_IRECV(t(:,:,n_local+2),kx*ky,MPI_REAL8,back,tag1,MPI_COMM_WORLD,request(3),ierr)
		call MPI_WAITALL(num_req, request(:), status, ierr) 
		call GetResidual(aw,ae,an,as,ap,af,ab,b,t,r)
  	ELSE
  		DO j = jl,jh
				DO i = il,ih 
					do m = front_surface,back_surface
						r(i,j,m) = r(i,j,m) - alpha*q(i,j,m)
					enddo
				ENDDO
			ENDDO		
		ENDIF
		


		IF (precon) THEN
			! Apply a Preconditioner to get 's'
			CALL JacobiPrecon(aw,ae,an,as,ap,af,ab,r,s,k)
		ELSE
			s = r
		ENDIF		

		! Update 'delta_o'
		delta_o = delta_all
		

		! Update 'delta'
		CALL DotProduct(s,r,delta)
		! Reduce local delta to the global
		call MPI_REDUCE(delta, delta_root, 1,MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
		call MPI_BCAST(delta_root,1,MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
		delta_all = delta_root/nprocs
		
		! Get Beta
		beta = delta_all/delta_o

		
		! Update 'd'
		DO j = jl,jh
			DO i = il,ih
				do m = front_surface,back_surface 
					d(i,j,m) = s(i,j,m) + beta*d(i,j,m)
				enddo
			ENDDO
		ENDDO	
			
		! Calculate overall residual
		rcurrent = SUM(SUM(SUM(ABS(r(il:ih,jl:jh,ml:n_local)),1),1),1) / ((kx-2)*(ky-2))*(n_local-2)
		call MPI_REDUCE(rcurrent, r0, 1,MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
		
		!Brodcast all residual to each processor and do the calculation for the next step
		call MPI_BCAST(r0,1,MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
	  
		r_all = r0/nprocs	
  		iter = iter + 1

	  
    	WRITE(*,'(a40)') '-----------------------------------------'
  		WRITE(*,'(a20,i20.1)') '# Current Iteration = ', iter
  		WRITE(*,'(a20,E20.6)') '# Overall Residual = ', r_all
		WRITE(*,'(a20,E20.6)') '# Local Residual   = ', rcurrent
  		WRITE(*,'(a20,f20.10)') '# Current Max Temp = ', MAXVAL(t)

  	ENDDO
  CALL MPI_Gather(t(:,:,2), n_local*kx*ky, MPI_REAL8, t_final(:,:,:), n_local*kx*ky, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
  WRITE(*,'(a20,f20.10)') '# Final Max Temp = ', MAXVAL(t_final)
  call MPI_FINALIZE(ierr)
END SUBROUTINE CGSOLVE

!-------------------------------------------------
! Subroutine to get Residual Matrix
SUBROUTINE GetResidual(aw,ae,an,as,ap,af,ab,b,x,r)
	IMPLICIT NONE
	REAL(KIND = 8), DIMENSION(kx,ky,n_local+2), INTENT(IN) :: aw,ae,an,as,ap,af,ab,b,x
	REAL(KIND = 8), DIMENSION(kx,ky,n_local+2), INTENT(OUT) :: r
	
	INTEGER :: i,j,m
	
	r(:,:,:) = 0.0

	DO j = jl,jh
		DO i = il,ih
			do m = front_surface,back_surface
				r(i,j,m) = b(i,j,m) - (aw(i,j,m)*x(i-1,j,m) + ae(i,j,m)*x(i+1,j,m) &
						 + an(i,j,m)*x(i,j+1,m) + as(i,j,m)*x(i,j-1,m) &
						 + ap(i,j,m)*x(i,j,m) + af(i,j,m)*x(i,j,m-1) + ab(i,j,m)*x(i,j,m+1))
			enddo
		ENDDO
	ENDDO		
END SUBROUTINE GetResidual


!-------------------------------------------------
! Subroutine to get Preconditoner Matrix 'd' using Jacobi
SUBROUTINE JacobiPrecon(aw,ae,an,as,ap,af,ab,b,x,k)
	IMPLICIT NONE
	REAL(KIND = 8), DIMENSION(kx,ky,n_local+2), INTENT(IN) :: aw,ae,an,as,ap,b,k,af,ab
	REAL(KIND = 8), DIMENSION(kx,ky,n_local+2), INTENT(INOUT) :: x

	REAL(KIND = 8), DIMENSION(kx,ky,n_local+2) :: r
	INTEGER :: i,j,m,ii

	DO j = jl,jh
		DO i = il,ih
			do m = front_surface,back_surface
				x(i,j,m) = b(i,j,m)*k(i,j,m)
			enddo 
		ENDDO
  ENDDO	

	DO ii = 1,niter_precon
    
		! Get Residual of Current system Ax = b
		CALL GetResidual(aw,ae,an,as,ap,af,ab,b,x,r)
		
		! Update Solution
		DO j = jl,jh
			DO i = il,ih 
				do m = front_surface,back_surface
					x(i,j,m) = x(i,j,m) + r(i,j,m)*k(i,j,m)
				enddo 
			ENDDO
		ENDDO	
			
	ENDDO
	
END SUBROUTINE JacobiPrecon
	

! Subroutine to get Dot Product of two solution matrices
SUBROUTINE DotProduct(l,m,dp)
	IMPLICIT NONE
	REAL(KIND = 8), DIMENSION(kx,ky,n_local+2), INTENT(IN) :: l,m
	REAL(KIND = 8), INTENT(OUT)	:: dp	
	INTEGER :: i,j,n
	
	dp = 0.0
		! Calculate dot product
		DO j = jl,jh
			DO i = il,ih
				do n = front_surface,back_surface
					dp = dp + l(i,j,n)*m(i,j,n)
				enddo
			ENDDO
		ENDDO	

END SUBROUTINE DotProduct
!-------------------------------------------------
! Subroutine to get Matrix Multiplication of two solution matrices
SUBROUTINE MatrixMultiply(aw,ae,an,as,ap,af,ab,l,m)
	IMPLICIT NONE
	REAL(KIND = 8), DIMENSION(kx,ky,n_local+2), INTENT(IN) :: aw,ae,an,as,ap,l,af,ab
	REAL(KIND = 8), DIMENSION(kx,ky,n_local+2), INTENT(OUT) :: m
	
	INTEGER :: i,j,n
		
		m = 0.0
		
		! Get m = A*l	

		DO j = jl,jh
			DO i = il,ih
				do n = front_surface,back_surface
					m(i,j,n) = ap(i,j,n)*l(i,j,n) + as(i,j,n)*l(i,j-1,n) + an(i,j,n)*l(i,j+1,n) & 
						+ ae(i,j,n)*l(i+1,j,n) + aw(i,j,n)*l(i-1,j,n) + af(i,j,n)*l(i,j,n-1) + &
						ab(i,j,n)*l(i,j,n+1)
				enddo
			ENDDO
	  ENDDO		 
	

END SUBROUTINE MatrixMultiply

END PROGRAM CGSolver3D
