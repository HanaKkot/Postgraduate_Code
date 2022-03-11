!----Week 3 Assignment-----!
!-Steady State heated Rod--!
!-Conjugate gradient method!
PROGRAM CGSolver3D
	IMPLICIT NONE
	
	INTEGER :: kx,ky,kz,il,ih,jl,jh,ml,mh
	INTEGER :: i,niter_precon,j,m
	
	Logical :: precon 
	
	REAL(KIND=8), PARAMETER :: pi = 3.14159265358979323846 
	REAL(KIND=8) :: time1,time2,timetotal
	REAL(KIND=8) :: meshlx,meshly,meshlz,gam,rcurrent,rmax,dx,dy,dz
	REAL(KIND=8), ALLOCATABLE :: aw(:,:,:),ae(:,:,:),ap(:,:,:),an(:,:,:),as(:,:,:),t(:,:,:),b(:,:,:),x(:),y(:),z(:)
	REAL(KIND=8), ALLOCATABLE :: af(:,:,:),ab(:,:,:)
	! Data IO
    INTEGER :: UNIT 
    CHARACTER(LEN = 80) :: filename 

	!##########################################################################################!
	! Initialise MPI here.
	! Temperature field decomposition here.
	! Symbol :: xl ---> local x low index 
	!			xh ---> local x high index
	!			yl ---> local y low index 
	!			yh ---> local y high index
	!##########################################################################################!	

	!##########################################################################################!	
	! If current PID lies on the left-most boundary, then it has one Send/Recv targetting at 
	! MPI_PROC_NULL and vice versa. It is good to have a subroutine to determine the destination
	! of current PID 
	!##########################################################################################!	

	
	!-----------------------------------------
	! USER INPUT------------------------------
	kx = 100; ky = 100; kz = 100 					! MESHSIZE	
	meshlx = 1.0; meshly = 1.0; meshlz = 1.0				! DOMAIN SIZE
	rmax = 0.0001 							! MAX RESIDUAL 	  
  	gam = 100.0 							! HEAT FORCING TERM
	niter_precon = 5 						! NUMBER OF PRECON ITERATIONS
	precon = .TRUE.  						! USE A PRECONDITIONER
	
	
	
	!-----------------------------------------
	! Allocate memory
	CALL AllocateMemory()
	
	
	!-----------------------------------------
	! Initialize Solution Parameters
	WRITE(*,'(a)') "# Initializing Simulation - CG Solver for Steady-State Heat "	
	CALL SolInit()
	WRITE(*,'(a20,i14.1)') "# Nx = ", kx
	WRITE(*,'(a20,i14.1)') "# Ny = ", ky
	WRITE(*,'(a20,i14.1)') "# Nz = ", kz
	WRITE(*,'(a20,f14.1)') "# Lx = ", meshlx
	WRITE(*,'(a20,f14.1)') "# Ly = ", meshly
	WRITE(*,'(a20,f14.1)') "# Lz = ", meshlz
	WRITE(*,'(a20,E20.10)') "# Max Residual = ", rmax
	IF (precon) THEN
		WRITE(*,'(a,i10.1)') "# Number of Jacobi Precon Iterations =  ", niter_precon
  ENDIF
  WRITE(*,'(a)') "# Solution Initialized "
  
  
  !-----------------------------------------
  ! Solver Start
  
	! Start Timer
	CALL CPU_TIME(time1)
	
	! Solve Using CG method
	CALL CGSolve(aw,ae,an,as,ap,af,ab,b,t,rmax)
	
	! Finnish Timer
  	CALL CPU_TIME(time2)
	filename = 'cg_3d.dat'
  	open(unit=2,file= filename)
			write(2,*) 'title="Tecplot"'
			write(2,*) 'variables="x" "y" "z" "Tem"'
			write(2,*) 'zone t="Fortran"  i=' ,kx,' j=',ky,' k=', kz ,' f=point'
		
			do j=1,kx
				do i=1,ky
					do m=1,kz
						write(2,*) j,i,m,real(t(i,j,m))
					end do
				end do
			end do
		close(2)
	
  
	
Contains

!-------------------------------------------------
! Subroutine to allocate memory
SUBROUTINE AllocateMemory()
	IMPLICIT NONE
	
	ALLOCATE(aw(kx,ky,kz))
	ALLOCATE(ae(kx,ky,kz))
	ALLOCATE(an(kx,ky,kz))
  	ALLOCATE(as(kx,ky,kz))
  	ALLOCATE(ap(kx,ky,kz))
	ALLOCATE(af(kx,ky,kz))
	ALLOCATE(ab(kx,ky,kz))
	ALLOCATE(t(kx,ky,kz))
	ALLOCATE(b(kx,ky,kz))
	ALLOCATE(x(kx))
	ALLOCATE(y(ky))
	ALLOCATE(z(kz))

END SUBROUTINE AllocateMemory

!-------------------------------------------------
! Subroutine to Initialize Solution and establish Variables
SUBROUTINE SolInit()
	IMPLICIT NONE
	INTEGER :: i,j,m
	
	!##########################################################################################!
	! When initialising A matrix, make sure that the dimension aligns with domain decomposition!
	! PID that works on the boundary needs to be treated carefully. Depending on numerical     !
	! stencil used, first or last few entries of matrix A maybe different from other entries.  !
	!##########################################################################################!

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
	
	! A and b matrices
	!##########################################################################################!
	! Pay special attention on PID and boundary conditions. We can do smth like this
	!##########################################################################################!
	! IF (PID HAS AT LEAST ONE SIDE on the BOUNDARY) then
	!	b(first few row,first few col) = BC
	! 	ae(first few row,first few col) = coefficient correspond to BC	(Probs 1 for most cases)
	!   aw(first few row,first few col) = coefficient correspond to BC  (Probs 1 for most cases)
	!		........etc...............
	! END IF 
	!##########################################################################################!
	DO j = 1,ky
		DO i = 1,kx
			do m = 1,kz
				b(i,j,m) = gam*x(i)*y(j)*z(m)*sin(x(i)*pi/meshlx)*sin(y(j)*pi/meshly)*sin(z(m)*pi/meshlz)  
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

!-------------------------------------------------
! Conjugent Gradient Solver with Jacobi Preconditioner
SUBROUTINE CGSolve(aw,ae,an,as,ap,af,ab,b,t,rmax)
  IMPLICIT NONE
  REAL(KIND = 8), DIMENSION(kx,ky,kz), INTENT(IN) :: aw,ae,an,as,ap,af,ab,b
  REAL(KIND = 8), DIMENSION(kx,ky,kz), INTENT(INOUT) :: t
  REAL(KIND = 8) :: rmax,rcurrent
  
  REAL(KIND = 8), DIMENSION(kx,ky,kz) :: d,q,s,r,k  
  REAL(KIND = 8) :: delta,delta_o,alpha,beta,gamma
  
  INTEGER :: i,j,iter,jj,m
  
  ! Initialize the solution constants and vectors
  delta = 1.0 ; delta_o = 1.0
  alpha = 1.0; beta = 1.0; gamma = 1.0
  d = 0.0; q(:,:,:) = 0.0; s(:,:,:) = 0.0
  
  ! Iteration Counter
  iter = 0
  
  ! Initialiaze 1/ap Matrix
  DO j = 1,ky
	DO i = 1,kx
		do m = 1,kz
			k(i,j,m) = 1/ap(i,j,m)
		enddo
	ENDDO
  ENDDO
  
  ! Initialize Precondtioned search vector
  !##########################################################################################!
  ! Change looping indexes to xl,xh,yl and yh 
  !##########################################################################################!
  IF (precon) then
  	DO j = jl,jh
			DO i = il,ih
				do m = ml,mh
					d(i,j,m) = r(i,j,m)*k(i,j,m)
				enddo
			ENDDO
  	ENDDO
  ENDIF
  
  ! Get Initial Residual
  !##########################################################################################!
  ! Matrix b here maybe different from each PID (ie. Boundary PID and internal PID)
  !##########################################################################################!
  CALL GetResidual(aw,ae,an,as,ap,af,ab,b,t,r)
  
  ! Calculate Domain averaged residual for stopping critterion
  !##########################################################################################!
  ! Perform MPI_ALLREDUCE here to check for CG residual
  !##########################################################################################!  
  rcurrent = SUM(SUM(SUM(ABS(r(il:ih,jl:jh,ml:mh)),1),1),1) / ((kx-2)*(ky-2))*(kz-2)
  write(*,*) "rcurrent = ",rcurrent
  WRITE(*,'(a20,E20.6)') '# Current Residual = ', rcurrent
  
	IF (precon) then
  	! Get Preconditoned Matrix 'd'  
  	CALL JacobiPrecon(aw,ae,an,as,ap,af,ab,r,d,k)
  ELSE
    ! Or else search vector is the residual
  	d = r
  ENDIF

  ! Get delta and delta_o  
  CALL DotProduct(d,r,delta)
  delta_o = delta
  write(*,*) "delta =",delta
  DO WHILE (rcurrent.GT.rmax)		 
    
  	! Get 'q' Matrix
	call MatrixMultiply(aw,ae,an,as,ap,af,ab,d,q)
  	
  	! Get 'gamma' = d^T*q
	call DotProduct(q,d,gamma)

	write(*,*) "gamma = ",gamma
  	! Get 'alpha'
  	alpha = delta/gamma 
  	write(*,*) "alpha = ",alpha

  	! Update 't' 
  		DO j = jl,jh
			DO i = il,ih
				do m = ml,mh 
					t(i,j,m) = t(i,j,m) + alpha*d(i,j,m)
				enddo
			ENDDO
		ENDDO	

		!##########################################################################################!
		! MPI SEND/RECV here to pass t neighouring PIDs
		! MPI_Waitall after SEND and RECV
		! Example:
		!		MPI_SEND t(xl,:) of PID n to t(xh+1,:) of PID n - 1
		!		MPI_Waitall all request
		!##########################################################################################!

  
  	! Update Residual
  	IF (MOD(iter,50).eq.0) THEN
		call GetResidual(aw,ae,an,as,ap,af,ab,b,t,r)
  	ELSE
		!##########################################################################################!
		! Change looping indexs to xl,xh,yl and yh 
		!##########################################################################################!
  		DO j = jl,jh
				DO i = il,ih 
					do m = ml,mh 
						r(i,j,m) = r(i,j,m) - alpha*q(i,j,m)
					enddo
				ENDDO
			ENDDO		
		ENDIF
		!##########################################################################################!
		! MPI SEND/RECV here to pass r neighouring PIDs
		! MPI_Waitall after SEND and RECV
		! Example:
		!		MPI_SEND t(xl,:) of PID n to t(xh+1,:) of PID n - 1
		!		MPI_Waitall all request
		!##########################################################################################!
		
		IF (precon) THEN
			! Apply a Preconditioner to get 's'
			CALL JacobiPrecon(aw,ae,an,as,ap,af,ab,r,s,k)
		ELSE
			s = r
		ENDIF		
		
		! Update 'delta_o'
		delta_o = delta
		
		! Update 'delta'
		CALL DotProduct(s,r,delta)
		
		! Get Beta
		beta = delta/delta_o

		
		! Update 'd'
		DO j = jl,jh
			DO i = il,ih
				do m = ml,mh 
					d(i,j,m) = s(i,j,m) + beta*d(i,j,m)
				enddo
			ENDDO
		ENDDO	
			
		! Update iteration counter and residual check	
	    !##########################################################################################!
        ! Perform MPI_ALLREDUCE here to check for CG residual
        !##########################################################################################!  	
		rcurrent = SUM(SUM(SUM(ABS(r(il:ih,jl:jh,ml:mh)),1),1),1) / ((kx-2)*(ky-2))*(kz-2)	
  	iter = iter + 1

	  
    WRITE(*,'(a40)') '-----------------------------------------'
  	WRITE(*,'(a20,i20.1)') '# Current Iteration = ', iter
  	WRITE(*,'(a20,E20.6)') '# Current Residual = ', rcurrent
  	WRITE(*,'(a20,f20.10)') '# Current Max Temp = ', MAXVAL(t)

  
  ENDDO
  
END SUBROUTINE CGSOLVE

!-------------------------------------------------
! Subroutine to get Residual Matrix
SUBROUTINE GetResidual(aw,ae,an,as,ap,af,ab,b,x,r)
	IMPLICIT NONE
	REAL(KIND = 8), DIMENSION(kx,ky,kz), INTENT(IN) :: aw,ae,an,as,ap,af,ab,b,x
	REAL(KIND = 8), DIMENSION(kx,ky,kz), INTENT(OUT) :: r
	
	INTEGER :: i,j,m
	
	r(:,:,:) = 0.0
	!##########################################################################################!
	! Change looping indexes to xl,xh,yl and yh 
	!##########################################################################################!
	DO j = jl,jh
		DO i = il,ih
			do m = ml,mh
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
	REAL(KIND = 8), DIMENSION(kx,ky,kz), INTENT(IN) :: aw,ae,an,as,ap,b,k,af,ab
	REAL(KIND = 8), DIMENSION(kx,ky,kz), INTENT(INOUT) :: x

	REAL(KIND = 8), DIMENSION(kx,ky,kz) :: r
	INTEGER :: i,j,m,ii
	!##########################################################################################!
	! Change looping indexes to xl,xh,yl and yh 
	!##########################################################################################!
	DO j = jl,jh
		DO i = il,ih
			do m = ml,mh
				x(i,j,m) = b(i,j,m)*k(i,j,m)
			enddo 
		ENDDO
  ENDDO	

	DO ii = 1,niter_precon
    
		! Get Residual of Current system Ax = b
		CALL GetResidual(aw,ae,an,as,ap,af,ab,b,x,r)
		
		!##########################################################################################!
		! Change looping indexes to xl,xh,yl and yh 
		!##########################################################################################!
		! Update Solution
		DO j = jl,jh
			DO i = il,ih 
				do m = ml,mh
					x(i,j,m) = x(i,j,m) + r(i,j,m)*k(i,j,m)
				enddo 
			ENDDO
		ENDDO	
			
	ENDDO
	
END SUBROUTINE JacobiPrecon
	
!-------------------------------------------------
! Subroutine to get Dot Product of two solution matrices
SUBROUTINE DotProduct(l,m,dp)
	IMPLICIT NONE
	REAL(KIND = 8), DIMENSION(kx,ky,kz), INTENT(IN) :: l,m
	REAL(KIND = 8), INTENT(OUT)	:: dp	
	INTEGER :: i,j,n
	
	dp = 0.0
	
		! Update Dot Product
		!##########################################################################################!
		! We can have openMP here
		!##########################################################################################!
		DO j = jl,jh
			DO i = il,ih
				do n = ml,mh
					dp = dp + l(i,j,n)*m(i,j,n)
				enddo
			ENDDO
		ENDDO	

END SUBROUTINE DotProduct
!-------------------------------------------------
! Subroutine to get Matrix Multiplication of two solution matrices
SUBROUTINE MatrixMultiply(aw,ae,an,as,ap,af,ab,l,m)
	IMPLICIT NONE
	REAL(KIND = 8), DIMENSION(kx,ky,kz), INTENT(IN) :: aw,ae,an,as,ap,l,af,ab
	REAL(KIND = 8), DIMENSION(kx,ky,kz), INTENT(OUT) :: m
	
	INTEGER :: i,j,n
		
		m = 0.0
		
		! Get m = A*l	
		!##########################################################################################!
		! We can have openMP here
		!##########################################################################################!
		DO j = jl,jh
			DO i = il,ih
				do n = ml,mh
					m(i,j,n) = ap(i,j,n)*l(i,j,n) + as(i,j,n)*l(i,j-1,n) + an(i,j,n)*l(i,j+1,n) & 
						+ ae(i,j,n)*l(i+1,j,n) + aw(i,j,n)*l(i-1,j,n) + af(i,j,n)*l(i,j,n-1) + &
						ab(i,j,n)*l(i,j,n+1)
				enddo
			ENDDO
	  ENDDO		 
	

END SUBROUTINE MatrixMultiply

END PROGRAM CGSolver3D
