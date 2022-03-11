program assignment2

! Assignment 1
! SID : 500633928
! Name : Yihua Hu
! In order to compile this file, use command gfortran -o exe Assignment2.f90 -I/usr/include/ -lfftw3, and ./exe
use, intrinsic :: iso_c_binding
implicit none
include 'fftw3.f03'

!--------------------------------------------------------------------------------------------------------
integer(C_INT), parameter :: Nx = 64
integer(C_INT), parameter :: Ny = 64
integer :: i,j,k, snapshot

type(C_PTR) :: planA,planB,planC,planD,planE,planF
complex(C_DOUBLE_COMPLEX), dimension(Nx, Ny) :: T, Tf, c_kx, c_ky, Tf1, Tf2,T1,T2, T_temp1, T_temp2, T_temp1f, T_temp2f

real(kind = 8) :: U(Nx, Ny) ! U table
real(kind = 8) :: V(Nx, Ny) ! V table
real(kind = 8) :: dx,dy,M,Lx,Ly,tmax,dt,time,alpha,kx1,ky1,start,end
real(kind = 8) :: kx(1:Nx),ky(1:Ny)
real(kind = 8), dimension(Nx, Ny) :: kxx, kyy
real(kind = 8), parameter		:: PI = 4.*ATAN(1.0d0)

call cpu_time(start) ! wall time
Lx = 2.0*pi ! Length
Ly = Lx ! Height
tmax = 10.0 ! Maximum time
dt = 0.005 ! Time resolution
time = 0.0 
alpha = 0.005
dx = Lx/Nx ! Grid size
dy = Ly/Ny ! Grid size dy=dx
snapshot = 0
! Same index algorithm in Assignment1
do j = 1, Nx	!x index
	do i = 1, Ny !y index
		U(i, j) = COS(2*pi*((j-1)*dx)) * SIN(2*pi*(Ny-i)*dy)
		V(i, j) =-SIN(2*pi*((j-1)*dx)) * COS(2*pi*(Ny-i)*dy)
	enddo
enddo



! Initial condition
do j = 1, Nx
	do i = 1, Ny
		T(i, j) = complex(sin((Ny-i)*dy), 0.)
	enddo
enddo
!-------------------------------------------------------------------------------------------------------- 
kx1 = 2.0_8*PI / Lx			! First wavenumber
ky1 = 2.0_8*PI / Ly			! First wavenumber

! Define wave numbers using FFTW storage order [ n= 0, n=1,2,3...N/2-1, N/2, -N/2+1,-N/2+2,...,-1]
do k=1,Nx/2+1
		kx(k) = kx1*(k - 1)
end do 
do k=Nx/2+2,Nx
		kx(k) = kx1*(-Nx - 1 + k )	
end do  
do k=1,Ny/2+1
		ky(k) = ky1*(k - 1)
end do 
do k=Ny/2+2,Ny
		ky(k) = ky1*(-Ny - 1 + k )	
end do  

! Create meshgrid of kx and ky

do i=1,Ny
	kxx(i,:) = kx*kx
	kyy(:,i) = ky*ky
end do

! Create meshgrid of complex kx and ky
do j=1,Nx
    do i=1,Ny
		c_kx(i,j) = complex(0.,kx(j))
		c_ky(i,j) = complex(0.,ky(i))
	end do 
end do

! INITIALISE PLANS		  
!see http://www.fftw.org/fftw3_doc/Overview-of-Fortran-interface.html#Overview-of-Fortran-interface  
planA = fftw_plan_dft_2d(Ny, Nx, T, Tf, FFTW_FORWARD, FFTW_ESTIMATE) !Swap rows and columns ( C to fortran memory layout)
planB = fftw_plan_dft_2d(Ny, Nx, Tf, T, FFTW_BACKWARD, FFTW_ESTIMATE) !Swap rows and columns ( C to fortran memory layout)

! INITIALISE DATA
!T(:,:) = COMPLEX(1., 0.)

call fftw_execute_dft(planA, T, Tf)

Tf = Tf 

! Iteration loop
do while(time<tmax)
	time = time + dt
	snapshot = snapshot + 1
	Tf1 = Tf*c_kx
	Tf2 = Tf*c_ky
	call fftw_execute_dft(planB, Tf1, T1)
	
	T_temp1 = real(T1)*U / (Nx * Ny)

	call fftw_execute_dft(planA, T_temp1, T_temp1f)

	call fftw_execute_dft(planB, Tf2, T2)
	
	T_temp2 = real(T2)*V / (Nx * Ny)
	
	call fftw_execute_dft(planA, T_temp2, T_temp2f)
	
	Tf = (-alpha*dt*(kxx+kyy) + 1.0_8) * Tf - T_temp1f*dt- T_temp2f*dt 
	
	!if (mod(snapshot,200)==0) then ! output snapshot every 200 iteration
	!	call fftw_execute_dft(planB, Tf, T)
	!	open(unit=2,file="Temp_assignment2.dat")
	!	write(2,*) 'title="Tecplot"'
	!	write(2,*) 'variables="x" "y" "Tem"'
	!	write(2,*) 'zone t="Fortran"  j=' ,Nx,' i=',Ny,' f=point'
	
	!	do j=1,Nx
	!		do i=1,Ny
	!			write(2,*) j,i,real(T(i,j))
	!		end do
	!	end do
	!endif
end do

call fftw_execute_dft(planB, Tf, T)
T = T/ (Nx * Ny)

call fftw_destroy_plan(planA)
call fftw_destroy_plan(planB)
call fftw_destroy_plan(planC)
call fftw_destroy_plan(planD)
call fftw_destroy_plan(planE)
call fftw_destroy_plan(planF)

open(unit=2,file="Temp_assignment2.dat")
	write(2,*) 'title="Tecplot"'
	write(2,*) 'variables="x" "y" "Tem"'
	write(2,*) 'zone t="Fortran"  j=',Nx,' i=',Ny,' f=point'
do j=1,Nx
    do i=1,Ny
		write(2,*) j,i,real(T(i,j))
    end do
end do
close(2) 
call cpu_time(end)
write(*,*) "T(middle)= ", T(Nx/2, Ny/2) ! Temperature at middle point
write(*,*) "Programme runtime = ", end-start ! Calculation time

end program assignment2

 
