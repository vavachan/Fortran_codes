module function_set
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
implicit none 
integer :: nofatoms,time,flag2,coun
real , dimension (3) :: box
integer ::i,j,k,p,IN=10,IN2=11,IN3=12,IN4=13,t
integer , dimension(:,:), allocatable :: break
real , dimension (:,:) , allocatable :: pos,initial
real , dimension (:,:) , allocatable :: vel
real , dimension (:,:) , allocatable :: acc
real , dimension (3) :: temp
real rr,alpha,temperature,kin,poten,temp1,temp2 
real , parameter :: delta = .0005 !Modify this to change step size(time)
real :: R(3)
real :: ivel,Rcut=2.5,randgauss!Rcut is the cut-off radius for LJ potential
real , dimension (:) , allocatable :: pot,ti,mea,diff
real :: velocity,phi,dphi,msd2
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine poscheck
!This subroutine implements the boundary conditions
do i=1,nofatoms
		
	do	
		flag2=0
		do k=1,3
			if (pos(i,k)>=box(k)/2) then
				pos(i,k)=pos(i,k)-box(k)
				break(i,k)=break(i,k)+1
			else if (pos(i,k)<-box(k)/2)then
				pos(i,k)=pos(i,k)+box(k)
				break(i,k)=break(i,k)-1
			else 
				flag2=flag2+1
			end if 
	
		end do 
		if (flag2==3) exit
	
	end do 
end do 
end subroutine poscheck
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine gauss
real :: v1,v2,rsq
integer, dimension(:), allocatable :: seed
integer ::m,clock
call random_seed (size=m)
allocate (seed(m))
do i=1,m
	call system_clock (Count=clock)
	seed(i)=clock
end do
call random_seed (PUT=seed) 
do 
	call random_number(v1)
	v1=2*v1-1
	call random_number(v2)
	v2=2*v2-1
	rsq=v1*v1+v2*v2

	if (rsq .LE. 1.0 .AND. rsq .GE. 0.0) then
	
	exit
	end if 
end do 
randgauss=v1*SQRT(-2.0*LOG(rsq)/rsq)

end subroutine gauss 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine inivel(vel)
!This subroutine initialzes the velocity vector, the velocities are randomly selected and its magnitudes are scaled according to the temperature.
real , dimension (:,:) , allocatable :: vel
integer :: i
real :: rt(3),sum1
do i=1,nofatoms
	sum1=0
	do k=1,3
		call gauss
		rt(k)=randgauss
		sum1=sum1+rt(k)**2
	end do
	sum1=sum1**(.5)
	 
	do k=1,3
	rt(k)=rt(k)/sum1
	vel(i,k)=ivel*rt(k)
	end do 
	
end do
end subroutine inivel
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine  inipos(pos) 
!Subroutine to initialize position. Particles are randomly placed in a 3D lattice with so that density is uniform. 
!
real , dimension (:,:) ,allocatable :: pos
integer :: k,i,j,nb(3),flag=1,IN=10,OUT=11,clock,m,meshsize
real :: density,volume,try,r(3)
real :: step
real :: gap(3)
integer, dimension(:), allocatable :: seed
real , dimension (:,:,:) , allocatable :: mesh 
call random_seed (size=m)
allocate (seed(m))
allocate (initial(nofatoms,3))
do i=1,m
	call system_clock (Count=clock)
	seed(i)=clock
end do

call random_seed (PUT=seed) 

volume =1
try=1
do k=1,3
	volume=volume*box(k)
end do
do i=1,nofatoms

	if(i**3>=nofatoms) then 
		meshsize=i
		exit 
	end if 

end do  
density = nofatoms/volume
write (*,*) 'density' , density 
allocate(mesh(meshsize,meshsize,meshsize))
do i=1,meshsize
	do j=1,meshsize
		do k=1,meshsize
			mesh(i,j,k)=0
		end do
	end do
end do 
do i=1,nofatoms
	do 
		do k=1,3
			call random_number(try)
			nb(k)=try*(meshsize)
!			pos(i,k)=nb(k)*box(k)/meshsize
		end do 
		if (mesh(nb(1),nb(2),nb(3))==0) then
			mesh(nb(1),nb(2),nb(3))=1
			do k=1,3
				pos(i,k)=nb(k)*box(k)/meshsize
			end do
			exit
		end if
	end do 
end do 
!The matrix inital has the data of inipital position of particles , which will be needed when we calculate Mean Square Displacement.
do i=1,nofatoms
	do k=1,3
		initial(i,k)=pos(i,k)
		break(i,k)=0
	end do 
end do 
end subroutine inipos
!Subroutines to display velocities, position and acceleration of particles.
!######################################################################
subroutine disvel
do i=1,nofatoms
	write (*,*) 'Particle',i,'at(', vel(i,1),',',vel(i,2),',',vel(i,3),')',' ',vel(i,1)**2+vel(i,2)**2+vel(i,3)**2
end do
end subroutine disvel
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1!
subroutine dispos
do i=1,nofatoms
	write (*,*) 'Particle',i,'at(', pos(i,1),',',pos(i,2),',',pos(i,3),')'
end do
end subroutine dispos 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine disacc
do i=1,nofatoms
	write (*,*) 'Particle',i,'at(', acc(i,1),',',acc(i,2),',',acc(i,3),')'
end do
end subroutine disacc
!#########################################################################
!Subroutine to diplay everything
subroutine displayall
do i=1,nofatoms
	write (*,*) 'position',i,'at(', pos(i,1),',',pos(i,2),',',pos(i,3),')'
	write (*,*) 'accelera',i,'at(', acc(i,1),',',acc(i,2),',',acc(i,3),')'
	write (*,*) 'velocity',i,'at(', vel(i,1),',',vel(i,2),',',vel(i,3),')'
end do
end subroutine displayall
!#################################################################################
!# 	Subroutines to calulate force.					         
!#In this Simulation we are using Lennard-Jones Potential. 
!#################################################################################
!Subroutine to find certatin constants needed for calculating the force.The input 
!variable is the distance between atoms
!#################################################################################
subroutine calc(r)
real :: r 
phi=4*(1/r**12-1/r**6)
dphi=-48/r*(1/r**12-0.5/r**6)

end subroutine calc
!To Find the Force.
subroutine force
do i=1,nofatoms
		do k=1,3
			acc(i,k)=0
		end do 
end do
poten=0
do i=1,nofatoms-1
	do j=i+1,nofatoms!We don't need to consider all pairs because of newtons second law,ie force on first particle by second =(-1) X force on second by first
		rr=0
		do k=1,3
			R(k)=pos(i,k)-pos(j,k)
		
			if(R(k)>=box(k)/2) then 
				R(k)=R(k)-box(k)!Since we are using periodic boundary condition while calculating the distances we consider the box around i'th particle
			else if (R(k)< - box(k)/2) then
				R(k)=R(k)+box(k)
			end if 
				rr=rr+R(k)**2
		end do
		rr=rr**(0.5)
		if (rr<Rcut) then !implementing the cut-off radius.
			call calc(Rcut)
			temp1=phi
			temp2=dphi
			call calc(rr)
			alpha=(-1/rr)*(dphi)
			do k=1,3
				acc(i,k)=acc(i,k)+R(k)*(alpha+(1/rr)*temp2) !Updating the acceleration 
				acc(j,k)=acc(j,k)-R(k)*(alpha+(1/rr)*temp2) !Newton's Third Law!
			end do 
			poten=poten+(phi-temp1-(rr-Rcut)*temp2)!Poten is potential energy , given the adjustments because we are using a cutoff.
		end if
		
	end do
end do 
end subroutine force 
!Velocity update (We are using velocity-verlet method).
subroutine velupdate
do i=1,nofatoms
	do k=1,3
		vel(i,k)=vel(i,k)+delta/2*acc(i,k)
	end do
end do 
end subroutine velupdate

subroutine posupdate

do i=1,nofatoms
	do k=1,3
		pos(i,k)=pos(i,k)+vel(i,k)*delta
	end do
end do

end subroutine posupdate
!#########################################################
!Subroutne to calculate mean-squre displacement
subroutine means
real :: msd1
msd2=0
do i=1,nofatoms
	msd1=0
	do k=1,3
		msd1=msd1+(pos(i,k)-initial(i,k)+break(i,k)*box(k))**2
	end do
	msd2=msd2+msd1
end do 
end subroutine means 
!############################################################
subroutine calckin
!Calculating kinetic energy
do i=1,nofatoms
	velocity=0
	do k=1,3
		velocity=velocity+vel(i,k)*vel(i,k)
	end do 
	pot(i)=velocity**(0.5)
	kin=kin+velocity
end do 
end subroutine calckin
!###########################################################
end module function_set
!###########################################################
!###########################################################
!BEGIN MAIN PROGRAM
program moledyn
use function_set
implicit none 
real:: potcut
!initialisations
write (*,*) 'NO of atoms'
read*, nofatoms
write (*,*) 'Specify the Dimensions of the Box'
t=0
do i=1,3
read*, box(i)
end do 
write (*,*) 'Enter the Temperature'
read*, temperature
write (*,*) 'Enter time'
read*, time
ivel=(3*temperature)**(0.5)

allocate (pos(nofatoms,3))
allocate (vel(nofatoms,3))
allocate (acc(nofatoms,3))
allocate (break(nofatoms,3))
allocate (ti(time))
allocate (mea(time))
allocate (diff(time))

!initial values
call inipos(pos)
call inivel(vel)
open (unit=IN,file='position.dat',ACTION='READWRITE')
open (unit=IN2,file='energy.dat',ACTION='READWRITE')
open (unit=IN4,file='MSD.dat', ACTION='READWRITE')
allocate (pot(nofatoms))
!MAIN LOOP
do t=1,time
	kin=0
	poten=0
	call force
	call velupdate
	call posupdate		
	call poscheck
	call force 
	call velupdate
	call calckin
	call means 
	write (IN2,*) t,kin/(2*nofatoms),poten/nofatoms,(kin/(2)+poten)/nofatoms,kin/(3*nofatoms)	
	write (IN4,*) t,msd2/(nofatoms)
end do
end program moledyn 
