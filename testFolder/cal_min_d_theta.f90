module para
implicit none
integer,parameter::sp=kind(1.0)
integer,parameter::dp=kind(1.0d0)
integer,public,save::Theta_length
real(dp),public,save,allocatable::Theta(:,:)
real(dp),public,save,allocatable::Distance(:,:)
real(dp),public,save,allocatable::Theta_distance(:,:)
end module para
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program cal_angle
use para
implicit none
integer,parameter::length=50
integer,parameter::m=16
integer,parameter::Totn_Theta=800
real(dp),parameter::eps_theta=1.0d-8
real(dp)::theta_change
real(dp)::distance_change
integer :: i
integer :: j
if(.not.allocated(Theta)) allocate(Theta(length,m))
if(.not.allocated(Distance)) allocate(Distance(length,m))
if(.not.allocated(Theta_distance)) allocate(Theta_distance(Totn_Theta,2))
write(*,*)"start reading"
Open(11,file="inputs.txt")
do i=1,50
read(11,*)Theta(i,1),Distance(i,1),Theta(i,2),Distance(i,2),&
         &Theta(i,3),Distance(i,3),Theta(i,4),Distance(i,4),&
         &Theta(i,5),Distance(i,5),Theta(i,6),Distance(i,6),&
         &Theta(i,7),Distance(i,7),Theta(i,8),Distance(i,8),&
         Theta(i,9),Distance(i,9),Theta(i,10),Distance(i,10),&
         &Theta(i,11),Distance(i,11),Theta(i,12),Distance(i,12),&
         &Theta(i,13),Distance(i,13),Theta(i,14),Distance(i,14),&
         &Theta(i,15),Distance(i,15),Theta(i,16),Distance(i,16)
end do
close(11)
do i=1,50
   Theta_distance(i,1)=Theta(i,1)
   Theta_distance(i,2)=Distance(i,1)
end do
do i=51,100
   Theta_distance(i,1)=Theta(i-50,2)
   Theta_distance(i,2)=Distance(i-50,2)
end do
do i=101,150
   Theta_distance(i,1)=Theta(i-100,3)
   Theta_distance(i,2)=Distance(i-100,3)
end do
do i=151,200
   Theta_distance(i,1)=Theta(i-150,4)
   Theta_distance(i,2)=Distance(i-150,4)
end do
do i=201,250
   Theta_distance(i,1)=Theta(i-200,5)
   Theta_distance(i,2)=Distance(i-200,5)
end do
do i=251,300
   Theta_distance(i,1)=Theta(i-250,6)
   Theta_distance(i,2)=Distance(i-250,6)
end do
do i=301,350
   Theta_distance(i,1)=Theta(i-300,7)
   Theta_distance(i,2)=Distance(i-300,7)
end do
do i=351,400
   Theta_distance(i,1)=Theta(i-350,8)
   Theta_distance(i,2)=Distance(i-350,8)
end do
do i=401,450
   Theta_distance(i,1)=Theta(i-400,9)
   Theta_distance(i,2)=Distance(i-400,9)
end do
do i=451,500
   Theta_distance(i,1)=Theta(i-450,10)
   Theta_distance(i,2)=Distance(i-450,10)
end do
do i=501,550
   Theta_distance(i,1)=Theta(i-500,11)
   Theta_distance(i,2)=Distance(i-500,11)
end do
do i=551,600
   Theta_distance(i,1)=Theta(i-550,12)
   Theta_distance(i,2)=Distance(i-550,12)
end do
do i=601,650
   Theta_distance(i,1)=Theta(i-600,13)
   Theta_distance(i,2)=Distance(i-600,13)
end do
do i=651,700
   Theta_distance(i,1)=Theta(i-650,14)
   Theta_distance(i,2)=Distance(i-650,14)
end do
do i=701,750
   Theta_distance(i,1)=Theta(i-700,15)
   Theta_distance(i,2)=Distance(i-700,15)
end do
do i=751,800
   Theta_distance(i,1)=Theta(i-750,16)
   Theta_distance(i,2)=Distance(i-750,16)
end do
do i=1,Totn_Theta
write(*,*)Theta_distance(i,1),Theta_distance(i,2)
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do j=1,Totn_Theta-1
   do i=1,Totn_Theta-1
      if(Theta_distance(i,1)>Theta_distance(i+1,1))then
        theta_change=Theta_distance(i+1,1)
        distance_change=Theta_distance(i+1,2)
        Theta_distance(i+1,1)=Theta_distance(i,1)
        Theta_distance(i+1,2)=Theta_distance(i,2)
        Theta_distance(i,1)=theta_change
        Theta_distance(i,2)=distance_change
      end if
   end do
end do
Open(12,file="arrange.txt")
do i=1,Totn_Theta
   write(12,*)Theta_distance(i,1),Theta_distance(i,2)
end do
close(12)
Open(13,file="useful.txt")
do i=2,Totn_Theta
   if(abs(Theta_distance(i,1)-Theta_distance(i-1,1))>eps_theta)then
     write(13,*)Theta_distance(i-1,1),Theta_distance(i-1,2)
   else if(abs(Theta_distance(i,1)-Theta_distance(i-1,1))<=eps_theta.and.Theta_distance(i,2)>Theta_distance(i-1,2))then
        theta_change=Theta_distance(i-1,1)
        distance_change=Theta_distance(i-1,2)
        Theta_distance(i,1)=Theta_distance(i-1,1)
        Theta_distance(i,2)=Theta_distance(i-1,2)
        Theta_distance(i-1,1)=theta_change
        Theta_distance(i-1,2)=distance_change    
   end if
end do
if(abs(Theta_distance(Totn_Theta,1)-Theta_distance(Totn_Theta-1,1))>eps_theta)then
write(13,*)Theta_distance(Totn_Theta,1),Theta_distance(Totn_Theta,2)
end if
close(13)
end