program cal_poscar
implicit none
integer,parameter::dp=kind(1.0d0)
integer,parameter::Nxyz=3
integer,parameter::Natom=28
real(dp),parameter::d=2.7d0
real(dp),parameter::superlattice_size=6.5085482243d0
integer::i
integer::j
real(dp)::L1
real(dp)::L2
real(dp)::pi
real(dp)::basis(Nxyz,Nxyz)
real(dp)::fractional_coor(Natom,Nxyz)
Open(12,file="fractional_coor.txt")
    do i=1,Natom
       read(12,"(2f14.10)")fractional_coor(i,1),fractional_coor(i,2)
    end do
close(12)
pi=dasin(1.0d0)*2.0d0
!!!!!!!!!!!!--edit lattice size--!!!!!!!!!!!!!
basis(1,1)=superlattice_size
basis(1,2)=0.0d0
basis(1,3)=0.0d0
basis(2,1)=-superlattice_size*dsin(pi/6.0d0)
basis(2,2)=superlattice_size*dcos(pi/6.0d0)
basis(2,3)=0.0d0
basis(3,1)=0.0d0
basis(3,2)=0.0d0
basis(3,3)=20.0d0+d
L1=0.5d0-d/2.0d0/(20.0d0+d)
L2=0.5d0+d/2.0d0/(20.0d0+d)
Open(13,file="POSCAR")
write(13,"(A)")"Moire_BLG"
write(13,"(f3.1)")1.0
write(13,"(3f20.10)")basis(1,1),basis(1,2),basis(1,3)
write(13,"(3f20.10)")basis(2,1),basis(2,2),basis(2,3)
write(13,"(3f20.10)")basis(3,1),basis(3,2),basis(3,3)
write(13,"(A)")"C"
write(13,"(I5)")Natom
write(13,"(A)")"Direct"
    do i=1,Natom/2
       write(13,"(3f16.10)")fractional_coor(i,1),fractional_coor(i,2),L1
    end do
    do i=Natom/2+1,Natom
       write(13,"(3f16.10)")fractional_coor(i,1),fractional_coor(i,2),L2
    end do
close(13)
end