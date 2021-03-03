module para
implicit none
integer,parameter::sp=kind(1.0)
integer,parameter::dp=kind(1.0d0)
real(dp),public,save::pi
real(dp),public,save::a1(2)
real(dp),public,save::a2(2)
real(dp),public,save,allocatable::frac_coor_L1(:,:,:)
real(dp),public,save,allocatable::frac_coor_L2(:,:,:)
real(dp),public,save,allocatable::xoy_frac_coor_L1(:,:,:)
real(dp),public,save,allocatable::xoy_frac_coor_L2(:,:,:)
real(dp),public,save,allocatable::Rxoy_frac_coor_L1(:,:,:)
real(dp),public,save,allocatable::Rxoy_frac_coor_L2(:,:,:)
real(dp),public,save,allocatable::Rfrac_coor_L1(:,:,:)
real(dp),public,save,allocatable::Rfrac_coor_L2(:,:,:)
real(dp),public,save,allocatable::coordinate_L1(:,:,:)
real(dp),public,save,allocatable::coordinate_L2(:,:,:)
end module para
program cal_coordinate
use para
implicit none
integer::Ntranslation_a1
integer::Ntranslation_a2
integer::n_cell
integer::i
integer::j
integer::j_atom
real(dp)::rotational_theta
real(dp)::superlattice_size
real(dp)::superlattice_vector1(2)
real(dp)::superlattice_vector2(2)
real(dp)::atom1(2)
real(dp)::atom2(2)
real(dp)::Rotation_matrix(2,2)
real(dp)::INV_Rotation_matrix(2,2)
complex(dp)::CRotation_matrix(2,2)
complex(dp)::CINV_Rotation_matrix(2,2)
pi=dasin(1.0d0)*2.0d0
a1(1)=2.46d0
a1(2)=0.0d0
a2(1)=-2.46d0*sin(pi/6.0d0)
a2(2)=2.46d0*cos(pi/6.0d0)
Ntranslation_a1=20
Ntranslation_a2=20
n_cell=0
do i=-Ntranslation_a1,Ntranslation_a1
   do j=-Ntranslation_a2,Ntranslation_a2
      n_cell=n_cell+1
   end do
end do
write(*,*)n_cell
atom1(1)=0.0d0
atom1(2)=0.0d0
atom2(1)=1.0d0/3.0d0
atom2(2)=2.0d0/3.0d0
if(.not.allocated(frac_coor_L1)) allocate(frac_coor_L1(n_cell,2,2))
if(.not.allocated(frac_coor_L2)) allocate(frac_coor_L2(n_cell,2,2))
if(.not.allocated(xoy_frac_coor_L1)) allocate(xoy_frac_coor_L1(n_cell,2,2))
if(.not.allocated(xoy_frac_coor_L2)) allocate(xoy_frac_coor_L2(n_cell,2,2))
if(.not.allocated(Rxoy_frac_coor_L1)) allocate(Rxoy_frac_coor_L1(n_cell,2,2))
if(.not.allocated(Rxoy_frac_coor_L2)) allocate(Rxoy_frac_coor_L2(n_cell,2,2))
if(.not.allocated(Rfrac_coor_L1)) allocate(Rfrac_coor_L1(n_cell,2,2))
if(.not.allocated(Rfrac_coor_L2)) allocate(Rfrac_coor_L2(n_cell,2,2))
if(.not.allocated(coordinate_L1)) allocate(coordinate_L1(n_cell,2,2))
if(.not.allocated(coordinate_L2)) allocate(coordinate_L2(n_cell,2,2))

do i=1,n_cell
   do j=1,2
      do j_atom=1,2
         frac_coor_L1(i,j,j_atom)=0.0d0
         frac_coor_L2(i,j,j_atom)=0.0d0
         xoy_frac_coor_L1(i,j,j_atom)=0.0d0
         xoy_frac_coor_L2(i,j,j_atom)=0.0d0
         Rxoy_frac_coor_L1(i,j,j_atom)=0.0d0
         Rxoy_frac_coor_L2(i,j,j_atom)=0.0d0
         Rfrac_coor_L1(i,j,j_atom)=0.0d0
         Rfrac_coor_L2(i,j,j_atom)=0.0d0
         coordinate_L1(i,j,j_atom)=0.0d0
         coordinate_L2(i,j,j_atom)=0.0d0
      end do
   end do
end do
n_cell=0
do i=-Ntranslation_a1,Ntranslation_a1
   do j=-Ntranslation_a2,Ntranslation_a2
      n_cell=n_cell+1
      frac_coor_L1(n_cell,1,1)=atom1(1)-dble(i)
      frac_coor_L1(n_cell,2,1)=atom1(2)-dble(j)
      frac_coor_L1(n_cell,1,2)=atom2(1)-dble(i)
      frac_coor_L1(n_cell,2,2)=atom2(2)-dble(j)
      frac_coor_L2(n_cell,1,1)=atom1(1)-dble(i)
      frac_coor_L2(n_cell,2,1)=atom1(2)-dble(j)
      frac_coor_L2(n_cell,1,2)=atom2(1)-dble(i)
      frac_coor_L2(n_cell,2,2)=atom2(2)-dble(j)
   end do
end do
!!!!!!!!!!!!!--Define the superlattice_size--!!!!!!!!!!!!!!
superlattice_size=SQRT(39.0d0)*2.46d0/SQRT(3.0d0)
!!!!!!!!!!!!!--Define the rotational_theta--!!!!!!!!!!!!!!!
rotational_theta=2.0d0*dasin(3.0d0/(2.0d0*SQRT(39.0d0)))
write(*,*)"rotational angle:",rotational_theta*180/pi
write(*,*)"superlattice size:",superlattice_size
!!!!!!!!!!!!!!!!!!!!--unit cell vector--!!!!!!!!!!!!!!!!!!!
superlattice_vector1(1)=a1(1)*superlattice_size/2.46d0
superlattice_vector1(2)=a1(2)*superlattice_size/2.46d0
superlattice_vector2(1)=a2(1)*superlattice_size/2.46d0
superlattice_vector2(2)=a2(2)*superlattice_size/2.46d0
write(*,*)"unrotated superlattice vector-clockwise:"
write(*,*)superlattice_vector1(1),superlattice_vector1(2)
write(*,*)superlattice_vector2(1),superlattice_vector2(2)
!!!!!!!!!!!!!!!!!--fractional coordinate--!!!!!!!!!!!!!!!!!
do i=1,n_cell
   frac_coor_L1(i,1,1)=frac_coor_L1(i,1,1)*2.46d0/superlattice_size
   frac_coor_L1(i,2,1)=frac_coor_L1(i,2,1)*2.46d0/superlattice_size
   frac_coor_L1(i,1,2)=frac_coor_L1(i,1,2)*2.46d0/superlattice_size
   frac_coor_L1(i,2,2)=frac_coor_L1(i,2,2)*2.46d0/superlattice_size
   frac_coor_L2(i,1,1)=frac_coor_L2(i,1,1)*2.46d0/superlattice_size
   frac_coor_L2(i,2,1)=frac_coor_L2(i,2,1)*2.46d0/superlattice_size
   frac_coor_L2(i,1,2)=frac_coor_L2(i,1,2)*2.46d0/superlattice_size
   frac_coor_L2(i,2,2)=frac_coor_L2(i,2,2)*2.46d0/superlattice_size
end do
!!!!!!!!!!!!!!!--rotational unit cell vector--!!!!!!!!!!!!!
write(*,*)"========================================================"
write(*,*)"========================================================"
write(*,*)"===========---for anticlockwise rotation---============="
do i=1,2
   do j=1,2
      Rotation_matrix(i,j)=0.0d0
      INV_Rotation_matrix(i,j)=0.0d0
      CRotation_matrix(i,j)=dcmplx(0.0_dp,0.0_dp)
      CINV_Rotation_matrix(i,j)=dcmplx(0.0_dp,0.0_dp)     
   end do
end do
write(*,*)"rotational angle:",-rotational_theta/2.0d0*180/pi
Rotation_matrix(1,1)=1.0d0
Rotation_matrix(1,2)=-dsin(pi/6.0d0)
Rotation_matrix(2,1)=0.0d0
Rotation_matrix(2,2)=dcos(pi/6.0d0)
!!!!!!!!!!!--rotate fractional coordinates--!!!!!!!!!!!!!!
do i=1,2
   do j=1,2
      CRotation_matrix(i,j)=CRotation_matrix(i,j)+Rotation_matrix(i,j)
   end do
end do
call INVERT(2,CRotation_matrix,CINV_Rotation_matrix)
do i=1,2
   do j=1,2
      INV_Rotation_matrix(i,j)=Real(CINV_Rotation_matrix(i,j))
   end do
end do
write(*,*)"rotate matrix:"
write(*,*)Rotation_matrix(1,1),Rotation_matrix(1,2)
write(*,*)Rotation_matrix(2,1),Rotation_matrix(2,2)
write(*,*)"Invert rotate matrix:"
write(*,*)INV_Rotation_matrix(1,1),INV_Rotation_matrix(1,2)
write(*,*)INV_Rotation_matrix(2,1),INV_Rotation_matrix(2,2)
do i=1,n_cell
   do j_atom=1,2
      xoy_frac_coor_L1(i,1,j_atom)=Rotation_matrix(1,1)*frac_coor_L1(i,1,j_atom)+Rotation_matrix(1,2)*frac_coor_L1(i,2,j_atom)
      xoy_frac_coor_L1(i,2,j_atom)=Rotation_matrix(2,1)*frac_coor_L1(i,1,j_atom)+Rotation_matrix(2,2)*frac_coor_L1(i,2,j_atom)
   end do
end do
do i=1,n_cell
   do j_atom=1,2
      Rxoy_frac_coor_L1(i,1,j_atom)=xoy_frac_coor_L1(i,1,j_atom)*dcos(-rotational_theta/2.0d0)+&
                                   &+xoy_frac_coor_L1(i,2,j_atom)*dsin(-rotational_theta/2.0d0)
      Rxoy_frac_coor_L1(i,2,j_atom)=-xoy_frac_coor_L1(i,1,j_atom)*dsin(-rotational_theta/2.0d0)+&
                                   &+xoy_frac_coor_L1(i,2,j_atom)*dcos(-rotational_theta/2.0d0)
   end do
end do
do i=1,n_cell
   do j_atom=1,2
      Rfrac_coor_L1(i,1,j_atom)=INV_Rotation_matrix(1,1)*Rxoy_frac_coor_L1(i,1,j_atom)+&
                               &INV_Rotation_matrix(1,2)*Rxoy_frac_coor_L1(i,2,j_atom)
      Rfrac_coor_L1(i,2,j_atom)=INV_Rotation_matrix(2,1)*Rxoy_frac_coor_L1(i,1,j_atom)+&
                               &INV_Rotation_matrix(2,2)*Rxoy_frac_coor_L1(i,2,j_atom)
   end do
end do
write(*,*)"========================================================"
write(*,*)"========================================================"
write(*,*)"=============---for clockwise rotation---==============="
do i=1,2
   do j=1,2
      Rotation_matrix(i,j)=0.0d0
      INV_Rotation_matrix(i,j)=0.0d0
      CRotation_matrix(i,j)=dcmplx(0.0_dp,0.0_dp)
      CINV_Rotation_matrix(i,j)=dcmplx(0.0_dp,0.0_dp)     
   end do
end do
write(*,*)"rotational angle:",rotational_theta/2.0d0*180/pi
Rotation_matrix(1,1)=1.0d0
Rotation_matrix(1,2)=-dsin(pi/6.0d0)
Rotation_matrix(2,1)=0.0d0
Rotation_matrix(2,2)=dcos(pi/6.0d0)
!!!!!!!!!!!--rotate fractional coordinates--!!!!!!!!!!!!!!
do i=1,2
   do j=1,2
      CRotation_matrix(i,j)=CRotation_matrix(i,j)+Rotation_matrix(i,j)
   end do
end do
call INVERT(2,CRotation_matrix,CINV_Rotation_matrix)
do i=1,2
   do j=1,2
      INV_Rotation_matrix(i,j)=Real(CINV_Rotation_matrix(i,j))
   end do
end do
write(*,*)"rotate matrix:"
write(*,*)Rotation_matrix(1,1),Rotation_matrix(1,2)
write(*,*)Rotation_matrix(2,1),Rotation_matrix(2,2)
write(*,*)"Invert rotate matrix:"
write(*,*)INV_Rotation_matrix(1,1),INV_Rotation_matrix(1,2)
write(*,*)INV_Rotation_matrix(2,1),INV_Rotation_matrix(2,2)
do i=1,n_cell
   do j_atom=1,2
      xoy_frac_coor_L2(i,1,j_atom)=Rotation_matrix(1,1)*frac_coor_L2(i,1,j_atom)+&
                                  &Rotation_matrix(1,2)*frac_coor_L2(i,2,j_atom)
      xoy_frac_coor_L2(i,2,j_atom)=Rotation_matrix(2,1)*frac_coor_L2(i,1,j_atom)+&
                                  &Rotation_matrix(2,2)*frac_coor_L2(i,2,j_atom)
   end do
end do
do i=1,n_cell
   do j_atom=1,2
Rxoy_frac_coor_L2(i,1,j_atom)=xoy_frac_coor_L2(i,1,j_atom)*dcos(rotational_theta/2.0d0)+&
                             &+xoy_frac_coor_L2(i,2,j_atom)*dsin(rotational_theta/2.0d0)
Rxoy_frac_coor_L2(i,2,j_atom)=-xoy_frac_coor_L2(i,1,j_atom)*dsin(rotational_theta/2.0d0)+&
                             &+xoy_frac_coor_L2(i,2,j_atom)*dcos(rotational_theta/2.0d0)
   end do
end do
do i=1,n_cell
   do j_atom=1,2
      Rfrac_coor_L2(i,1,j_atom)=INV_Rotation_matrix(1,1)*Rxoy_frac_coor_L2(i,1,j_atom)+&
                               &INV_Rotation_matrix(1,2)*Rxoy_frac_coor_L2(i,2,j_atom)
      Rfrac_coor_L2(i,2,j_atom)=INV_Rotation_matrix(2,1)*Rxoy_frac_coor_L2(i,1,j_atom)+&
                               &INV_Rotation_matrix(2,2)*Rxoy_frac_coor_L2(i,2,j_atom)
   end do
end do
Open(11,file="useful.txt")
do i=1,n_cell
   do j_atom=1,2
      if(Rfrac_coor_L1(i,1,j_atom)>=0.0d0.and.Rfrac_coor_L1(i,1,j_atom)<1.0d0.and.&
         &Rfrac_coor_L1(i,2,j_atom)>=0.0d0.and.Rfrac_coor_L1(i,2,j_atom)<1.0d0)then
         coordinate_L1(i,1,j_atom)=Rfrac_coor_L1(i,1,j_atom)*superlattice_vector1(1)+&
                                  &Rfrac_coor_L1(i,2,j_atom)*superlattice_vector2(1)
         coordinate_L1(i,2,j_atom)=Rfrac_coor_L1(i,1,j_atom)*superlattice_vector1(2)+&
                                  &Rfrac_coor_L1(i,2,j_atom)*superlattice_vector2(2)
         write(11,"(2f14.10)")coordinate_L1(i,1,j_atom),coordinate_L1(i,2,j_atom)
      end if
   end do
end do
         write(11,"(2f14.10)")
do i=1,n_cell
   do j_atom=1,2
      if(Rfrac_coor_L2(i,1,j_atom)>=0.0d0.and.Rfrac_coor_L2(i,1,j_atom)<1.0d0.and.&
         &Rfrac_coor_L2(i,2,j_atom)>=0.0d0.and.Rfrac_coor_L2(i,2,j_atom)<1.0d0)then
         coordinate_L2(i,1,j_atom)=Rfrac_coor_L2(i,1,j_atom)*superlattice_vector1(1)+&
                                  &Rfrac_coor_L2(i,2,j_atom)*superlattice_vector2(1)
         coordinate_L2(i,2,j_atom)=Rfrac_coor_L2(i,1,j_atom)*superlattice_vector1(2)+&
                                  &Rfrac_coor_L2(i,2,j_atom)*superlattice_vector2(2)
         write(11,"(2f14.10)")coordinate_L2(i,1,j_atom),coordinate_L2(i,2,j_atom)
      end if
   end do
end do
close(11)
Open(12,file="fractional_coor.txt")
do i=1,n_cell
   do j_atom=1,2
      if(Rfrac_coor_L1(i,1,j_atom)>=0.0d0.and.Rfrac_coor_L1(i,1,j_atom)<1.0d0.and.&
         &Rfrac_coor_L1(i,2,j_atom)>=0.0d0.and.Rfrac_coor_L1(i,2,j_atom)<1.0d0)then
         write(12,"(2f14.10)")Rfrac_coor_L1(i,1,j_atom),Rfrac_coor_L1(i,2,j_atom)
      end if
   end do
end do
         write(12,"(2f14.10)")
do i=1,n_cell
   do j_atom=1,2
      if(Rfrac_coor_L2(i,1,j_atom)>=0.0d0.and.Rfrac_coor_L2(i,1,j_atom)<1.0d0.and.&
         &Rfrac_coor_L2(i,2,j_atom)>=0.0d0.and.Rfrac_coor_L2(i,2,j_atom)<1.0d0)then
         write(12,"(2f14.10)")Rfrac_coor_L2(i,1,j_atom),Rfrac_coor_L2(i,2,j_atom)
      end if
   end do
end do
close(12)
end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>>> to calculat the inverse of the matrix for complex matrix !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     SUBROUTINE INVERT(N,A,C)
	 IMPLICIT DOUBLE PRECISION(A-H,O-Z)
     INTEGER     LWORK, INFO, LDA, N, M
	 COMPLEX*16  A(N,N),B(N,N),C(N,N),D(N,N)
	 INTEGER     IPIV(N)
	 PARAMETER   (LWRKMAX=1000)
	 COMPLEX*16,ALLOCATABLE:: WORK(:)
	 CHARACTER   OPT
!*
!*  -- LAPACK routine (version 3.2) --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     November 2006
!*
!*     .. Scalar Arguments ..
!*      INTEGER            INFO, LDA, LWORK, N
!*     ..
!*     .. Array Arguments ..
!*      INTEGER            IPIV( * )
!*      COMPLEX*16         A( LDA, * ), WORK( * )
!*     ..
!*
!*  Purpose
!*  =======
!*
!*  ZGETRI computes the inverse of a matrix using the LU factorization
!*  computed by ZGETRF.
!*
!*  This method inverts U and then computes inv(A) by solving the system
!*  inv(A)*L = inv(U) for inv(A).
!*
!*  Arguments
!*  =========
!*
!*  N       (input) INTEGER
!*          The order of the matrix A.  N >= 0.
!*
!*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
!*          On entry, the factors L and U from the factorization
!*          A = P*L*U as computed by ZGETRF.
!*          On exit, if INFO = 0, the inverse of the original matrix A.
!*
!*  LDA     (input) INTEGER
!*          The leading dimension of the array A.  LDA >= max(1,N).
!*
!*  IPIV    (input) INTEGER array, dimension (N)
!*          The pivot indices from ZGETRF; for 1<=i<=N, row i of the
!*          matrix was interchanged with row IPIV(i).
!*
!*  WORK    (workspace/output) COMPLEX*16 array, dimension (MAX(1,LWORK))
!*          On exit, if INFO=0, then WORK(1) returns the optimal LWORK.
!*
!*  LWORK   (input) INTEGER
!*          The dimension of the array WORK.  LWORK >= max(1,N).
!*          For optimal performance LWORK >= N*NB, where NB is
!*          the optimal blocksize returned by ILAENV.
!*
!*          If LWORK = -1, then a workspace query is assumed; the routine
!*          only calculates the optimal size of the WORK array, returns
!*          this value as the first entry of the WORK array, and no error
!*          message related to LWORK is issued by XERBLA.
!*
!*  INFO    (output) INTEGER
!*          = 0:  successful exit
!*          < 0:  if INFO = -i, the i-th argument had an illegal value
!*          > 0:  if INFO = i, U(i,i) is exactly zero; the matrix is
!*                singular and its inverse could not be computed.
!*
!* OPT      (INPUT)CHARACTER
!*          on entry,"YES" correspond to LWORK=-1 ,otherwise,LWORK=N    
!*  =====================================================================
!**********************************
     IF(N.EQ.1)THEN
	   C=DCMPLX(1.0D0,0.0D0)/A
	   RETURN
	 ENDIF
!**********************************
     DO I=1,N
	    DO J=1,N
	       B(I,J)=A(I,J)
	       D(I,J)=A(I,J)
	       C(I,J)=DCMPLX(0.0D0,0.0D0)
	    END DO
	 END DO

!C INITIALIZE THESE PARAMETERS
	 LDA=N;LWORK=N;M=N

	 OPT="NO"

	 CALL ZGETRF( M, N, A, LDA, IPIV, INFO )

      IF(OPT.EQ."YES")THEN
	  IF(INFO.EQ.0)THEN
          LWORK=-1
	    IF(.NOT.ALLOCATED(WORK))ALLOCATE(WORK(LWRKMAX))
	    CALL ZGETRI( N, B, LDA, IPIV, WORK, LWORK, INFO )
	    LWORK=MIN(INT(WORK(1)),LWRKMAX)
	    DEALLOCATE(WORK)
	    IF(.NOT.ALLOCATED(WORK))ALLOCATE(WORK(LWORK))

	    CALL ZGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
	    IF(INFO.NE.0)WRITE(*,*)"ERROR IN ZGETRI!"
	   ELSE
	            WRITE(*,*)"INFO=",INFO
	            STOP 'MATRIX INVERSE FAILED'
	   ENDIF
	 ELSE

	  IF(INFO.EQ.0)THEN
	  IF(.NOT.ALLOCATED(WORK))ALLOCATE(WORK(LWORK))
	     CALL ZGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
	    IF(INFO.NE.0) STOP "ERROR IN ZGETRI!"
	  DEALLOCATE(WORK)

	  ELSE
	            WRITE(*,*)"INFO=",INFO
	            STOP 'MATRIX INVERSE FAILED'
	  ENDIF
	 ENDIF
!C TO RESTORE THESE MATRICES

       DO I=1,N
	    DO J=1,N
	       C(I,J)=A(I,J)
	       A(I,J)=D(I,J)
	    END DO
	 END DO

50   CONTINUE
	 RETURN
	 END SUBROUTINE INVERT