module constant
implicit none
!-------------------------------integer constants------------------------------!
integer,public,parameter::sp=kind(1.0)
integer,public,parameter::dp=kind(1.0d0)
!------------------------------circumference ratio-----------------------------!
real(dp),public,parameter::pi=3.1415926535897932384626_dp
!--------------------------------real constants--------------------------------!
real(dp),public,parameter::th=1.0_dp
real(dp),public,parameter::zero=0.0_dp
real(dp),public,parameter::one =1.0_dp
!------------------------------complex constants-------------------------------!
complex(dp),public,parameter::czero =dcmplx(0.0_dp,0.0_dp)
complex(dp),public,parameter::cone  =dcmplx(1.0_dp,0.0_dp)
complex(dp),public,parameter::ione  =dcmplx(0.0_dp,1.0_dp)
end module constant
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module para
use constant
implicit none
integer,public,parameter::itermax=1000
integer,public,parameter::Nsite=28
integer,public,parameter::Nalpha=2
real(dp),public,parameter::max_eps_SE=1.0d-4
real(dp),public,parameter::max_eps_n=1.0d-5
real(dp),public,parameter::eta=1.0d-2
integer,public,save::Nw
integer,public,save::NK
integer,public,save::iter
integer,public,save::ifermi
real(dp),public,save::Maxerror_SE
real(dp),public,save::Maxerror_ne
real(dp),public,save::U
real(dp),public,save::Delta1
real(dp),public,save::Delta2
real(dp),public,save::Delta3
real(dp),public,save::Ef
real(dp),public,save::E_max
real(dp),public,save::E_min
real(dp),public,save::domega
real(dp),public,save::Tn0
real(dp),public,save::Tnt
real(dp),public,save,allocatable::error_ne(:)
real(dp),public,save,allocatable::error_SE(:,:)
real(dp),public,save,allocatable::DELTA(:)
real(dp),public,save,allocatable::ne0(:)
real(dp),public,save,allocatable::net(:)
real(dp),public,save,allocatable::DOS(:,:)
real(dp),public,save,allocatable::TDOS(:)
real(dp),public,save,allocatable::P_alpha(:,:)
complex(dp),public,save,allocatable::Omega(:)
complex(dp),public,save,allocatable::SE_old(:,:)
complex(dp),public,save,allocatable::SE_new(:,:)
complex(dp),public,save,allocatable::Cavity_G(:,:)
complex(dp),public,save,allocatable::Lattice_G(:,:)
complex(dp),public,save,allocatable::Impurity_G(:,:)
complex(dp),public,save,allocatable::Impurity_G_alpha(:,:,:)
complex(dp),public,save,allocatable::HK(:,:,:)
complex(dp),public,save,allocatable::HK_k(:,:,:)
complex(dp),public,save,allocatable::INVHK_k(:,:,:)
end module para
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program TBLG_CPA
use constant
use para
implicit none
integer i
integer j
integer k
integer i_site
integer j_site
integer i_alpha
! integer t3, t4
character(20)ch
if(.not.allocated(ne0)) allocate(ne0(Nsite))
if(.not.allocated(net)) allocate(net(Nsite))
if(.not.allocated(DELTA)) allocate(DELTA(Nsite))
if(.not.allocated(error_ne)) allocate(error_ne(Nsite))
do i_site=1,Nsite
   ne0(i_site)=zero
   error_ne(i_site)=one
   DELTA(i_site)=zero
end do
Open(11,file="inputs.txt")
read(11,"(A)")ch
read(11,*)U
read(11,"(A)")ch
read(11,*)Delta1
read(11,"(A)")ch
read(11,*)Delta2
read(11,"(A)")ch
read(11,*)Delta3
read(11,"(A)")ch
read(11,*)NK
read(11,"(A)")ch
read(11,*)Nw
read(11,"(A)")ch
do i_site=1,Nsite
read(11,*)ne0(i_site)
end do
close(11)
open(12,file='CPA_process.txt')
write(12,*) "On-site Coulomb interaction between electtons:",U
write(12,*) "On-site potential of cluster1:",Delta1
write(12,*) "On-site potential of cluster2:",Delta2
write(12,*) "On-site potential of cluster3:",Delta3
write(12,*) "Total energy mesh at the energy axis:",Nw
write(12,*) "Total k mesh at the first Brillouin zone:",NK
write(12,*) "Initial occupations of sites:"
write(12,"(7f12.8)")ne0(:)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!---allocate the matrix---!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if(.not.allocated(DOS)) allocate(DOS(Nsite,Nw))
if(.not.allocated(TDOS)) allocate(TDOS(Nw))
if(.not.allocated(error_SE)) allocate(error_SE(Nsite,Nw))
if(.not.allocated(SE_old)) allocate(SE_old(Nsite,Nw))
if(.not.allocated(SE_new)) allocate(SE_new(Nsite,Nw))
if(.not.allocated(Omega)) allocate(Omega(Nw))
if(.not.allocated(Cavity_G)) allocate(Cavity_G(Nsite,Nw))
if(.not.allocated(Lattice_G)) allocate(Lattice_G(Nsite,Nw))
if(.not.allocated(Impurity_G)) allocate(Impurity_G(Nsite,Nw))
if(.not.allocated(Impurity_G_alpha)) allocate(Impurity_G_alpha(Nalpha,Nsite,Nw))
if(.not.allocated(P_alpha)) allocate(P_alpha(Nalpha,Nsite))
if(.not.allocated(HK)) allocate(HK(Nsite,Nsite,NK))
if(.not.allocated(HK_k)) allocate(HK_k(Nsite,Nsite,NK))
if(.not.allocated(INVHK_k)) allocate(INVHK_k(Nsite,Nsite,NK))
!!!!!!!!!!!!!!!!!!!!!!!!!------To read HK-----!!!!!!!!!!!!!!!!!!!!!!!!!
open(13,file='HK.dat')
do k=1,NK
   do i_site=1,Nsite
      do j_site=1,Nsite
         read(13,"(2f18.10)")HK(j_site,i_site,k)
      end do
   end do
end do
close(13)
!!!!!!!!!!!!!!!!!!!!!!!!!------To add CDW-----!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!------You can change CDW here-----!!!!!!!!!!!!!!!!!!
DELTA(1)=Delta2
DELTA(2)=Delta3
DELTA(3)=Delta2
DELTA(4)=Delta3
DELTA(5)=Delta3
DELTA(6)=Delta2
DELTA(7)=Delta2
DELTA(8)=Delta2
DELTA(9)=Delta1
DELTA(10)=Delta3
DELTA(11)=Delta3
DELTA(12)=Delta3
DELTA(13)=Delta1
DELTA(14)=Delta2
DELTA(15)=Delta3
DELTA(16)=Delta3
DELTA(17)=Delta3
DELTA(18)=Delta2
DELTA(19)=Delta2
DELTA(20)=Delta3
DELTA(21)=Delta3
DELTA(22)=Delta3
DELTA(23)=Delta2
DELTA(24)=Delta2
DELTA(25)=Delta1
DELTA(26)=Delta1
DELTA(27)=Delta2
DELTA(28)=Delta2
!!!!!!!!!!!!!!!!!!!!!!!!------To cal omega-----!!!!!!!!!!!!!!!!!!!!!!!!
do i=1,Nw
   Omega(i)=zero
   do i_site=1,Nsite
      error_SE(i_site,i)=one
      SE_old(i_site,i)=czero
      SE_new(i_site,i)=czero
   end do
end do
E_max=16.0d0
E_min=-13.0d0
domega=(E_max-E_min)/DBLE(Nw-1)
Ef=U/2.0d0
do i=1,Nw
   Omega(i)=E_min+DBLE(i-1)*domega+ione*eta
end do
!!!!!!!!!!!!!!!!!!!!!!!--guess the initial self-energy--!!!!!!!!!!!!!!!!!!!!!!
do i=1,Nw
   do i_site=1,Nsite
      SE_old(i_site,i)=SE_old(i_site,i)+ne0(i_site)*U
      SE_old(i_site,i)=SE_old(i_site,i)-ione
   end do
end do
!!!!!!!!!!!!!!!!!!!!--initialization of the probabilities--!!!!!!!!!!!!!!!!!!
do i_site=1,Nsite
   do i_alpha=1,Nalpha
      P_alpha(i_alpha,i_site)=zero
   end do
end do
Maxerror_SE=one
Maxerror_ne=one
!!!!!!!!!!!!!!!!!!!!!!!!!!!--start SCF calculation--!!!!!!!!!!!!!!!!!!!!!!!!!!
iter=0
80 continue
iter=iter+1
write(12,*)"self-consistent steps",iter
!!!!!!!!!!!!!!!!!!!!!--initialization the greenfunctions--!!!!!!!!!!!!!!!!!!!!
do i=1,Nw
   TDOS(i)=zero
   do i_site=1,Nsite
      Cavity_G(i_site,i)=czero
      Impurity_G(i_site,i)=czero
      DOS(i_site,i)=zero
      Lattice_G(i_site,i)=czero
      do i_alpha=1,Nalpha
         Impurity_G_alpha(i_alpha,i_site,i)=czero
      end do
   end do
end do
!!!!!!!!!!!!!!!!!!!!--initialization of the final occupation--!!!!!!!!!!!!!!!!
do i_site=1,Nsite
      net(i_site)=zero
      error_ne(i_site)=zero
end do
! Call system_clock(t3)
!!!!!!!!!!!!!!!--calculate the effective medium greenfunctions--!!!!!!!!!!!!!
!$Omp parallel do num_threads(6) default(shared) private(k, i_site, j_site, i)&
!$Omp firstprivate(Omega, SE_old, DELTA, HK_k, INVHK_k, HK)&
!$Omp shared(NK,Ef)
!!$Omp do schedule(guided)
do k=1,NK
   do i_site=1,Nsite
      do j_site=1,Nsite
         HK_k(j_site,i_site,k)=-HK(j_site,i_site,k)
      end do
   end do
   do i=1,Nw
      do i_site=1,Nsite
         HK_k(i_site,i_site,k)=Omega(i)-HK(i_site,i_site,k)-SE_old(i_site,i)-Ef-DELTA(i_site)
      end do
      
      call INVERT(Nsite,HK_k(:,:,k),INVHK_k(:,:,k))
      do i_site=1,Nsite  
         Lattice_G(i_site,i)=Lattice_G(i_site,i)+1.0d0/(DBLE(NK))*INVHK_k(i_site,i_site,k)
      end do
   end do
end do
!$Omp end parallel do
! Call system_clock(t4)
! write(*, *) (t4-t3)/1000
!!!!!!!!!!!!!!!!--<<<Calculating fermi surface>>>--!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!--<<<Calculating correlation function>>>--!!!!!!!!!!!!!
do i=1,Nw
   do i_site=1,Nsite
      DOS(i_site,i)=DOS(i_site,i)-(1.0d0/pi)*AIMAG(Lattice_G(i_site,i))
   end do
end do
Tn0=zero
Tnt=zero
do i=2,Nw
   do i_site=1,Nsite
        Tn0=Tn0+0.5d0*(DOS(i_site,i)+DOS(i_site,i-1))*domega
   end do     
end do
ifermi=0
do i=2,Nw
   if(Tnt.LT.Tn0/2.0d0)then
     ifermi=ifermi+1
     do i_site=1,Nsite
        Tnt=Tnt+0.5d0*(DOS(i_site,i)+DOS(i_site,i-1))*domega
        net(i_site)=net(i_site)+0.5d0*(DOS(i_site,i)+DOS(i_site,i-1))*domega
     end do
   else  
     goto 81
   end if
end do
81 continue
Ef=REAL(omega(ifermi)+omega(ifermi-1))/2.0d0
write(*, *)iter,Ef
!!!!!!!!!!!!!!!!!!!!--construct the cavity greenfunctions--!!!!!!!!!!!!!!!!!!!
do i=1,Nw
   do i_site=1,Nsite
      Cavity_G(i_site,i)=1.0d0/((1.0d0/Lattice_G(i_site,i)+SE_old(i_site,i)))
   end do
end do
!!!!!!!!!!!!!!!!!!!!--The calculation of the probabilities--!!!!!!!!!!!!!!!!!!
do i_site=1,Nsite
   P_alpha(1,i_site)=ne0(i_site)
   P_alpha(2,i_site)=1.0d0-ne0(i_site)
end do
!!!!!!!!!!!!!!!!!!!--calculate the impurity green function--!!!!!!!!!!!!!!!!!!
do i=1,Nw
   do i_site=1,Nsite
      do i_alpha=1,Nalpha
         if(i_alpha==1)then
           Impurity_G_alpha(i_alpha,i_site,i)=1.0d0/(1.0d0/Cavity_G(i_site,i)+Ef-U)
         else
           Impurity_G_alpha(i_alpha,i_site,i)=1.0d0/(1.0d0/Cavity_G(i_site,i)+Ef-0)
         end if
      end do
      do i_alpha=1,Nalpha     
         Impurity_G(i_site,i)=Impurity_G(i_site,i)+P_alpha(i_alpha,i_site)*Impurity_G_alpha(i_alpha,i_site,i)
      end do
   end do
end do
!open(13,file='Test.txt')
!do i=1,Nw
!   write(13,"(4f10.6)")Impurity_G(1,i),Impurity_G(28,i)
!end do
!close(13)
!stop
!!!!!!!!!!!!--calculate the new self-energy according to Dyson EQ--!!!!!!!!!!!
do i=1,Nw
   do i_site=1,Nsite
      SE_new(i_site,i)=1.0d0/(Cavity_G(i_site,i))-1.0d0/(Impurity_G(i_site,i))
   end do
end do
write(12,*)"Total occupation:",Tn0
write(12,*)"Total half occupation:",Tnt
write(12,*)"=========================================================="
write(12,*)"Initial occupations:"
write(12,"(7f12.6)")ne0(:)
write(12,*)"Final occupations:"
write(12,"(7f12.6)")net(:)
!!!!!!!!--<<<Comparing initial and final occupation numbers>>>--!!!!!!!!
do i_site=1,Nsite
   error_ne(i_site)=ABS(ne0(i_site)-net(i_site))
end do
do i_site=2,Nsite
Maxerror_ne=MAX(error_ne(i_site),error_ne(i_site-1))
end do
!!!!!!!!!--<<<Comparing initial and final self-energy>>>--!!!!!!!!!!
do i=1,Nw
   do i_site=1,Nsite
      error_SE(i_site,i)=ABS(SE_new(i_site,i)-SE_old(i_site,i))
   end do
end do
Maxerror_SE=error_SE(1,1)
do i=1,Nw
   do i_site=1,Nsite
      Maxerror_SE=MAX(Maxerror_SE,error_SE(i_site,i))
   end do
end do
write(12,*)"Ef=",Ef,"Eferror=",REAL(omega(ifermi)+omega(ifermi-1))/2.0d0
write(12,*)"Maximal occupation numbers error:",Maxerror_ne
write(12,*)"Maximal self-energy error:",Maxerror_SE
write(12,*)"=========================================================="
write(12,*)
if((Maxerror_ne.gt.max_eps_n).OR.(Maxerror_SE.gt.max_eps_SE))then
  do i=1,Nw
     do i_site=1,Nsite
        SE_old(i_site,i)=SE_new(i_site,i)
     end do
  end do
  do i_site=1,Nsite
     ne0(i_site)=0.5d0*ne0(i_site)+0.5d0*net(i_site)
  end do
  if(iter.gt.itermax)then
    write(10,*)"convergence is not arrived."
    goto 82
  end if
  goto 80
else
  write(10,*)"convergence has arrived."
end if
82  continue
!!!!!!!!!!!!!!---re-calculate DOS of different sites and TDOS---!!!!!!!!!!!!!!
do i=1,Nw
   TDOS(i)=zero
   do i_site=1,Nsite
      DOS(i_site,i)=zero
   end do
end do
do i=1,Nw
   do i_site=1,Nsite
      DOS(i_site,i)=DOS(i_site,i)-(1.0d0/pi)*AIMAG(Lattice_G(i_site,i))
      TDOS(i)=TDOS(i)+DOS(i_site,i)
   end do
end do
open(14,file='TDOS.txt')
do i=1,Nw
   write(14,"(2f14.10)")REAL(omega(i))-Ef,TDOS(i)
end do 
close(14)
open(15,file='ATOM_DOS.txt')
do i=1,Nw
   write(15,"(29f10.6)")REAL(omega(i))-Ef,DOS(:,i)
end do
close(15)
end
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!------Subroutine------!!!!!!!!!!!!!!!!!!!!!
! !>>> to calculat the inverse of the matrix for complex matrix !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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