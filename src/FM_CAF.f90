      program oneband_HF
!      use msimsl
      IMPLICIT REAL(8) (a-h,o-z)
      PARAMETER  (Nsite=2, Nband=1,Nkx=4096,Nky=Nkx,Nk=4*Nkx*Nky)
      PARAMETER  (Bolzman=8.61734315d-5,Temp=1.0d-7)
!      PARAMETER (betaT=1.0d0/Bolzman/temp)
      PARAMETER (betaT=1.0d0/Temp)
! half filling: number of electrons = number of bands in each unit cell
      PARAMETER (Rfilling=2.4d0)
      dimension eigen_up(Nsite,Nk,Nband),EVAL(4),eigen_dn(Nsite,Nk,Nband)
!!!this group dimenson is to used to test
     ! real group1(NK)
     ! real group2(NK)
! define pi
      pi=dasin(1.0d0)*2.0d0
! initialize the hoping
      t1=-1.0d0
      t2=-0.8d0
      open(31,file="checking_T_CAFoneBandData_t20707_U_4_5_4_7_T057.txt")
      write(31,*)"t2=",t2
      write(31,*)"        Free_Energy            ","magnetization            "," metal_flag "
      Ustep=0.01d0
      do  i_U=100,800,20
      U=Ustep*dble(i_U)
! initialize magnetization on each site
! initialize density on each site
      val_mag_delta=0.5d0
      val_den_avg=0.6d0
! selfconsistency begins
! if converged, Nconvergence=1, if after 400 loops, unconverged, then Nconvergence=0, stop interation
      Nconvergence=1
! counting the interation loop
      iteration=0
99    iteration=iteration+1
! initialize Hamiltonian matrix
      stepx=pi/dble(Nkx)
      stepy=pi/dble(Nky)
      upmax=-1.0d10
      dnmax=1.0d10
      do i_kx=-Nkx,Nkx-1
       xk=stepx*dble(i_kx)
       cx=dcos(xk)
       cx2=dcos(xk/2.0d0)
       do i_ky=-Nky,Nky-1
        yk=stepy*dble(i_ky)
        cy=dcos(yk)
        cy2=dcos(yk/2.0d0)
! define epsilon
        epsi11=2.0d0*t2*(cx+cy)
        epsi12=4.0d0*t1*cx2*cy2 
!! spin up
        EVAL(1)=epsi11+U*(val_den_avg-val_mag_delta/2)-dsqrt(epsi12*epsi12)
        EVAL(2)=epsi11+U*(val_den_avg-val_mag_delta/2)+dsqrt(epsi12*epsi12)
        EVAL(3)=epsi11+U*(val_den_avg+val_mag_delta/2)-dsqrt(epsi12*epsi12)
        EVAL(4)=epsi11+U*(val_den_avg+val_mag_delta/2)+dsqrt(epsi12*epsi12)
! search maximum and minimum of the bands
        upmax=DMAX1(upmax,EVAL(1),EVAL(2),EVAL(3),EVAL(4))
        dnmax=DMIN1(dnmax,EVAL(1),EVAL(2),EVAL(3),EVAL(4))
! save all the eigenvalues
        i_pos=(i_kx+Nkx)*Nky*2+i_ky+Nky+1
         eigen_up(1,i_pos,1)=EVAL(1)
         eigen_up(2,i_pos,1)=EVAL(2)
         eigen_dn(1,i_pos,1)=EVAL(3)
         eigen_dn(2,i_pos,1)=EVAL(4)
       end do
      end do
! search the Fermi surface
6     Fermilevel=(upmax+dnmax)/2.0d0
      fillingup1=0.0d0
      fillingdn1=0.0d0
      do i_kx=1,2*Nkx
       do i_ky=1,2*Nky
        i_pos=(i_kx-1)*2*Nkx+i_ky
        do i_site=1,Nsite
! spin up
         if(betaT*(eigen_up(i_site,i_pos,1)-Fermilevel).gt.3.0d1)then
           fermi1=0.0d0
         else
          if(betaT*(eigen_up(i_site,i_pos,1)-Fermilevel).lt.-3.0d1)then
           fermi1=1.0d0
          else
           fermi1=1.0d0/(1.0d0+dexp(betaT*(eigen_up(i_site,i_pos,1)-Fermilevel)))
          end if
         end if
         fillingup1=fillingup1+fermi1/dble(2*Nkx*2*Nky)
! spin dn
         if(betaT*(eigen_dn(i_site,i_pos,1)-Fermilevel).gt.3.0d1)then
           fermi2=0.0d0
         else
          if(betaT*(eigen_dn(i_site,i_pos,1)-Fermilevel).lt.-3.0d1)then
           fermi2=2.0d0
          else
           fermi2=2.0d0/(1.0d0+dexp(betaT*(eigen_dn(i_site,i_pos,1)-Fermilevel)))
          end if
         end if
         fillingdn1=fillingdn1+fermi2/dble(2*Nkx*2*Nky)
        end do
       end do
      end do
! and here number of electrons in spin up	channel is equal to that in spin dn channel	fillingdn1=fillingup1
      filling=fillingup1+fillingdn1
      if((filling-Rfilling).lt.-1.0d-11)then
       dnmax=Fermilevel
       goto 6
      else
       if((filling-Rfilling).gt.1.0d-11)then
        upmax=Fermilevel
        goto 6
       end if
      end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
! begin to calculate the magnetization on each site
      val_mag_delta_new=0.0d0
      stepx=pi/dble(Nkx)
      stepy=pi/dble(Nky)
      fillingup1=0.0d0
      fillingdn1=0.0d0
      do i_kx=-Nkx,Nkx-1
       xk=stepx*dble(i_kx)
       cx=dcos(xk)
       cx2=dcos(xk/2.0d0)
       do i_ky=-Nky,Nky-1
        yk=stepy*dble(i_ky)
        cy=dcos(yk)
        cy2=dcos(yk/2.0d0)
! define epsilon
        epsi11=2.0d0*t2*(cx+cy)
        epsi12=4.0d0*t1*cx2*cy2 
! spin up
        EVAL(1)=epsi11+U*(val_den_avg-val_mag_delta)-dsqrt(epsi12*epsi12)
        EVAL(2)=epsi11+U*(val_den_avg-val_mag_delta)+dsqrt(epsi12*epsi12)
! spin dn
        EVAL(3)=epsi11+U*(val_den_avg+val_mag_delta)-dsqrt(epsi12*epsi12)
        EVAL(4)=epsi11+U*(val_den_avg+val_mag_delta)+dsqrt(epsi12*epsi12)
        do i_site=1,Nsite
! spin up
         if(betaT*(EVAL(i_site)-Fermilevel).gt.3.0d1)then
              fermi1=0.0d0
         else
          if(betaT*(EVAL(i_site)-Fermilevel).lt.-3.0d1)then
           fermi1=1.0d0
          else
           fermi1=1.0d0/(1.0d0+dexp(betaT*(EVAL(i_site)-Fermilevel)))
          end if
         end if
      fillingup1=fillingup1+fermi1
!!!!!!!!!!!!!!progress mark
!         if(dabs(epsi12).eq.0.0d0.and.val_mag_delta.eq.0.0d0)then
!          val_mag_delta_new=val_mag_delta
!         else
!          if(i_site.eq.1)then
!            val_mag_delta_new=val_mag_delta_new+val_mag_delta/dsqrt(val_mag_delta*val_mag_delta+epsi12*epsi12)*fermi1
!          else
!           if(i_site.eq.2)then
!            val_mag_delta_new=val_mag_delta_new-val_mag_delta/dsqrt(val_mag_delta*val_mag_delta+epsi12*epsi12)*fermi1
!           end if
!          end if
!         end if
        end do
        do i_site=3,4
! spin dn
         if(betaT*(EVAL(i_site)-Fermilevel).gt.3.0d1)then
              fermi2=0.0d0
         else
          if(betaT*(EVAL(i_site)-Fermilevel).lt.-3.0d1)then
           fermi2=1.0d0
          else
           fermi2=1.0d0/(1.0d0+dexp(betaT*(EVAL(i_site)-Fermilevel)))
          end if
         end if
      fillingdn1=fillingdn1+fermi2
        end do
       end do
      end do
	val_mag_delta_new=fillingup1-fillingdn1



! calculate the delta
      val_mag_delta_new_rec=val_mag_delta_new
      errinmag1=dabs(val_mag_delta_new-val_mag_delta)
      if(errinmag1.gt.1.0d-10)then
       val_mag_delta=val_mag_delta_new
       if(iteration.lt.400)then
        goto 99
       else
        Nconvergence=0
       end if 
      end if

!to see if it is a metal or what
  
    metal_flag=1  !!metal_flag��������ǲ��ǽ���������ǣ���ֵΪ1������Ϊ0  
!  a_max_site1=eigen_up(1,1,1)
!  a_min_site2=eigen_up(2,1,1)
!     do i_kx=-Nkx,Nkx-1
!       do i_ky=-Nky,Nky-1
!          i_pos=(i_kx+Nkx)*Nky*2+i_ky+Nky+1
!          if(a_max_site1<eigen_up(1,i_pos,1))then 
!             a_max_site1=eigen_up(1,i_pos,1)
!          end if
!          if(a_min_site2>eigen_up(2,i_pos,1))then
!             a_min_site2=eigen_up(2,i_pos,1)
!          end if
!       end do
!    end do 
!   if(a_max_site1>=a_min_site2)then
!     metal_flag=1
!   end if
!   if(a_max_site1<a_min_site2)then
!     metal_flag=0
!   end if
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
 !  metal_flag=0
 !  do i_kx=-Nkx,Nkx-1
 !   do i_ky=-Nky,Nky-1
 !    i_pos=(i_kx+Nkx)*Nky*2+i_ky+Nky+1
 !    group1(i_pos)=eigen_up(1,i_pos,1)
 !    group2(i_pos)=eigen_up(2,i_pos,1)
 !   end do
 !  end do
!  max_site1=maxval(group1)
!  min_site2=minval(group2)
!  if(max_site1>=min_site2)then
!    metal_flag=1
!  end if
!  if(max_site1<min_site2)then
!    metal_flag=0
!  end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


 
  

! write out the DOS and band structure
  call cal_BSDOS(eigen_up,eigen_dn,Fermilevel,upmax,dnmax)
! calculate the free energy
! first from grand canonical potential:
      pot_grand=0.0d0
      do i_kx=1,2*Nkx
       do i_ky=1,2*Nky
        i_pos=(i_kx-1)*2*Nkx+i_ky
! ln(1+exp(-\beta*(e(k)-\mu)))
        do i_site=1,Nsite
! spin up
         if(betaT*(eigen_up(i_site,i_pos,1)-Fermilevel).gt.3.0d1)then
           potkinetic1=0.0d0
         else
          if(betaT*(eigen_up(i_site,i_pos,1)-Fermilevel).lt.-3.0d1)then
           potkinetic1=-betaT*(eigen_up(i_site,i_pos,1)-Fermilevel)
          else
           potkinetic1=dlog(1.0d0+dexp(-betaT*(eigen_up(i_site,i_pos,1)-Fermilevel)))
          end if
         end if
         pot_grand=pot_grand+potkinetic1
        end do
        do i_site=1,Nsite
! spin dn
         if(betaT*(eigen_dn(i_site,i_pos,1)-Fermilevel).gt.3.0d1)then
           potkinetic2=0.0d0
         else
          if(betaT*(eigen_dn(i_site,i_pos,1)-Fermilevel).lt.-3.0d1)then
           potkinetic2=-betaT*(eigen_dn(i_site,i_pos,1)-Fermilevel)
          else
           potkinetic2=dlog(1.0d0+dexp(-betaT*(eigen_dn(i_site,i_pos,1)-Fermilevel)))
          end if
         end if
         pot_grand=pot_grand+potkinetic2
        end do
       end do
      end do
      pot_grand=-pot_grand/betaT/dble(2*Nkx*2*Nky)
! second: from rest terms
      pot_res=-2.0d0*U*val_den_avg*val_den_avg+0.5*U*(val_mag_delta_new_rec*val_mag_delta_new_rec)
! third: from chemical potential
      pot_chem=Fermilevel*Rfilling
! total potential
      pot_tot=pot_grand+pot_res+pot_chem
      write(31,*)pot_tot/2,val_mag_delta_new_rec,metal_flag
      end do
      close(31) 
      end program oneband_HF !end program
     
     
     
     
      SUBROUTINE ResetLocation(N,EVAL,LDEVEC,EVEC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION   EVAL(N)
      complex(8)  EVEC(LDEVEC,N),y
      DO I=1,N-1
       DO J=I+1,N
        IF(EVAL(J).lt.EVAL(I))THEN
         x=EVAL(I)
         EVAL(I)=EVAL(J)
         EVAL(J)=x
         DO K=1,ldevec
          y=EVEC(K,I)
          EVEC(K,I)=EVEC(K,J)
          EVEC(K,J)=y
         end do
        END IF
       end do
      end do
      return
      end

      subroutine cal_BSDOS(eigen_up,eigen_dn,Fermilevel,upmax,dnmax)
      IMPLICIT REAL(8) (a-h,o-z)
      PARAMETER  (Nsite=2, Nband=1,Nkx=4096,Nky=Nkx,Nk=4*Nkx*Nky)
      PARAMETER  (Bolzman=8.61734315d-5,Temp=1.0d-5)
      PARAMETER (betaT=1.0d0/Temp)
      dimension eigen_up(Nsite,Nk,Nband)
      ! REAL *8, allocatable :: eigen_up(:, :, :)
      pi=dasin(1.0d0)*2.0d0
      stepx=pi/dble(Nkx)
      stepy=pi/dble(Nky)
! plot 3D band surface
      open(30, file="BSsurface_t2025_U0_0.dat")
      do i_kx=1,2*Nkx
       do i_ky=1,2*Nky
        write(30,30)i_kx, i_ky,(eigen_up(i_site,(i_kx-1)*2*Nkx+i_ky,1)-Fermilevel,i_site=1,Nsite) 
       end do
      end do
      close(30)
30    format(2I5, 2f10.5)
! plot DOS with common Fermi level
      open(80,file="DOS_t2025_U0_0.dat")
      do omiga=-20,20,0.01d0
       dos1=0.0d0
       dos2=0.0d0
       do i_kx=1,2*Nkx
        do i_ky=1,2*Nky
         i_pos=(i_kx-1)*2*Nkx+i_ky
         do i_site=1,Nsite
          dos1=dos1+dimag(dcmplx(1.0d0,0.0d0)/(dcmplx(omiga+Fermilevel-eigen_up(i_site,i_pos,1),0.01d0)))
          dos2=dos2+dimag(dcmplx(1.0d0,0.0d0)/(dcmplx(omiga+Fermilevel-eigen_dn(i_site,i_pos,1),0.01d0)))
         end do
        end do
       end do
! factor of 2 means spin up and spin dn
       write(80,80)omiga,-dos1/Pi/dble(2*Nkx*2*Nky),-dos2/Pi/dble(2*Nkx*2*Nky),-(dos1+dos2)/Pi/dble(2*Nkx*2*Nky)
      end do
80    format(2f10.5)
      close(80)
! plot band structure
      open(60,file="BSline_t2025_U0_0.dat")
      open(50,file="specialk_t2025_U0_0.dat")
      i_step=0
! from -\pi,0 -> 0,0
60    format(3f10.5)
      do i_kx=1,Nkx+1
       i_pos=(i_kx-1)*2*Nkx+Nkx+1
       r_pos=dble(i_kx-1)*stepx
       write(60,60)r_pos,(eigen_up(i_site,i_pos,1)-Fermilevel,i_site=1,Nsite)
      end do
      r_pos_rec=r_pos
      write(50,60)r_pos_rec,upmax,dnmax 
! from 0,0 -> -\pi,-pi
      do i_kx=Nkx,1,-1
       i_ky=i_kx
       i_pos=(i_kx-1)*2*Nkx+i_ky
       r_pos=r_pos_rec+dble(Nkx+1-i_kx)*dsqrt(stepx*stepx+stepy*stepy)
       write(60,60)r_pos,(eigen_up(i_site,i_pos,1)-Fermilevel,i_site=1,Nsite)
      end do
      r_pos_rec=r_pos
      write(50,60)r_pos_rec,upmax,dnmax 
! from -\pi,-pi -> -\pi,0
      do i_ky=2,Nky+1
       i_pos=i_ky
       r_pos=r_pos_rec+dble(i_ky-1)*stepy
       write(60,60)r_pos,(eigen_up(i_site,i_pos,1)-Fermilevel,i_site=1,Nsite)
      end do
      r_pos_rec=r_pos
      write(50,60)r_pos_rec,upmax,dnmax
      close(50) 
      close(60)
      return
      end 
