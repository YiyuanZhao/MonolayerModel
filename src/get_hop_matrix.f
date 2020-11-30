c************************************************************************!!!
c the main parameters used in the calculations !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c************************************************************************!!!
      module param
	implicit double precision(a-h,o-z)
      integer  nwann,nhalfwan,NRdg
	complex*16,allocatable::ham_r    (:,:,:)
	integer   ,allocatable::Rnspac   (:,:)
	character(256) filename,casename
	integer   korigin,Kwien2k,korigextend
	real*8,allocatable::bandklist(:,:)
	real*8,allocatable::bandklist_extend(:,:) 
c the fractional translation vector !!!
      real*8    vec_trans(3),unfoldvec(3)
c the transformation matrix in orbital space !!!
	real*8,allocatable::W_Ci(:,:)
c band structure in full brillouin zone !!!
c	integer   NKX,NKY,NKZ
      integer   NKX,NKY,DNKY,TNKY
      parameter(NKX=1201,NKY=600,DNKY=1201,TNKY=DNKY+NKY)
	real*8    pi,twopi
	real*8    veclat(3)
	end
c************************************************************************!!!
c************************************************************************!!!
	program get_hop_matrix
	use param
	implicit double precision(a-h,o-z)
	logical exist1

	filename="grap.dat"
c constants !!!
      pi=dasin(1.0d0)*2.0d0
      twopi=2.0d0*pi

c to get hopping matrix t_mn(R) !!!
	call read_hr()
	call generat_k()

c to get K_path !!!

	Kwork=korigin
	call cal_origin_band(bandklist,Kwork)
c
      call cal_hk_allzone()

	end
c************************************************************************!!!
c************************************************************************!!!
c to read hamitonial in wannier representations !!!!!!!!!!!!!!!!!!!!!!!!!!!!
c************************************************************************!!!
c************************************************************************!!!
	subroutine read_hr()
	use param
	implicit double precision(a-h,o-z)
	character(30) ch
	integer,allocatable::ndegen(:)
      open(10,file=filename)
	    read(10,*)ch
	    read(10,*)nwann
	    write(*,*)"nwann=",nwann
	    nhalfwan=nwann/2
	    read(10,*)NRdg
	    write(*,*)"NRdg=",NRdg
c to read out degenerate of every real point!
          if(.not.allocated(ndegen))allocate(ndegen(NRdg))
	    read(10,"(15i5)")(ndegen(isit),isit=1,NRdg)
          if(allocated(ndegen))deallocate(ndegen)
	    if(.not.allocated(Rnspac))allocate(Rnspac(3,NRdg))
	    if(.not.allocated(ham_r ))allocate(ham_r(nwann,nwann,NRdg))
c to read out hamitonial in real space !
          do isit=1,NRdg
             do iorb=1,nwann
	          do jorb=1,nwann
 	             read(10,"(5i5,2f12.6)")NRx,NRy,NRz,korb,lorb,
     $                                    ham_r(jorb,iorb,isit)
       if(dabs(dimag(ham_r(jorb,iorb,isit))).lt.2.d-6)then
      ham_r(jorb,iorb,isit)=dcmplx(DBLE(ham_r(jorb,iorb,isit)),0.0d0)
                        endif
	          end do
		   end do
		   Rnspac(1,isit)=NRx
		   Rnspac(2,isit)=NRy
		   Rnspac(3,isit)=NRz
		end do
	close(10)
	write(*,*)"have read case_hr.dat successfully..."
	return
	end
c************************************************************************!!!
c************************************************************************!!!
c************************************************************************!!!
c************************************************************************!!!
	subroutine cal_origin_band(bandkpath,Ktot)
	use param
	implicit double precision(a-h,o-z)
      real*8,allocatable:: eigenval(:,:)
	real*8     eval(nwann),bandkpath(3,Ktot)
	complex*16 ham_work(nwann,nwann),hak(nwann,nwann)
	real*8     wk(3)

	       
	    if(.not.allocated(eigenval))allocate(eigenval(nwann,Ktot))
	         ktotal=Ktot
		   do kloop=1, ktotal 
	          do korb=1,nwann
	             eigenval(korb,kloop)=0.0d0
			  end do
	       end do

	       do kloop=1,ktotal
 
                do ik=1,3
	             wk(ik)=bandkpath(ik,kloop)
		      end do

	   call cal_spectrum(wk,hak,ham_work,eval,nwann,ham_r,Rnspac,NRdg)

             do iband=1,nwann
	          eigenval(iband,kloop)=eval(iband)
		   end do

		   end do

	       open(10,file="spectrum_orig.dat")
	       do iband=1,nwann
	          do kloop=1,ktotal
	             write(10,"(i5,f12.6)")kloop,eigenval(iband,kloop)
			  end do
	          write(10,*)
		   end do
	       close(10)
	return
	end
c**********************************************************************************!!!
c**********************************************************************************!!!
        subroutine cal_hk_allzone()
        use param
        implicit double precision(a-h,o-z)
        real*8  wk(3),XK,YK
        complex*16,allocatable :: HK(:,:,:,:)
        complex*16 ham_work(nwann,nwann),hak(nwann,nwann)
        A=2.46D0/sqrt(3.0D0)                  !A value: relaxed C-C or B-N bond
        pi=dasin(1.0d0)*2.0d0
c********one sixth Brilliouin zone****************************************
        DKx=(2.0d0*pi/(3.0*A))/dble(NKX-1)
        DKy=(2.0d0*pi/(3.0*SQRT(3.0D0)*A))/dble(NKY)
c********one sixth Brilliouin zone****************************************
        if(.not.allocated(HK))allocate(HK(nwann,nwann,1:TNKY,1:(NKX+1)))
        OPEN(11,file='HK.dat')
        OPEN(18,file='nk.dat')
        OPEN(17,file='kx+ky.dat')
         NKSUM=0
         DO KX=1,NKX+1
           XK=dble(KX-1)*Dkx
           DO KY=1,TNKY
              ndy=KY-2*NKY-1
              YK=dble(KY-2.0d0*NKY-1)*DKy
              NKSUM=NKSUM+1
              wk(1)=XK
              wk(2)=YK
           write(17,"(2f12.6)")XK,YK
                  wk(3)=0.0d0
        call cal_spectrum(wk,hak,ham_work,eval,nwann,ham_r,Rnspac,NRdg)
           write(17,"(2f12.6)")XK,YK
                  do iorb=1,nwann
                     do jorb=1,nwann
                        HK(jorb,iorb,KY,KX)=hak(jorb,iorb)
                        write(11,"(2f12.6)") HK(jorb,iorb,KY,KX)
                     enddo
                  enddo
            enddo
         enddo
          write(11,*)"The-End!"
          write(*,*)"total k points",NKSUM
        close(17)
        close(18)
        close(11)
        return
        end
c********all Brilliouin zone****************************************
c**********************************************************************************!!!
c**********************************************************************************!!!
cto calculate bandstructure of given point based on unit cell used in the calclation !
c**********************************************************************************!!!
c**********************************************************************************!!!
      subroutine cal_spectrum(wk,hak,hamvec,eig,nwann,ham_r,Rnspac,NRdg)
	implicit double precision(a-h,o-z)
      real*8     wk(3),Rn(3)
	complex*16 structphase
	complex*16 hamk_wk(nwann,nwann),ham_work(nwann,nwann)
	complex*16 hak(nwann,nwann)
	complex*16 hamvec(nwann,nwann),ham_r(nwann,nwann,NRdg)
	real*8     eigval(nwann),A1(3),A2(3),A3(3)
	real*8                   B1(3),B2(3),B3(3)
	real*8     temp(3),newk(3)
	integer    Rnspac(3,NRdg)
	character(30) ch_state


      pi=dasin(1.0d0)*2.0d0
      twopi=2.0d0*pi
C     relaxed lattice constant need to modify!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	A1(1)=2.1304225d0
	A1(2)=-1.23d0
	A1(3)=0.0D0

	A2(1)=0.0D0
	A2(2)=2.46D0
	A2(3)=0.0D0

	A3(1)=0.0D0
	A3(2)=0.0D0
	A3(3)=20.3D0

	B1(1)=0.469390463d0*twopi
	B1(2)=0.0D0
	B1(3)=0.0D0

	B2(1)=0.234695231d0*twopi
	B2(2)=0.406504065d0*twopi
	B2(3)=0.0D0

	B3(1)=0.0d0
	B3(2)=0.0d0
	B3(3)=0.0492610837d0*twopi


	               do iorb=1,nwann
	                  do jorb=1,nwann
	                     hamk_wk(jorb,iorb)=dcmplx(0.0d0,0.0d0)
					  end do
				   end do

				  do isit=1,NRdg
c for the translation of crystal based used unit cell !!!
                      
c for the position of nearest cell !!!
					 do jsit=1,3
	                    temp(jsit)=0.0d0
					 end do

c for A1
					 do jsit=1,3
						temp(jsit)=temp(jsit)+
     $                               dble(Rnspac(1,isit))*A1(jsit)
					 end do
c for A2
					 do jsit=1,3
						temp(jsit)=temp(jsit)+
     $                               dble(Rnspac(2,isit))*A2(jsit)
					 end do
c for A3
					 do jsit=1,3
						temp(jsit)=temp(jsit)+
     $                               dble(Rnspac(3,isit))*A3(jsit)
					 end do

	                 do jsit=1,3
	                    Rn(jsit)=temp(jsit)
					 end do

c for k vectors !
					 do jsit=1,3
	                    temp(jsit)=0.0d0
					 end do
c                   
c for B1
C					 do jsit=1,3
C						temp(jsit)=temp(jsit)+
C     $                               dble(wk(1))*B1(jsit)
C					 end do
c for B2
C					 do jsit=1,3
C						temp(jsit)=temp(jsit)+
C     $                               dble(wk(2))*B2(jsit)
C					 end do
c for B3
C					 do jsit=1,3
C						temp(jsit)=temp(jsit)+
C     $                               dble(wk(3))*B3(jsit)
C					 end do
	                    
	                 do jsit=1,3
	                    newk(jsit)=temp(jsit)+wk(jsit)
					 end do

						dotproduct=0.0d0
	                 do jsit=1,3
	                     dotproduct=dotproduct+newk(jsit)*Rn(jsit)
					 end do

c the K vector using 2pi/a , 2pi/b and 2pi/c as unit !
c the R vector using a , b and c as unit !!!!!!!!!!!!!

                   dotproduct=dotproduct         
	             structphase=dcmplx(dcos(dotproduct),dsin(dotproduct))

c for original part alpha !
c to take the original phase factor !

	               do iorb=1,nwann
	                  do jorb=1,nwann
	                     hamk_wk(jorb,iorb)=hamk_wk(jorb,iorb)+
     $       	             ham_r(jorb,iorb,isit)*structphase
					  end do
				   end do

c *******************************!
c have obtained H_wk !!!!!!!!!!!!!
c *******************************!
				  end do

	     do iorb=1,nwann
	        if(dabs(dimag(hamk_wk(iorb,iorb))).lt.1.d-15)then
	         hamk_wk(iorb,iorb)=dcmplx(dble(hamk_wk(iorb,iorb)),0.0d0)
			endif
	     end do
c to get ham_work !!!

           do iorb=1,nwann
	        do jorb=1,nwann
	           ham_work(jorb,iorb)=dcmplx(0.0d0,0.0d0)
			end do
		 end do

           do iorb=1,nwann
	        do jorb=1,nwann
	     ham_work(jorb,iorb)=ham_work(jorb,iorb)+hamk_wk(jorb,iorb)
			end do
		 end do
C
	     do iorb=1,nwann
	        do jorb=1,nwann
	     hak(jorb,iorb)=ham_work(jorb,iorb)
			end do
		 end do
c to diagonalize hamk_work !
	     ch_state='V'
           call cal_eigenVS(nwann,ham_work,eig,ch_state)

	     do lorb=1,nwann
	        do korb=1,nwann
	           hamvec(korb,lorb)=ham_work(korb,lorb)
			end do
		 end do

	return
	end
c************************************************************************!!!
c************************************************************************!!!
c to calculate band structure with given K-path !!!!!!!!!!!!!!!!!!!!!!!!!!!!
c************************************************************************!!!
c************************************************************************!!!

c************************************************************************!!!
c************************************************************************!!!
c to calculate band structure with given K-path !!!!!!!!!!!!!!!!!!!!!!!!!!!!
c************************************************************************!!!
c************************************************************************!!!

c************************************************************************!!!
c************************************************************************!!!
c************************************************************************!!!
      subroutine generat_k()
	use param
	implicit double precision(a-h,o-z) 
C	parameter(A=1.43590)
	parameter(NSUBR1=85)
	parameter(NSUBR2=42)
	parameter(NSUBR3=73)
c original G point
	real*8 DG(3)
c original X point
      real*8 DK(3)
	real*8 DM(3)
      A=2.46D0/sqrt(3.0D0)                     !A value: relaxed C-C or B-N bond
	pi=dasin(1.0d0)*2.0d0

c to determine korigin
        
	DG(1)= 0.00d0
	DG(2)= 0.00d0
	DG(3)= 0.00d0

	DK(1)= 2.0d0*pi/(3.0d0*A)
	DK(2)= 2.0d0*pi/3.0d0/sqrt(3.0d0)/A          
	DK(3)= 0.00d0

	DM(1)= 2.0d0*pi/(3.0d0*A)
	DM(2)= 0.00d0
	DM(3)= 0.00d0


      korigin=NSUBR1+NSUBR2+NSUBR3+1
	if(.not.allocated(bandklist))allocate(bandklist(3,korigin))


	 kcount=0
       stepkx=(DK(1)-DG(1))/(NSUBR1)
       stepky=(DK(2)-DG(2))/(NSUBR1)
       stepkz=(DK(3)-DG(3))/(NSUBR1)
C G-->K
       do  ik=1,NSUBR1
	     xk=DG(1)+(ik-1)*stepkx
	     yk=DG(2)+(ik-1)*stepky
	     zk=DG(3)+(ik-1)*stepkz
	 kcount=kcount+1
	 bandklist(1,kcount)=xk
	 bandklist(2,kcount)=yk
	 bandklist(3,kcount)=zk
	 end do

C K-->M
       stepkx=(DM(1)-DK(1))/(NSUBR2)
       stepky=(DM(2)-DK(2))/(NSUBR2)
       stepkz=(DM(3)-DK(3))/(NSUBR2)
       do  ik=1,NSUBR2
	     xk=DK(1)+(ik-1)*stepkx
	     yk=DK(2)+(ik-1)*stepky
	     zk=DK(3)+(ik-1)*stepkz
	 kcount=kcount+1
	 bandklist(1,kcount)=xk
	 bandklist(2,kcount)=yk
	 bandklist(3,kcount)=zk
	 end do

C M-->G
       stepkx=(DG(1)-DM(1))/(NSUBR3)
       stepky=(DG(2)-DM(2))/(NSUBR3)
       stepkz=(DG(3)-DM(3))/(NSUBR3)
       do  ik=1,NSUBR3+1
	     xk=DM(1)+(ik-1)*stepkx
	     yk=DM(2)+(ik-1)*stepky
	     zk=DM(3)+(ik-1)*stepkz
	 kcount=kcount+1
	 bandklist(1,kcount)=xk
	 bandklist(2,kcount)=yk
	 bandklist(3,kcount)=zk
	 end do


      write(*,*)"have obtained K-path successfully..."
	
	open(10,file='kpath.dat')
	write(10,*)"original_k_path:"
	write(10,"(3f12.6)")((bandklist(isit,knum),isit=1,3),
     $                                     knum=1,korigin)
	write(10,*)
      write(*,*)"have obtained new K-path successfully..."

	return
	end
c************************************************************************!!!
c************************************************************************!!!
c the subroutine used to diagonalize complex matrix !!!!!!!!!!!!!!!!!!!!!!!!
c************************************************************************!!!
c************************************************************************!!!
!!!!!!!!!!!!!!!!!!!!!use lapack!!!!!!!!!!!!!!!!!!!!!!!!      
      subroutine cal_eigenVS(Nsite,ham,EVAL,ch_state)
      IMPLICIT REAL(8) (a-h,o-z)
      PARAMETER ( LWMAX = 1000 )
      character(1) ch_state
      COMPLEX*16 ham(Nsite,Nsite),EVEC(Nsite,Nsite)
      dimension EVAL(Nsite)
      dimension RWORK( 3*Nsite-2 )
      COMPLEX(8) WORK( LWMAX ),c_check
!!!!!!!!!!!!!!!! use lapack !!!!!!!!!!!!!!!!!!!!!!!!!
!Compute all of the eigenvalues and eigenvectors of a complex Hermitian matrix.
!     Query the optimal workspace.
      if(ch_state.eq."V")then
       do i_site=1,Nsite
        do j_site=1,Nsite
         EVEC(j_site,i_site)=ham(j_site,i_site)
        end do
       end do
       LWORK = -1
       CALL ZHEEV('V','U',Nsite,ham,Nsite,EVAL,WORK,LWORK,RWORK,
     $            INFO)
       do i_site=1,Nsite
        do j_site=1,Nsite
         ham(j_site,i_site)=EVEC(j_site,i_site)
        end do
       end do
       LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
!      Solve eigenproblem.
       CALL ZHEEV('V','U',Nsite,ham,Nsite,EVAL,WORK,LWORK,RWORK,
     $            INFO )
!      Check for convergence.
       IF( INFO.GT.0 ) THEN
          WRITE(*,*)'The algorithm failed to compute eigenvalues.'
          STOP
       END IF
      else
       do i_site=1,Nsite
        do j_site=1,Nsite
         EVEC(j_site,i_site)=ham(j_site,i_site)
        end do
       end do
       LWORK = -1
       CALL ZHEEV('N','U',Nsite,ham,Nsite,EVAL,WORK,LWORK,RWORK,
     $            INFO)
       do i_site=1,Nsite
        do j_site=1,Nsite
         ham(j_site,i_site)=EVEC(j_site,i_site)
        end do
       end do
       LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
!      Solve eigenproblem.
       CALL ZHEEV('N','U',Nsite,ham,Nsite,EVAL,WORK,LWORK,RWORK,
     $            INFO )
!      Check for convergence.
       IF( INFO.GT.0 ) THEN
          WRITE(*,*)'The algorithm failed to compute eigenvalues.'
          STOP
       END IF
      end if
      return
      end 
c************************************************************************!!!
c************************************************************************!!!
      subroutine typetransform(kz,cfname,casename)
	implicit double precision(a-h,o-z)
	integer  iu,id,ih,it
	character(256) chunit,chdecd,chhund,cfname,casename

	         ih=kz/100

	         it=mod(kz,100)
	         id=it/10

	         it=mod(kz,100)
	         it=mod(it,10)
	         iu=it
	       
	         ih=ih+48
	         id=id+48
	         iu=iu+48

	         chhund=char(ih)
	         chdecd=char(id)
	         chunit=char(iu)
	         
	         cfname=trim(casename)

	         cfname=trim(cfname)//"_eigen"//trim(chhund)
     $                             //trim(chdecd) //trim(chunit)
	return
	end
c************************************************************************!!!
c************************************************************************!!!
