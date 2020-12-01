!// ForQuill v1.01 Beta www.fcode.cn
!************************************************************************!!!
! the main parameters used in the calculations !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!************************************************************************!!!
Module param
    ! Implicit Double Precision (A-H, O-Z)
    Integer nwann, nhalfwan, nrdg
    Complex *16, Allocatable :: ham_r(:, :, :)
    Integer, Allocatable :: rnspac(:, :)
    Character (256) filename, casename
    Integer korigin, kwien2k, korigextend, kwork
    Real *8, Allocatable :: bandklist(:, :)
    Real *8, Allocatable :: bandklist_extend(:, :)
  ! the fractional translation vector !!!
    Real *8 vec_trans(3), unfoldvec(3)
  ! the transformation matrix in orbital space !!!
    Real *8, Allocatable :: w_ci(:, :)
  ! band structure in full brillouin zone !!!
  !	integer   NKX,NKY,NKZ
    Integer nkx, nky, dnky, tnky
    Parameter (nkx=151, nky=75, dnky=151, tnky=dnky+nky)
    Real *8 pi, twopi
    Real *8 veclat(3)
  End Module param
  !************************************************************************!!!
  !************************************************************************!!!
  Program get_hop_matrix
    Use param
    Implicit None
  ! Timer Settings
    real(kind=4)    :: t1,t2
    character(len=10)   :: time1,time2
    character(len=8)    :: thedate
    call date_and_time(thedate,time1)

  !  Timer Settings End
    filename = '../data/grap.dat'
  ! constants !!!
    pi = dasin(1.0D0)*2.0D0
    twopi = 2.0D0*pi
  
  ! to get hopping matrix t_mn(R) !!!
    Call read_hr()
    Call generat_k()
  
  ! to get K_path !!!
  
    kwork = korigin
    Call cal_origin_band(bandklist, kwork)
  !
    Call cal_hk_allzone()
    ! Timer Settings Begin
    call date_and_time(thedate,time2)
    read(time1,*) t1
    read(time2,*) t2
    write(*,*) 'time:', t2-t1
    ! Timer Settings End
  End Program get_hop_matrix
  !************************************************************************!!!
  !************************************************************************!!!
  ! to read hamitonial in wannier representations !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !************************************************************************!!!
  !************************************************************************!!!
  Subroutine read_hr()
    Use param
    ! Implicit Double Precision (A-H, O-Z)
    Character (30) ch
    Integer, Allocatable :: ndegen(:)
    Open (10, File=filename)
    Read (10, *) ch
    Read (10, *) nwann
    Write (*, *) 'nwann=', nwann
    nhalfwan = nwann/2
    Read (10, *) nrdg
    Write (*, *) 'NRdg=', nrdg
  ! to read out degenerate of every real point!
    If (.Not. allocated(ndegen)) Allocate (ndegen(nrdg))
    Read (10, '(15i5)')(ndegen(isit), isit=1, nrdg)
    If (allocated(ndegen)) Deallocate (ndegen)
    If (.Not. allocated(rnspac)) Allocate (rnspac(3,nrdg))
    If (.Not. allocated(ham_r)) Allocate (ham_r(nwann,nwann,nrdg))
  ! to read out hamitonial in real space !
    Do isit = 1, nrdg
      Do iorb = 1, nwann
        Do jorb = 1, nwann
          Read (10, '(5i5,2f12.6)') nrx, nry, nrz, korb, lorb, ham_r(jorb, iorb, isit)
          If (dabs(dimag(ham_r(jorb,iorb,isit)))<2.D-6) Then
            ham_r(jorb, iorb, isit) = dcmplx(dble(ham_r(jorb,iorb,isit)), 0.0D0)
          End If
        End Do
      End Do
      rnspac(1, isit) = nrx
      rnspac(2, isit) = nry
      rnspac(3, isit) = nrz
    End Do
    Close (10)
    Write (*, *) 'have read case_hr.dat successfully...'
    Return
  End Subroutine read_hr
  !************************************************************************!!!
  !************************************************************************!!!
  !************************************************************************!!!
  !************************************************************************!!!
  Subroutine cal_origin_band(bandkpath, ktot)
    Use param
    ! Implicit Double Precision (A-H, O-Z)
    Real *8, Allocatable :: eigenval(:, :)
    Real *8 eval(nwann), bandkpath(3, ktot)
    Complex *16 ham_work(nwann, nwann), hak(nwann, nwann)
    Real *8 wk(3)
  
  
    If (.Not. allocated(eigenval)) Allocate (eigenval(nwann,ktot))
    ktotal = ktot
    Do kloop = 1, ktotal
      Do korb = 1, nwann
        eigenval(korb, kloop) = 0.0D0
      End Do
    End Do
  
    Do kloop = 1, ktotal
  
      Do ik = 1, 3
        wk(ik) = bandkpath(ik, kloop)
      End Do
  
      Call cal_spectrum(wk, hak, ham_work, eval, nwann, ham_r, rnspac, nrdg)
  
      Do iband = 1, nwann
        eigenval(iband, kloop) = eval(iband)
      End Do
  
    End Do
  
    Open (10, File='../data/spectrum_orig.dat')
    Do iband = 1, nwann
      Do kloop = 1, ktotal
        Write (10, '(i5,f12.6)') kloop, eigenval(iband, kloop)
      End Do
      Write (10, *)
    End Do
    Close (10)
    Return
  End Subroutine cal_origin_band
  !**********************************************************************************!!!
  !**********************************************************************************!!!
  Subroutine cal_hk_allzone()
    Use param
    ! Implicit Double Precision (A-H, O-Z)
    Real *8 wk(3), xk, yk, a, dkx, dky, eval(nwann)
    Complex *16, Allocatable :: hk(:, :, :, :)
    Complex *16 ham_work(nwann, nwann), hak(nwann, nwann)
    a = 2.46D0/sqrt(3.0D0) !A value: relaxed C-C or B
    pi = dasin(1.0D0)*2.0D0
  !********one sixth Brilliouin zone****************************************
    
    
  !********one sixth Brilliouin zone****************************************
    If (.Not. allocated(hk)) Allocate (hk(nwann,nwann,1:tnky,1:(nkx+1)))
    Open (11, File='../data/HK.dat')
    Open (17, File='../data/kx+ky.dat')
    nksum = 0
    Do kx = 1, nkx + 1
      dkx = (2.0D0*pi/(3.0*a))/dble(nkx-1)
      xk = dble(kx-1)*dkx
      Do ky = 1, tnky
        ndy = ky - 2*nky - 1
        dky = (2.0D0*pi/(3.0*sqrt(3.0D0)*a))/dble(nky)
        yk = dble(ky-2.0D0*nky-1)*dky
        nksum = nksum + 1
        wk(1) = xk
        wk(2) = yk
        Write (17, '(2f12.6)') xk, yk
        wk(3) = 0.0D0
        Call cal_spectrum(wk, hak, ham_work, eval, nwann, ham_r, rnspac, nrdg)
        Write (17, '(2f12.6)') xk, yk
        Do iorb = 1, nwann
          Do jorb = 1, nwann
            hk(jorb, iorb, ky, kx) = hak(jorb, iorb)
            Write (11, '(2f12.6)') hk(jorb, iorb, ky, kx)
          End Do
        End Do
      End Do
    End Do
    Write (11, *) 'The-End!'
    Write (*, *) 'total k points', nksum
    Close (17)
    Close (11)
    Return
  End Subroutine cal_hk_allzone
  !********all Brilliouin zone****************************************
  !**********************************************************************************!!!
  !**********************************************************************************!!!
  !to calculate bandstructure of given point based on unit cell used in the calclation !
  !**********************************************************************************!!!
  !**********************************************************************************!!!
  Subroutine cal_spectrum(wk, hak, hamvec, eig, nwann, ham_r, rnspac, nrdg)
    ! Implicit Double Precision (A-H, O-Z)
    Real *8 wk(3), rn(3), eig(nwann)
    Complex *16 structphase
    Complex *16 hamk_wk(nwann, nwann), ham_work(nwann, nwann)
    Complex *16 hak(nwann, nwann)
    Complex *16 hamvec(nwann, nwann), ham_r(nwann, nwann, nrdg)
    Real *8 a1(3), a2(3), a3(3)
    Real *8 b1(3), b2(3), b3(3)
    Real *8 temp(3), newk(3), pi, twopi
    Integer rnspac(3, nrdg)
    Character (30) ch_state
    Double precision dotproduct
  
    pi = dasin(1.0D0)*2.0D0
    twopi = 2.0D0*pi
  !     relaxed lattice constant need to modify!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    a1(1) = 2.1304225D0
    a1(2) = -1.23D0
    a1(3) = 0.0D0
  
    a2(1) = 0.0D0
    a2(2) = 2.46D0
    a2(3) = 0.0D0
  
    a3(1) = 0.0D0
    a3(2) = 0.0D0
    a3(3) = 20.3D0
  
    b1(1) = 0.469390463D0*twopi
    b1(2) = 0.0D0
    b1(3) = 0.0D0
  
    b2(1) = 0.234695231D0*twopi
    b2(2) = 0.406504065D0*twopi
    b2(3) = 0.0D0
  
    b3(1) = 0.0D0
    b3(2) = 0.0D0
    b3(3) = 0.0492610837D0*twopi
  
  
    Do iorb = 1, nwann
      Do jorb = 1, nwann
        hamk_wk(jorb, iorb) = dcmplx(0.0D0, 0.0D0)
      End Do
    End Do
  
    Do isit = 1, nrdg
  ! for the translation of crystal based used unit cell !!!
  
  ! for the position of nearest cell !!!
      Do jsit = 1, 3
        temp(jsit) = 0.0D0
      End Do
  
  ! for A1
      Do jsit = 1, 3
        temp(jsit) = temp(jsit) + dble(rnspac(1,isit))*a1(jsit)
      End Do
  ! for A2
      Do jsit = 1, 3
        temp(jsit) = temp(jsit) + dble(rnspac(2,isit))*a2(jsit)
      End Do
  ! for A3
      Do jsit = 1, 3
        temp(jsit) = temp(jsit) + dble(rnspac(3,isit))*a3(jsit)
      End Do
  
      Do jsit = 1, 3
        rn(jsit) = temp(jsit)
      End Do
  
  ! for k vectors !
      Do jsit = 1, 3
        temp(jsit) = 0.0D0
      End Do
  !
  ! for B1
  !					 do jsit=1,3
  !						temp(jsit)=temp(jsit)+
  !     $                               dble(wk(1))*B1(jsit)
  !					 end do
  ! for B2
  !					 do jsit=1,3
  !						temp(jsit)=temp(jsit)+
  !     $                               dble(wk(2))*B2(jsit)
  !					 end do
  ! for B3
  !					 do jsit=1,3
  !						temp(jsit)=temp(jsit)+
  !     $                               dble(wk(3))*B3(jsit)
  !					 end do
  
      Do jsit = 1, 3
        newk(jsit) = temp(jsit) + wk(jsit)
      End Do
  
      dotproduct = 0.0D0
      Do jsit = 1, 3
        dotproduct = dotproduct + newk(jsit)*rn(jsit)
      End Do
  
  ! the K vector using 2pi/a , 2pi/b and 2pi/c as unit !
  ! the R vector using a , b and c as unit !!!!!!!!!!!!!
  
      ! dotproduct = dotproduct
      structphase = dcmplx(dcos(dotproduct), dsin(dotproduct))
  
  ! for original part alpha !
  ! to take the original phase factor !
  
      Do iorb = 1, nwann
        Do jorb = 1, nwann
          hamk_wk(jorb, iorb) = hamk_wk(jorb, iorb) + ham_r(jorb, iorb, isit)*structphase
        End Do
      End Do
  
  ! *******************************!
  ! have obtained H_wk !!!!!!!!!!!!!
  ! *******************************!
    End Do
  
    Do iorb = 1, nwann
      If (dabs(dimag(hamk_wk(iorb,iorb)))<1.D-15) Then
        hamk_wk(iorb, iorb) = dcmplx(dble(hamk_wk(iorb,iorb)), 0.0D0)
      End If
    End Do
  ! to get ham_work !!!

    ! ham_work = dcmplx(0.0D0, 0.0D0)
    ! ! Do iorb = 1, nwann
    ! !   Do jorb = 1, nwann
    ! !     ham_work(jorb, iorb) = dcmplx(0.0D0, 0.0D0)
    ! !   End Do
    ! ! End Do
  
    !  Do iorb = 1, nwann
    !   Do jorb = 1, nwann
    !     ham_work(jorb, iorb) = ham_work(jorb, iorb) + hamk_wk(jorb, iorb)
    !   End Do
    ! End Do
    ham_work = hamk_wk;
    hak = ham_work;

    ! Do iorb = 1, nwann
    !   Do jorb = 1, nwann
    !     hak(jorb, iorb) = ham_work(jorb, iorb)
    !   End Do
    ! End Do
  ! to diagonalize hamk_work !
    ch_state = 'V'
    Call cal_eigenvs(nwann, ham_work, eig, ch_state)

    hamvec = ham_work;
    ! Do lorb = 1, nwann
    !   Do korb = 1, nwann
    !     hamvec(korb, lorb) = ham_work(korb, lorb)
    !   End Do
    ! End Do
  
    Return
  End Subroutine cal_spectrum
  !************************************************************************!!!
  !************************************************************************!!!
  ! to calculate band structure with given K-path !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !************************************************************************!!!
  !************************************************************************!!!
  
  !************************************************************************!!!
  !************************************************************************!!!
  ! to calculate band structure with given K-path !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !************************************************************************!!!
  !************************************************************************!!!
  
  !************************************************************************!!!
  !************************************************************************!!!
  !************************************************************************!!!
  Subroutine generat_k()
    Use param
    ! Implicit Double Precision (A-H, O-Z)
  !	parameter(A=1.43590)
    Parameter (nsubr1=85)
    Parameter (nsubr2=42)
    Parameter (nsubr3=73)
  ! original G point
    Real *8 dg(3)
  ! original X point
    Real *8 dk(3)
    Real *8 dm(3)
    Real *8 a, stepkx, stepky, stepkz, xk, yk, zk
    a = 2.46D0/sqrt(3.0D0) !A value: relaxed C-C or 
    pi = dasin(1.0D0)*2.0D0
  
  ! to determine korigin
  
    dg(1) = 0.00D0
    dg(2) = 0.00D0
    dg(3) = 0.00D0
  
    dk(1) = 2.0D0*pi/(3.0D0*a)
    dk(2) = 2.0D0*pi/3.0D0/sqrt(3.0D0)/a
    dk(3) = 0.00D0
  
    dm(1) = 2.0D0*pi/(3.0D0*a)
    dm(2) = 0.00D0
    dm(3) = 0.00D0
  
  
    korigin = nsubr1 + nsubr2 + nsubr3 + 1
    If (.Not. allocated(bandklist)) Allocate (bandklist(3,korigin))
  
  
    kcount = 0
    stepkx = (dk(1)-dg(1))/(nsubr1)
    stepky = (dk(2)-dg(2))/(nsubr1)
    stepkz = (dk(3)-dg(3))/(nsubr1)
  ! G-->K
    Do ik = 1, nsubr1
      xk = dg(1) + (ik-1)*stepkx
      yk = dg(2) + (ik-1)*stepky
      zk = dg(3) + (ik-1)*stepkz
      kcount = kcount + 1
      bandklist(1, kcount) = xk
      bandklist(2, kcount) = yk
      bandklist(3, kcount) = zk
    End Do
  
  ! K-->M
    stepkx = (dm(1)-dk(1))/(nsubr2)
    stepky = (dm(2)-dk(2))/(nsubr2)
    stepkz = (dm(3)-dk(3))/(nsubr2)
    Do ik = 1, nsubr2
      xk = dk(1) + (ik-1)*stepkx
      yk = dk(2) + (ik-1)*stepky
      zk = dk(3) + (ik-1)*stepkz
      kcount = kcount + 1
      bandklist(1, kcount) = xk
      bandklist(2, kcount) = yk
      bandklist(3, kcount) = zk
    End Do
  
  ! M-->G
    stepkx = (dg(1)-dm(1))/(nsubr3)
    stepky = (dg(2)-dm(2))/(nsubr3)
    stepkz = (dg(3)-dm(3))/(nsubr3)
    Do ik = 1, nsubr3 + 1
      xk = dm(1) + (ik-1)*stepkx
      yk = dm(2) + (ik-1)*stepky
      zk = dm(3) + (ik-1)*stepkz
      kcount = kcount + 1
      bandklist(1, kcount) = xk
      bandklist(2, kcount) = yk
      bandklist(3, kcount) = zk
    End Do
  
  
    Write (*, *) 'have obtained K-path successfully...'
  
    Open (10, File='../data/kpath.dat')
    Write (10, *) 'original_k_path:'
    Write (10, '(3f12.6)')((bandklist(isit,knum),isit=1,3), knum=1, korigin)
    Write (10, *)
    Write (*, *) 'have obtained new K-path successfully...'
  
    Return
  End Subroutine generat_k
  !************************************************************************!!!
  !************************************************************************!!!
  ! the subroutine used to diagonalize complex matrix !!!!!!!!!!!!!!!!!!!!!!!!
  !************************************************************************!!!
  !************************************************************************!!!
  !!!!!!!!!!!!!!!!!!!!!use lapack!!!!!!!!!!!!!!!!!!!!!!!!
  Subroutine cal_eigenvs(nsite, ham, eval, ch_state)
    ! Implicit Real (8)(A-H, O-Z)
    Parameter (lwmax=1000)
    Character (1) ch_state
    Character UPLO
    Complex *16 ham(nsite, nsite), evec(nsite, nsite)
    ! Dimension eval(nsite)
    Real *8 rwork(3*nsite-2)
    Complex (Kind=8) work(lwmax)
    Real *8 eval(nsite)
    Integer lwork, NB
  !!!!!!!!!!!!!!!! use lapack !!!!!!!!!!!!!!!!!!!!!!!!!
  !Compute all of the eigenvalues and eigenvectors of a complex Hermitian matrix.
  !     Query the optimal workspace.
    If (ch_state=='V') Then
      ! evec = ham;
      ! Do i_site = 1, nsite
      !   Do j_site = 1, nsite
      !     evec(j_site, i_site) = ham(j_site, i_site)
      !   End Do
      ! End Do
      ! lwork = -1
      ! Call zheev('V', 'U', nsite, ham, nsite, eval, work, lwork, rwork, info)
      ! ham = evec;
      ! Do i_site = 1, nsite
      !   Do j_site = 1, nsite
      !     ham(j_site, i_site) = evec(j_site, i_site)
      !   End Do
      ! End Do
      NB = ILAENV( 1, 'ZHETRD', UPLO, nsite, -1, -1, -1 )
      LWKOPT = MAX( 1, ( NB+1 )*nsite )
      WORK( 1 ) = LWKOPT
      lwork = min(lwmax, int(work(1)))
  !      Solve eigenproblem.
      Call zheev('V', 'U', nsite, ham, nsite, eval, work, lwork, rwork, info)
  !      Check for convergence.
      If (info>0) Then
        Write (*, *) 'The algorithm failed to compute eigenvalues.'
        Stop
      End If
    Else
      evec = ham
      Do i_site = 1, nsite
        Do j_site = 1, nsite
          evec(j_site, i_site) = ham(j_site, i_site)
        End Do
      End Do
      lwork = -1
      Call zheev('N', 'U', nsite, ham, nsite, eval, work, lwork, rwork, info)
      Do i_site = 1, nsite
        Do j_site = 1, nsite
          ham(j_site, i_site) = evec(j_site, i_site)
        End Do
      End Do
      lwork = min(lwmax, int(work(1)))
  !      Solve eigenproblem.
      Call zheev('N', 'U', nsite, ham, nsite, eval, work, lwork, rwork, info)
  !      Check for convergence.
      If (info>0) Then
        Write (*, *) 'The algorithm failed to compute eigenvalues.'
        Stop
      End If
    End If
    Return
  End Subroutine cal_eigenvs
  !************************************************************************!!!
  !************************************************************************!!!
  Subroutine typetransform(kz, cfname, casename)
    ! Implicit Double Precision (A-H, O-Z)
    Integer iu, id, ih, it
    Character (256) chunit, chdecd, chhund, cfname, casename
  
    ih = kz/100
  
    it = mod(kz, 100)
    id = it/10
  
    it = mod(kz, 100)
    it = mod(it, 10)
    iu = it
  
    ih = ih + 48
    id = id + 48
    iu = iu + 48
  
    chhund = char(ih)
    chdecd = char(id)
    chunit = char(iu)
  
    cfname = trim(casename)
  
    cfname = trim(cfname) // '_eigen' // trim(chhund) // trim(chdecd) // trim(chunit)
    Return
  End Subroutine typetransform
  !************************************************************************!!!
  !************************************************************************!!!
  