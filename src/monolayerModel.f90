Module Parameter
    integer atomNumber;
    real *8, dimension(:, :), allocatable :: kpointsRecip;
    real *8, dimension(:), allocatable :: hoppingPara, Delta;
    complex *16, allocatable :: Hkin(:, :);
    integer, allocatable :: hoppingMatrix(:, :, :, :);
    real *8, dimension(3) :: a1 = (/  3.3291,  -0.0000,  0.0000 /);
    real *8, dimension(3) :: a2 = (/ -1.6645,   2.8831, -0.0000 /);
    real *8, dimension(3) :: a3 = (/  0.0000,  -0.0000, 23.1180 /);
end module Parameter

Program monolayerModel
    use Parameter
    implicit none
    integer :: nsite1NumIdx = 0, nsite2NumIdx = 0, kpointsRecipNumIdx = 0, arrayNumIdx = 0;
    integer :: kpointsRecipLength = 0, hoppingMatrixLength = 0, readBaseLine = 1;
    integer, allocatable :: hoppingMatrixSplitSites(:, :, :);
    integer, dimension(3) :: hoppingNeighbourNumIdx = (/ 1, 1, 1 /);
    real *8, dimension(3, 3) :: transMatA;
    real *8, dimension(3) :: kpointsRecipSinglePoint;
    complex *16 :: KineticEnergy = 0D0;

    interface 
        subroutine calculateKineticEnergy(transMatA, kpointsRecipSinglePoint, hoppingPara, delta, hoppingMatrixSplitSites, Hkin)
        implicit none 
        real *8, dimension(3), intent(in) :: kpointsRecipSinglePoint;
        real *8, dimension(:), allocatable, intent(in) :: hoppingPara, delta;
        real *8, dimension(3, 3), intent(in) :: transMatA;
        integer, allocatable, intent(in) :: hoppingMatrixSplitSites(:, :, :);
        complex *16, intent(inout) :: Hkin;
        end subroutine calculateKineticEnergy
    end interface   

    call readKpointsRecip();
    call readHoppingMatrix();
    call readHoppingParameter();
    transMatA(:, 1) = a1(:);
    transMatA(:, 2) = a2(:);
    transMatA(:, 3) = a3(:);
    kpointsRecipLength = size(kpointsRecip, 1);
    do nsite2NumIdx = 1, atomNumber
        do nsite1NumIdx = 1, atomNumber           
            hoppingMatrixLength = size(hoppingMatrix(:, :, nsite1NumIdx, nsite2NumIdx), 2);
            if (.Not. allocated(hoppingMatrixSplitSites)) allocate(hoppingMatrixSplitSites(6, 12, size(hoppingPara)));
            hoppingMatrixSplitSites = 0D0;
            do readBaseLine = 1, hoppingMatrixLength
                if ( hoppingMatrix(1, readBaseLine, nsite1NumIdx, nsite2NumIdx) /= hoppingNeighbourNumIdx(1) ) then
                    hoppingNeighbourNumIdx(1) = hoppingNeighbourNumIdx(1) + 1;
                    hoppingNeighbourNumIdx(2) = 1;
                end if
                forall (arrayNumIdx = 1: 6)
                    hoppingMatrixSplitSites(arrayNumIdx, hoppingNeighbourNumIdx(2), hoppingNeighbourNumIdx(1)) &
                    = hoppingMatrix(arrayNumIdx, readBaseLine, nsite1NumIdx, nsite2NumIdx);
                end forall
                hoppingNeighbourNumIdx(2) = hoppingNeighbourNumIdx(2) + 1;
            end do
            do kpointsRecipNumIdx = 1, kpointsRecipLength
                kpointsRecipSinglePoint = kpointsRecip(:, kpointsRecipNumIdx);
                call calculateKineticEnergy(transMatA, kpointsRecipSinglePoint, hoppingPara, delta,&
                 hoppingMatrixSplitSites, KineticEnergy);
            end do            
        end do
    end do
    
    call destructor();
end Program monolayerModel

subroutine readKpointsRecip()
    use Parameter
    implicit none
    integer kpointsRecipLength, numIdx;
    open(10, file = '../data/kpointsRecip.dat');
    read(10, *) kpointsRecipLength;
    allocate(kpointsRecip(3, kpointsRecipLength));
    do numIdx = 1, kpointsRecipLength
        read(10, *) kpointsRecip(1, numIdx), kpointsRecip(2, numIdx), kpointsRecip(3, numIdx);
    end do
    close(10);
    
end subroutine readKpointsRecip  


subroutine readHoppingMatrix()
    use Parameter
    integer readLength, numIdx, nsite1, nsite2;
    integer readStatus;
    readStatus = 0;
    numIdx = 1;
    open(20, file = '../data/hoppingMat.dat');
    read(20, *) readLength, atomNumber;
    allocate(hoppingMatrix(6, readLength, atomNumber, atomNumber));
    do while (readStatus == 0)
        do nsite2 = 1, atomNumber
            do nsite1 = 1, atomNumber
                read (20, *, iostat = readStatus) hoppingMatrix(1, numIdx, nsite1, nsite2),&
                 hoppingMatrix(2, numIdx, nsite1, nsite2), hoppingMatrix(3, numIdx, nsite1, nsite2),&
                 hoppingMatrix(4, numIdx, nsite1, nsite2), hoppingMatrix(5, numIdx, nsite1, nsite2),&
                 hoppingMatrix(6, numIdx, nsite1, nsite2);
                numIdx = numIdx + 1;
            end do
        end do
    end do
    close(20)
end subroutine readHoppingMatrix

subroutine readHoppingParameter()
    use Parameter
    implicit none
    Integer readLine, numIdx;
    integer readStatus;
    open (30, file = '../data/hoppingPara.dat');
    read(30, *) readLine;
    allocate(hoppingPara(readLine - 1));
    allocate(delta(atomNumber));
    select case (atomNumber)
    case (1)
        read(30, *) delta(1);
    case (2)
        read(30, *) delta(1), delta(2);
    case (3)
        read(30, *) delta(1), delta(2), delta(3);
    case (4)
        read(30, *) delta(1), delta(2), delta(3), delta(4);
    case default
        write(*, *) "atom exceeded, change origin Code in subroutine readHoppingParameter";
        stop
    end select
    numIdx = 1;
    do
        read(30, *,iostat = readStatus) hoppingPara(numIdx);
        if (readStatus /= 0) then
            exit
        end if
        numIdx = numIdx + 1;
    end do

    close(30);
end subroutine readHoppingParameter

subroutine calculateKineticEnergy(transMatA, kpointsRecipSinglePoint, hoppingPara, delta, hoppingMatrixSplitSites, Hkin)
!   Declaration of passing variables
implicit none
real *8, dimension(3), intent(in) :: kpointsRecipSinglePoint;
real *8, dimension(:), allocatable, intent(in) :: hoppingPara, delta;
real *8, dimension(3, 3), intent(in) :: transMatA;
integer, allocatable, intent(in) :: hoppingMatrixSplitSites(:, :, :);
complex *16, intent(inout) :: Hkin;
!   Declaration of the temperary variables
integer :: neighbourLoopLength = 0, innerLoopLength = 0;
integer :: neighbourNumIdx = 0, innerNumIdx = 0, numIdx = 0;
complex *16 :: dotpart = 0D0, sum = 0D0;
complex *16, dimension(3) :: a1part = 0D0, a2part = 0D0, a3part = 0D0;
!   Subroutine Content
Hkin = 0D0;
neighbourLoopLength = size(hoppingPara);

do neighbourNumIdx = 1, neighbourLoopLength
    innerLoopLength = size(hoppingMatrixSplitSites, 2);
    do innerNumIdx = 1, innerLoopLength
        forall (numIdx = 1: 3)
            a1part(numIdx) = hoppingMatrixSplitSites(neighbourNumIdx, innerNumIdx, numIdx) * transMatA(numIdx, 1);
            a2part(numIdx) = hoppingMatrixSplitSites(neighbourNumIdx, innerNumIdx, numIdx) * transMatA(numIdx, 2);
            a3part(numIdx) = hoppingMatrixSplitSites(neighbourNumIdx, innerNumIdx, numIdx) * transMatA(numIdx, 3);
        end forall
        dotpart = dot_product(kpointsRecipSinglePoint, a1part + a2part + a3part);
    end do
end do

write(*, *) kpointsRecipSinglePoint, hoppingPara(1), delta(1), hoppingMatrixSplitSites(1, 1, 1), Hkin;

end subroutine calculateKineticEnergy

subroutine destructor()
    use Parameter
    deallocate(kpointsRecip);
    deallocate(hoppingMatrix);
    deallocate(delta);
    deallocate(Hkin);
end subroutine destructor