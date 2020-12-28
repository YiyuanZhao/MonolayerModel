Module Parameter
    integer atomNumber;
    real *8, dimension(:, :), allocatable :: kpointsRecip;
    real *8, dimension(:), allocatable :: hoppingPara, delta;
    complex *16, allocatable :: Hkin(:, :);
    integer, allocatable :: hoppingMatrix(:, :, :, :);   
    real *8, dimension(3) :: a1 = (/  3.3291,  -0.0000,  0.0000 /);
    real *8, dimension(3) :: a2 = (/ -1.6645,   2.8831, -0.0000 /);
    real *8, dimension(3) :: a3 = (/  0.0000,  -0.0000, 23.1180 /);
end module Parameter

Program monolayerModel
    use Parameter
    implicit none
    real *8, dimension(3, 3) :: transMatA;  

    interface 
        subroutine calculateKineticEnergy(transMatA, kpointsRecip, hoppingPara, delta, hoppingMatrix, atomNumber, Hkin)
        implicit none 
        real *8, dimension(:, :), allocatable, intent(in) :: kpointsRecip;
        real *8, dimension(:), allocatable, intent(in) :: hoppingPara, delta;
        real *8, dimension(3, 3), intent(in) :: transMatA;
        integer, intent(in) :: atomNumber;
        integer, allocatable, intent(in) :: hoppingMatrix(:, :, :, :);
        complex *16, allocatable, intent(inout) :: Hkin(:, :);
        end subroutine calculateKineticEnergy
    end interface   

    call readKpointsRecip();
    call readHoppingMatrix();
    call readHoppingParameter();
    transMatA(1, :) = a1(:);
    transMatA(2, :) = a2(:);
    transMatA(3, :) = a3(:);
    call calculateKineticEnergy(transMatA, kpointsRecip, hoppingPara, delta, hoppingMatrix, atomNumber, Hkin);
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
    allocate(hoppingMatrix(atomNumber, atomNumber, 6, readLength));
    do while (readStatus == 0)
        do nsite2 = 1, atomNumber
            do nsite1 = 1, atomNumber
                read (20, *, iostat = readStatus) hoppingMatrix(nsite1, nsite2, 1, numIdx),&
                 hoppingMatrix(nsite1, nsite2, 2, numIdx), hoppingMatrix(nsite1, nsite2, 3, numIdx),&
                 hoppingMatrix(nsite1, nsite2, 4, numIdx), hoppingMatrix(nsite1, nsite2, 5, numIdx),&
                 hoppingMatrix(nsite1, nsite2, 6, numIdx);
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

subroutine calculateKineticEnergy(transMatA, kpointsRecip, hoppingPara, delta, hoppingMatrix, atomNumber, Hkin)
implicit none
real *8, dimension(:, :), allocatable, intent(in) :: kpointsRecip;
real *8, dimension(:), allocatable, intent(in) :: hoppingPara, delta;
real *8, dimension(3, 3), intent(in) :: transMatA;
integer, intent(in) :: atomNumber;
integer, allocatable, intent(in) :: hoppingMatrix(:, :, :, :);
complex *16, allocatable, intent(inout) :: Hkin(:, :);

allocate(Hkin(atomNumber, atomNumber));
Hkin(1, 1) = 1;
write(*, *) kpointsRecip(1, 1), hoppingPara(1), delta(1), hoppingMatrix(1,1,1,1), Hkin(1, 1);

end subroutine calculateKineticEnergy

subroutine destructor()
    use Parameter
    deallocate(kpointsRecip);
    deallocate(hoppingMatrix);
    deallocate(delta);
    deallocate(Hkin);
end subroutine destructor