Module Parameter
    real *8, dimension(:, :), allocatable :: kpointsRecip;
    real *8, dimension(:), allocatable :: hoppingPara, delta;
    integer, allocatable :: hoppingMatrix(:, :);
    integer atomNumber;
end Module Parameter

Program monolayerModel
    use Parameter
    implicit none       
    call readKpointsRecip();
    call readHoppingMatrix();
    call readHoppingParameter();
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
    integer readLength, numIdx;
    integer readStatus;
    atomNumber = 0;
    numIdx = 1;
    open(20, file = '../data/hoppingMat.dat');
    read(20, *) readLength;
    allocate(hoppingMatrix(6, readLength));
    do
        read (20, *, iostat = readStatus) hoppingMatrix(1, numIdx), hoppingMatrix(2, numIdx),&
         hoppingMatrix(3, numIdx), hoppingMatrix(4, numIdx), hoppingMatrix(5, numIdx), hoppingMatrix(6, numIdx);
        if (readStatus /= 0) then
            exit
        end if
        if (hoppingMatrix(6, numIdx) > atomNumber) then
            atomNumber =  hoppingMatrix(6, numIdx)
        end if
        numIdx = numIdx + 1;        
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
        write(*, *) "atom exceeded, change origin Code";
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

subroutine destructor()
    use Parameter
    deallocate(kpointsRecip);
    deallocate(hoppingMatrix);
    deallocate(delta);
end subroutine destructor