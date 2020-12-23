Module Parameter
    real *8, dimension(:, :), allocatable :: kpointsRecip;
    integer, allocatable :: hoppingMatrix(:, :);
end Module Parameter

Program monolayerModel
    use Parameter
    implicit none    
    call readKpointsRecip();
    call readHoppingMatrix();
    Call destructor();
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
    integer readLength, numIdx
    integer readStatus;
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
        numIdx = numIdx + 1;        
    end do
    close(20)
end subroutine readHoppingMatrix

subroutine destructor()
    use Parameter
    deallocate(kpointsRecip);
    deallocate(hoppingMatrix);
end subroutine destructor