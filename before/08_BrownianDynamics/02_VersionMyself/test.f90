program test
    use module_common
    implicit none

    integer,parameter     :: arrayNum = 10
    integer               :: num(arrayNum,m)
    real(8)               :: rangeVelocity(arrayNum)
    real(8)               :: t1 ! time recording
    real(8)               :: t2 ! time recording
    real(8)               :: tmp
    real(8)               :: xtmp

    call getInitialParameter()
    do i = 1,10
        print *,"i = ",i
    enddo !cycle ends
     
    filename = "data.txt"
    open(10,file = filename)
    !test the gauss random number
    
    do i = 1,arrayNum
        rangeVelocity(i) = (i*2.0 - arrayNum)*veloMax/dble(arrayNum)
    enddo !cycle ends
    num  = 0
    do i = 1,n
        do k = 1,m
            do j = 1,arrayNum
                if ( velocity(i,k) < rangeVelocity(j) )then
                    num(j,k) = num(j,k) + 1
                    exit
                endif ! if ends
            enddo !cycle ends
        enddo !cycle ends
    enddo !cycle ends


    ! write the data
    do i = 1,arrayNum
        write(10,*)i,num(i,:)
    enddo !cycle ends
    print *,sum(velocity(:,1)),sum(velocity(:,2))
    close(10)

    call cpu_time(t1)
    do i = 1,n
        call random_number(xtmp)
    enddo !cycle ends
    call cpu_time(t2)
    print *,"random_number time cost is ",t2 - t1
    call cpu_time(t1)
    xtmp = 0.3
    do i = 1,n
        tmp = sin(xtmp)
    enddo !cycle ends
    call cpu_time(t2)
    print *,"random_number time cost is ",t2 - t1
     

end program test
