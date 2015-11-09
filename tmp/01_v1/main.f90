program main
    implicit none
    integer,parameter    :: n    = 48             ! months
    integer,parameter    :: s    = 24000          ! total money
    real,parameter       :: rate = 5E-3           ! rate
    real                 :: x                     ! money every month
    real                 :: y                     ! fluctuation
    real                 :: y1                    ! random number
    real                 :: s1                    ! 
    real                 :: s2                    !
    real                 :: tmp                   !
    integer              :: i                     !  count
    integer              :: j                     !
    integer              :: num(100)
    character(20)        :: filename

    filename = "new.txt"
    open(10,file = filename,status = "old")

    call random_seed()
    s1 = 4.056
    s2 = 89.056
    x  = s1/s2

    write(10,*)"s1 = ",s1 
    write(10,*)"s2 = ",s2
    write(10,*)"x  = ",x
    close(10)
    do i = 1,100
        num(i) = 100 + i
    end do
    do i = 1,100
        write(filename,"('file',i3.3,'.txt')")i
        open(num(i),file = filename,status = "old")
    end do
    do i  = 1,100
        write(num(i),"(i8)")i
    end do
    do i = 1,100
        close(num(i))
    end do
    
    
end program main
