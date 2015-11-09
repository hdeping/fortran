program main
    use module_common
    implicit none
    
    table = getMultiplyTable()
    
    do i = 1,9
        do j = 1,9
            print *,i,j,table(i,j)
        end do
    end do
    pause
    filename = "data.dat"
    open(10,file = filename)
    filename = "result.dat"
    open(20,file = filename)
    ! initial the a(n) and b(m)
!initial a,b{{{
    call random_seed()
    ! get a(1) (not equal to 0)
    do 
       call random_number(x1)
       stmp = int(10*x1)
       if (stmp /= 0)exit
    enddo
    a(1) = stmp
    ! get b(1) (not equal to 0)
    do 
       call random_number(x1)
       stmp = int(10*x1)
       if (stmp /= 0)exit
    enddo
    b(1) = stmp
    ! get the rest values of a 
    do i = 2,n
       call random_number(x1)
       a(i) = int(10*x1)
    end do
    ! get the rest values of b
    do i = 2,m
       call random_number(x1)
       b(i) = int(10*x1)
    end do
!}}}
    ! get the results
    ! write into files
    write(10,"(<n>i1,'*',<m>i1)")a,b
    c = multiply(a,b,n,m)
    write(20,"(68i1)"),c 

    close(10)
end program main
