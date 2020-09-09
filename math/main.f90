program main
    use module_common
    implicit none
    real(8)               :: a
    real(8)               :: b
    real(8)               :: length

    filename = "data.txt"
    open(10,file = filename)

    !a = 1.0
    !b = 1.0
    !length = getlength(a,b)
    !print *,length*2.0
    do i = 1,n
        !call cpu_time(t1)
        do j = 1,n
            a      = 2D-2*dble(i)
            b      = 2D-2*dble(j)
            length = getlength(a,b)
            write(10,*)a,b,length
        enddo !cycle ends
        !call cpu_time(t2)
        !print *,"time cost is ",t2 - t1
    enddo !cycle ends
    

    close(10)
end program main
