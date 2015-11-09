program main
    use module_common
    implicit none

    filename = "code.txt"
    open(10,file = filename,status = "old",iostat = ierror)
    filename = "code_c.txt"
    open(20,file = filename,status = "old",iostat = ierror)

    call cpu_time(t1)
    j = 0
    do 
        read(10,*,iostat = ierror)i
        if ( ierror /= 0 )then
            exit
        endif ! if ends
        j = j + i
    enddo !cycle ends
    print *,j 
    j = 0
    do 
        read(20,*,iostat = ierror)i
        if ( ierror /= 0 )then
            exit
        endif ! if ends
        j = j + i
    enddo !cycle ends
    print *,j 
    call cpu_time(t2)
    print *,"time cost is ", t2 - t1


end program main
