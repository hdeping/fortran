program main
    use module_common
    implicit none
    integer                icount
    real(8)                vvalue
    
    filename = "data.txt"
    open(10,file = filename)
    call random_seed()
    ! initialize the configuration
    call init_config()
    eta = 2D1
    do j = 1,100
        eta = 1D-1 * dble(j)
        icount = 0
        vvalue = 0D0
        call cpu_time(t1)
        do i = 1,n
            call update
            if ( mod(i,fre) == 0 )then
                icount = icount + 1
                vMean  = getv(angle,num)
                vvalue = (vvalue*(icount - 1) + vMean)/dble(icount)
                !print *,"mean velocity is ",vvalue
                !pause
            endif ! if ends
        enddo !cycle ends
        write(10,*)eta,vvalue
        call cpu_time(t2)
        print *,"time cost is ",t2 - t1
    enddo !cycle ends
    close(10)
end program main
