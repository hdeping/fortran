program main
    use module_common
    implicit none
    integer                icount
    real(8)                vvalue
    
    filename = "data.txt"
    open(10,file = filename)
    call random_seed()
    !do j = 1,50
        j   = 10
        eta = 1D-1 * dble(j)
        !eta = 0D0
        ! initialize the configuration
        call init_config()
        icount = 0
        vvalue = 0D0
        veloMeanTime = 0D0
        call cpu_time(t1)
        do i = 1,n
            call update()
            if ( mod(i,fre) == 0 )then
                icount = icount + 1
                vMean  = getv(velo)
                vvalue = (vvalue*(icount - 1) + vMean)/dble(icount)
                do k = 1,num
                    do kk = 1,2
                        veloMeanTime(k,kk) = (veloMeanTime(k,kk)*&
                                (icount - 1) + velo(k,kk))/dble(icount)
                    end do
                end do
                print *,"mean velocity is ",vvalue
                pause
            endif ! if ends
        enddo !cycle ends
        call getNormVelo(veloMeanTime)
        vMean = getv(veloMeanTime)
        write(10,*)eta,vvalue,vMean
        !write(10,*)eta,vMean
        call cpu_time(t2)
        print *,"time cost is ",t2 - t1
    !enddo !cycle ends
    close(10)
end program main
