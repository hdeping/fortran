program main
    !use module_fst 
    use module_evolution
    integer  itmp

    real(8)        :: a(n)
    real(8)        :: b(n)
    real(8)        :: c(n)
    

    filename = "g_r.txt"
    open(10,file = filename)
    filename = "hk.txt"
    open(20,file = filename)
    filename = "check.txt"
    open(40,file = filename)
    !filename = "hkmm.txt"
    !open(30,file = filename,status = "old",iostat = ierror)

    !filename = "g_rfm.txt"
    !open(20,file = filename)
    !filename = "g_rff.txt"
    !open(30,file = filename)
    !filename = "test.txt"
    !open(50,file = filename)
    !  initial the k and r
    do i = 1,n
        dk(i) = (i - 1)*deltak
        dr(i) = (i - 1)*deltar
    end do
    call random_seed()

    maymm = may(dmm)
    mayfm = may(dfm)
    mayff = may(dff)
    maysm = may(dsm)
    maysf = may(dsf)
    mayss = may(dss)
!********************************************************** 

    !  solve the equation
    times = 0
    call cpu_time(t1)
    ! get hkmm
    call evolution_mm()
    call evolution()
    !do i = 1,n
    !    write(20,"(4f18.9)")hkmm(i),hkfm(i),hkffc(i),hkffb(i)
    !end do
    !close(20)
    ! get single particle hk
    call evolution_particle()
    !call evolution_gr()
    call cpu_time(t2)
    print *,"time cost is ",t2 - t1

!print g(r){{{
    g_rmm(1) = 0.0

    do i = 2,n 
       g_rmm(i) = (crmm(i)  +  grmm(i) )/dr(i)  + 1.0  
       g_rfm(i) = (crfm(i)  +  grfm(i) )/dr(i)  + 1.0
       g_rff(i) = (crff(i)  +  grff(i) )/dr(i)  + 1.0
       g_rsm(i) = (crsm(i)  +  grsm(i) )/dr(i)  + 1.0
       g_rsf(i) = (crsf(i)  +  grsf(i) )/dr(i)  + 1.0
       g_rss(i) = (crss(i)  +  grss(i) )/dr(i)  + 1.0
    end do   !  i
    ! print the results
    !write(10,*)"r","g_rmm","g_rfm","g_rff"
    do i = 2,n
        write(10,"(7f10.4)")dr(i),g_rmm(i),g_rfm(i),g_rff(i),g_rsm(i),g_rsf(i),g_rss(i)
    end do
    print *,"time cost is ",t2 - t1

    close(10)
    close(20)
    close(30)
    close(40)
!}}}
    
end program main
    !do i = 1,n
    !    write(30,"(f18.9)")hkmm(i)
    !end do
    !close(30)
    
    ! get other hk
   
    !do i = 1,n
    !    read(30,*,iostat = ierror)hkmm(i)
    !    if(ierror /= 0)exit
    !end do
    !close(30)


