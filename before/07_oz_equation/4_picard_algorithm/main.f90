program main
    !use module_fst 
    use module_evolution
    integer  itmp

    real(8)        :: a(n)
    real(8)        :: b(n)
    real(8)        :: c(n)
    

    filename = "g_r.txt"
    open(10,file = filename)
    filename = "sk.txt"
    open(20,file = filename)
    filename = "check.txt"
    open(40,file = filename)

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

    maymm = may(dmm)
    mayfm = may(dfm)
    mayff = may(dff)
!********************************************************** 

    !  solve the equation
    call cpu_time(t1)
    ! get hkmm
    call evolution_mm()
    call evolution()
    !call evolution_gr()
    call cpu_time(t2)
    print *,"time cost is ",t2 - t1
    !************ check the result *************************
    !call check_cr()
    !call check_gr()
    !******************************************************
    !filename = "tmp_data.txt"
    !open(50,file = filename)
    !do i = 1,n
    !    write(50,"(5f15.6)")dk(i),gkmm(i),gkfm(i),gkff(i),gkffb(i)
    !end do
    !close(50)

!print g(r){{{
    print *,"lambda is ",lambda
    !  write ck to the file
    g_rmm(1) = 0.0


    do i = 2,n 
        g_rmm(i) = (crmm(i)  +  grmm(i) )/dr(i)  + 1.0 
        g_rfm(i) = (crfm(i)  +  grfm(i) )/dr(i)  + 1.0
        g_rff(i) = (crff(i)  +  grff(i) )/dr(i)  + 1.0
        skmm(i)  = trho*xrate1**2.0*hkmm(i)/dk(i)  + 1.0 
        skfm(i)  = trho*xrate1*xrate2*(ckfm(i)  +  gkfm(i) )/dk(i)
        skff(i)  = trho*xrate2**2.0*(ckff(i)  +  gkff(i) )/dk(i)  + 1.0 
    end do   !  i
    ! print the results
    write(10,*)"r","g_rmm","g_rfm","g_rff"
    do i = 2,n
        write(10,"(4f18.8)")dr(i),g_rmm(i),g_rfm(i),g_rff(i) 
        write(20,"(4f18.8)")dk(i),skmm(i),skfm(i),skff(i) 
    end do
    print *,"time cost is ",t2 - t1

    close(10)
    close(20)
!}}}
    
end program main
