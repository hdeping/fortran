program main
    !use module_fst 
    use module_evolution
    integer  itmp

    real(8)        :: a(n)
    real(8)        :: b(n)
    real(8)        :: c(n)
  

    filename = "gr.txt"
    open(10,file = filename)
    filename = "cr.txt"
    open(20,file = filename)
    filename = "gk.txt"
    open(30,file = filename)
    filename = "ck.txt"
    open(40,file = filename)
    !filename = "test.txt"
    !open(50,file = filename)
    !  initial the k and r
    do i = 1,n
        dk(i) = (i - 1)*deltak
        dr(i) = (i - 1)*deltar
    end do

    maymm = may(dmm)
!********************************************************** 

    !  solve the equation
    call cpu_time(t1)
    call evolution()
    !call evolution_gamma_r()
    call cpu_time(t2)
    print *,"time cost is ",t2 - t1
    !************ check the result *************************
    !call check_cr()
    call check_gr()
    !******************************************************

!print g(r){{{
    print *,"lambda is ",lambda
    !  write ck to the file
    g_rmm(1) = 0.0

    do i = 2,n
        grmm(i) = grmm(i)/dr(i)
        grfm(i) = grfm(i)/dr(i)
        grff(i) = grff(i)/dr(i)
        crmm(i) = grmm(i)/dr(i)
        crfm(i) = grfm(i)/dr(i)
        crff(i) = grff(i)/dr(i)
        gkmm(i) = grmm(i)/dk(i)
        gkfm(i) = grfm(i)/dk(i)
        gkff(i) = grff(i)/dk(i)
        ckmm(i) = grmm(i)/dk(i)
        ckfm(i) = grfm(i)/dk(i)
        ckff(i) = grff(i)/dk(i)
    g_rmm(i)= crmm(i) + grmm(i) + 1
        g_rfm(i)= crfm(i)  +  grfm(i) + 1
        g_rff(i)= crff(i)  +  grff(i) + 1
    end do
    ! print the results
    do i = 2,n
        write(10,"(4f10.4)")dr(i),g_rmm(i),g_rfm(i),g_rff(i)
        write(20,"(4f10.4)")dr(i),grmm(i), grfm(i), grff(i)
        write(30,"(4f10.4)")dr(i),gkmm(i), gkfm(i), gkff(i)
        write(40,"(4f10.4)")dr(i),ckmm(i), ckfm(i), ckff(i)
    end do
    print *,"time cost is ",t2 - t1

    close(10)
    close(20)
    close(30)
    close(40)
!}}}
    
end program main
