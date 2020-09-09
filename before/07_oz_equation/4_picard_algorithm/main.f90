program main
    !use module_fst 
    use module_evolution
    integer  itmp

    real(8)        :: a(n)
    real(8)        :: b(n)
    real(8)        :: c(n)
    character(7)   :: ch_b = "       "
    

    filename = "g_r.txt"
    open(10,file = filename)
    filename = "gr.txt"
    open(20,file = filename)
    filename = "cr.txt"
    open(30,file = filename)
    filename = "hr.txt"
    open(40,file = filename)
    filename = "check.txt"
    open(50,file = filename)
    filename = "sk.txt"
    open(60,file = filename)
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
    !call evolution_old()
    call cpu_time(t2)

    print *,"lambda is ",lambda
    !  write ck to the file
    g_rmm(1) = 0.0

    do i = 2,n 
       g_rmm(i) = (crmm(i)  +  grmm(i) )/dr(i)  + 1
    end do   !  i
    ! print the results
    do i = 2,n
        write(10,*)dr(i),g_rmm(i) 
        write(20,*)dr(i),grmm(i)/dr(i) 
        write(30,*)dr(i),crmm(i)/dr(i) 
        write(40,*)dr(i),(grmm(i) + crmm(i))/dr(i) 
        write(60,*)dk(i),1.0 + rhom*(gkmm(i) + ckmm(i))/dk(i)
    end do
    print *,"time cost is ",t2 - t1

    close(10)
    close(20)
    close(30)
    close(40)
    close(60)
    !************ check the result *************************
    write(50,*)ch_b,"test",ch_b,ch_b,"grmm"
    do i = 2,n
        write(50,"(2f18.10)")test(i)/dr(i),grmm(i)/dr(i)
    end do
    !******************************************************
    !print *,hkmm(10)
end program main
