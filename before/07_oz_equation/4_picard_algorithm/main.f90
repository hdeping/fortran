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
    filename = "others.txt"
    open(20,file = filename)
    filename = "sk.txt"
    open(30,file = filename)
    filename = "check.txt"
    open(50,file = filename)
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
    !call evolution()
    call evolution_old()
    call cpu_time(t2)

    !  write ck to the file
    g_rmm(1) = 0.0

    do i = 2,n 
        crmm(i) = crmm(i)/dr(i)
        grmm(i) = grmm(i)/dr(i)
        hrmm(i) = crmm(i) + grmm(i)
        hkmm(i) = (ckmm(i) + gkmm(i))/dk(i)
       g_rmm(i) = hrmm(i) + 1
    end do   !  i
    ! print the results
    do i = 2,n
        write(10,*)dr(i),g_rmm(i) 
        write(20,"(4f12.5)")dr(i),grmm(i) ,crmm(i) ,hrmm(i)
        ! get sk
        write(30,"(4f12.5)")dk(i),1.0 + rhom*hkmm(i)
    end do
    print *,"time cost is ",t2 - t1

    close(10)
    close(20)
    
    !************ check the result *************************
    call new_check()
    do i = 2,n
        write(50,"(3f12.5)")dr(i),hrmm(i),test(i)
    end do
    close(50)
    print *,"test lambda = ",lambda
    print *,"test lambda1 = ",lambda1
    !******************************************************
end program main
