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
    filename = "ck.txt"
    open(30,file = filename)
    filename = "hr.txt"
    open(40,file = filename)
    filename = "check.txt"
    open(50,file = filename)
    write(filename,"('sk',i3.3,'.txt')")int(100*rhom)
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
    call cpu_time(t2)

    print *,"lambda is ",lambda
    !  write ck to the file
    g_rmm(1) = 0.0

    do i = 2,n 
       g_rmm(i) = (crmm(i)  +  grmm(i) )/dr(i)  + 1
    end do   !  i
    ! print the results
    skmm(1) = 0.0
    times = 0
    do i  = 2,n
        skmm(i) = rhom*(gkmm(i) + ckmm(i))/dk(i) + 1.0
        if(i == 2)cycle
        if(dk(i) > 60)cycle
        tmp = skmm(i - 1)
        if(tmp > skmm(i - 2).and.tmp > skmm(i))then
            times = times + 1
            print *,times,dk(i)
        endif
    end do

    do i = 2,n
        write(10,*)dr(i),g_rmm(i) 
        write(20,*)dr(i),grmm(i)/dr(i) 
        write(30,*)dr(i),ckmm(i)/dk(i) 
        write(40,*)dr(i),(grmm(i) + crmm(i))/dr(i) 
        write(60,*)dk(i),skmm(i)
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
!test the fst program{{{
    !do i = 1,n
    !    a(i) = (i - 1)*1.0
    !end do
    !write(*,*)"before forward backward fft"
    !b = fst(a,n,l,1)
    !c = fst(b,n,l,-1)
    !do i = 1,n
    !    write(*,"(i3,3f18.9)")i,a(i),b(i),c(i)
    !end do

    !!*******************************************
    !!***test the answer*************************
    !lambda1 = 0
    !lambda2 = 0
    !do i = 1,n
    !    xtmp = hkmm(i) - ckmm(i) - rhom*ckmm(i)*hkmm(i)
    !    ytmp = crmm(i) + dr(i) + grmm(i)
    !    write(50,*)xtmp,ytmp
    !    lambda1 = lambda1 + abs(xtmp)
    !    lambda2 = lambda2 + abs(ytmp)
    !end do
    !lambda1 = lambda1/dble(n)
    !lambda2 = lambda2/dble(n)
    !write(50,*)"lambda1 = ",lambda1
    !write(50,*)"lambda2 = ",lambda2
    !!*******************************************
!}}}
