program main
    !use module_fst 
    !use module_common
    use module_evolution
    integer  itmp

    real(8)    :: a(n)
    real(8)    :: b(n)
    real(8)    :: c(n)
    

    !filename = "g_r.txt"
    !open(10,file = filename)
    !filename = "gr.txt"
    !open(20,file = filename)
    !filename = "cr.txt"
    !open(30,file = filename)
    !filename = "hr.txt"
    !open(40,file = filename)
    !filename = "ck.txt"
    !open(70,file = filename)
    !filename = "gk.txt"
    !open(80,file = filename)
    filename = "sk.txt"
    open(80,file = filename)

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
    do i = 2,n
        skmm(i) = rhom*(ckmm(i) + gkmm(i))/dk(i) + 1.0
        write(80,*)dk(i),skmm(i)
    end do

    close(80)

    
    !*********check grmm*************************
    filename = "check_gr.txt"
    open(50,file = filename)
    test = 0.0

    call check_gr()
    write(50,*)"    test    ","    grmm    "
    do i = 2,n
        write(50,"(3f12.6)")test(i)/dr(i),grmm(i)/dr(i)
    end do
    close(50)
    !*********check*************************
    !*********check c(k) *************************
    test = 0.0
    filename = "check_ck.txt"
    open(50,file = filename)

    call check_ck()
    write(50,*)"    test    ","    ckmm    "
    do i = 2,n
        write(50,"(3f12.6)")test(i)/dk(i),ckmm(i)/dk(i)
    end do
    close(50)
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
!}}}
!{{{
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
