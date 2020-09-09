program main
    !use module_fst 
    use module_common
    use module_evolution
    integer  itmp

    !real(8)    :: a(n)
    !real(8)    :: b(n)
    !real(8)    :: c(n)
    !print *,10
    

    filename = "data_g_rmm.txt"
    open(20,file = filename)
    filename = "data_g_rfm.txt"
    open(30,file = filename)
    filename = "data_g_rff.txt"
    open(40,file = filename)
    filename = "hkmm.txt"
    open(50,file = filename)

    !  initial the k and r
    do i = 1,n
        dk(i) = (i - 1)*deltak
        dr(i) = (i - 1)*deltar
    end do

    ! calculate crmm and crffb
! crmm and crffb{{{
    !   initial eta(packing fraction)
    eta = 0.3
    !   calculate crmm
    lambda1 = (1.0 + 2.0*eta)**2/(1 - eta)**4
    lambda1 = (1.0 + 0.5*eta)**2/(1 - eta)**4
    crmm(1) = 0.0
    do i = 2,n
        xtmp = dr(i)/dmm
        if(xtmp < 1.0)then
            crmm(i) = dr(i)*(1.0 - lambda1 - 6*eta*lambda2&
                *xtmp**3)
        else
            crmm(i) = 0.0
        endif
    end do
    !  fft to calculate ckmm
    ckmm = fst(crmm,1)
    !   calculate hkmm with oz equation
    do i = 1,n
        hkmm(i) = ckmm(i)/(dk(i) - rhom*ckmm(i))
    end do
    hrmm     = fst(hkmm,- 1)
    g_rmm(:) = crmm(:) + grmm(:) + dr(:) 
    test     = hkmm
    hkmm     = fst(hrmm,1)
    do  itmp = 1,n
        print "(i5,3f18.9)",itmp,test(itmp),hrmm(itmp),hkmm(itmp)
    end do
    !  calculate crffb 
    crffb(1) = 0.0
    do i = 2,n
        xtmp = dr(i)/dff
        if(xrmp < 1.0)then
            crffb(i) = g_rmm(i)
        else
            crffb(i) = crmm(i)
        endif
    end do
    ckffb = fst(crffb,1)
!}}}

    ! initial the value of c(four arrays)

         ckfm = 1.0
         ckff = 1.0
      ckfm(1) = 0.0
      ckff(1) = 0.0
    j = 0
    times = 0
    !  calculate the mayer function
    mayfm = may(dfm)
    mayff = may(dff)
!********************************************************** 
    call cpu_time(t1)
    !  solve the equation
    call evolution()
    call cpu_time(t2)

    print *,"lambda is ",lambda
    !  write ck to the file
     g_rmm(1) = 0.0 
     g_rfm(1) = 0.0 
     g_rff(1) = 0.0 
    do i = 2,n
       g_rmm(i) = (crmm(i)  +  grmm(i) )/dr(i)  + 1
       g_rfm(i) = (crfm(i)  +  grfm(i) )/dr(i)  + 1
       g_rff(i) = (crff(i)  +  grff(i) )/dr(i)  + 1
    end do   !  i
    do i = 1,n
        write(20,*)dr(i),g_rmm(i) 
        write(30,*)dr(i),g_rfm(i) 
        write(40,*)dr(i),g_rff(i)
    end do
    print *,"times = ",t2 - t1

    close(50)
    close(40)
    close(30)
    close(20)
    !print *,hkmm(10)
end program main
