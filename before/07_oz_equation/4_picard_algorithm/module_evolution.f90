module module_evolution
      use module_fst
      use module_algebra
      use module_judge
      contains
!   newton raphason
!subroutine evolution{{{
subroutine evolution()
    real(8)                :: a(3,3)
    real(8)                :: b(3)
    real(8)                :: deltax(3)
   print *,"other solutions" 
   !!!!!!!!***********************************
    !  get crffb
    do i = 1,n
        if(dr(i) < dff)then
            crffb(i) = dr(i) + crmm(i) + grmm(i)
        else
            crffb(i) = crmm(i)
        endif
    end do
    ckffb = fst(crffb, 1)
    !***********************************
    do  
        test  = ckfm 
        test1 = ckff
        gkfm(1) = 0.0
        gkff(1) = 0.0
        do   jj = 2,n 
           chik     = dk(jj) + rhom*hkmm(jj)
           ckffc(jj)= ckff(jj) - ckffb(jj)
           gkfm(jj) = ckfm(jj)*chik/(dk(jj) - & 
                      rhof*ckffc(jj)) - ckfm(jj)
           gkff(jj) = (dk(jj)**2.0*ckff(jj) + &
                      rhom*ckfm(jj)**2.0*chik -&
                      dk(jj)*rhof*ckffc(jj)**2.0)&
                      /(dk(jj) - rhof*ckffc(jj))**2.0 &
                      - ckff(jj)
           gkffb(jj)= (dk(jj)**2.0*ckff(jj) +& 
                      rhom*ckfm(jj)**2.0*chik)/&
                      (dk(jj) - rhof*ckffc(jj))**2.0 &
                      - ckffb(jj)
           ! newton raphason
           b(1)   = gkfm(jj)
           b(2)   = gkff(jj)
           b(3)   = gkffb(jj)
           a(1,1) = chik/(dk(jj) - rhof*ckffc(jj)) - 1.0
           a(1,2) = ckfm(jj)*(dk(jj) - hkmm(jj)*rhom)*rhof &
                    /(dk(jj) - rhof*ckffc(jj))**2.0
           a(1,3) = - a(1,2)
           a(2,1) = 2*ckfm(jj)*rhom*(dk(jj) + hkmm(jj)*&
                    rhom)*rhof/(dk(jj) - rhof*ckffc(jj))**2.0
           a(2,2) = (dk(jj)**3.0 - ckff(jj)*dk(jj)**2.0*rhof &
                    + 3.0*ckffb(jj)*dk(jj)**2.0*rhof + 2.0*&
                    ckfm(jj)**2.0*rhom*rhof*chik)/&
                    (dk(jj) - rhof*ckffc(jj))**3.0 - 1.0
           a(2,3) = - 2.0*(ckffb(jj)*dk(jj)**2.0+ ckfm(jj)**&
                    2.0*rhom*chik)*rhof/(dk(jj) - rhof*&
                    ckffc(jj))**3.0
           a(3,1) = 2.0*ckfm(jj)*rhom*chik/&
                    (dk(jj) - rhof*ckffc(jj))**2.0
           a(3,2) = - a(2,3)
           a(3,3) = (dk(jj)**3.0 - ckff(jj)*dk(jj)**2.0*rhof &
                    - ckffb(jj)*dk(jj)**2.0*rhof - 2.0*&
                    ckfm(jj)**2.0*rhom*rhof*chik)/&
                    (dk(jj) - rhof*ckffc(jj))**3.0 - 1.0
           deltax = sol_equ(a,b,3)

           gkfm(jj)   = gkfm(jj)  - deltax(1)
           gkff(jj)   = gkff(jj)  - deltax(2) 
           gkffb(jj)  = gkffb(jj) - deltax(3) 
        end do

        ! inverse fft to calculate grfm and grff
        grfm  = fst(gkfm, - 1)
        grff  = fst(gkff, - 1)
        !grffb = fst(gkffb,- 1)
        !  calculate crfm and crff with PY closure
        do i = 1,n
            crfm(i) = (dr(i) + grfm(i))*mayfm(i)
            crff(i) = (dr(i) + grff(i))*mayff(i)
        end do

        !   fft to calculate ckfm and ckff
        ckfm  = fst(crfm,1)
        ckff  = fst(crff,1)
        ! judge convergence
        
        lambda1 = conver(test,ckfm)
        lambda2 = conver(test1,ckff)
        if(mod(times,fre) == 0)then
            print *,"lambda = ",lambda
        endif
           if(lambda < error)exit
        if(lambda < error)exit
        !times = times + 1
        !if (mod(times,1000) == 0)then
        !    print *,"lambda = ",lambda
        !    pause
        !endif
    end do
end subroutine evolution
!}}}
!subroutine evolution_gr{{{
! iteration method to solve OZ equation
subroutine evolution_gr()
    times = 0
    ckmm = 1.0
    call cpu_time(t1)
    do 
        test    = grmm
        test_cr = crmm
        ! calculate crmm with PY closure
        do i = 1,n
            crmm(i) = (dr(i) + grmm(i))*maymm(i)
        end do
        ! fft to calculate ckmm
        ckmm = fst(crmm,1)
        gkmm(1) = 0.0
        do i = 2,n
           gkmm(i) = rhom*ckmm(i)**2.0/(dk(i) - rhom*ckmm(i))
        end do
        ! inverse fft to calculate grmm
        grmm = fst(gkmm,- 1)
        ! judge the convergence
        ! get a mean value of ckmm  (golden setion)
         times = times + 1

        lambda = conver(test,grmm)
        !if(times > int(1E4))exit
        if(mod(times,1000) == 0)then
            print *,"lambda = ",lambda
        endif
        if(lambda < error)exit
        !print *,"rhom = ",rhom
        !pause
    end do
end subroutine evolution_gr
!}}}
!!subroutine evolution_mm_approx{{{
!! iteration method to solve OZ equation
!subroutine evolution_mm()
!    print *,"get hkmm"
!    eta     = 0.5
!    lambda1 =  (1.0 + 2*eta)**2.0/(1.0 - eta)**4.0
!    lambda2 =  - (1.0 + 0.5*eta)**2.0/(1.0 - eta)**4.0
!    !print *,lambda1,lambda2
!    !pause
!    do i = 1,n
!        xtmp    = dr(i)/dmm
!        if(xtmp < 1.0)then
!        crmm(i) = dr(i)*(- lambda1 - 6.0*eta*lambda2*dr(i)-&
!                  0.5*eta*lambda1*dr(i)**3.0)
!        else
!            crmm(i) = 0.0
!        endif
!        !print *,i,crmm(i)
!    end do
!    
!    ckmm = fst(crmm,1)
!    !  get hkmm with OZ equation
!    hkmm(1) = 0
!    do i = 2,n
!        hkmm(i) = dk(i)*ckmm(i)/(dk(i) - rhom*ckmm(i))  
!    end do
!    ! get hrmm
!    hrmm = fst(hkmm,- 1)
!    !do i = 1,n
!    !    hkmm(i) = ckmm(i) + gkmm(i)
!    !end do
!    !print *,"iteration times is ",times
!end subroutine evolution_mm
!!}}}
!subroutine evolution_mm{{{
! iteration method to solve OZ equation
subroutine evolution_mm()
    print *,"get hkmm"
    times = 0
    ckmm = 1.0
    call cpu_time(t1)
    do 
        test    = ckmm
        test_cr = crmm
        gkmm(1) = 0.0
        do i = 2,n
            gkmm(i) = rhom*ckmm(i)**2.0/(dk(i) - rhom*ckmm(i))
        end do
        ! inverse fft to calculate grmm
        grmm = fst(gkmm,- 1)
        ! calculate crmm with PY closure
        do i = 1,n
            crmm(i) = (dr(i) + grmm(i))*maymm(i)
        end do
        ! fft to calculate ckmm
        ckmm = fst(crmm,1)
        ! judge the convergence
        ! get a mean value of ckmm  (golden setion)
         times = times + 1

        lambda = judge(test,ckmm)
        !if(times > int(1E4))exit
        if(mod(times,1000) == 0)then
            print *,"lambda = ",lambda
        endif
        if(lambda < error)exit
        !print *,"rhom = ",rhom
        !pause
    end do
    do i = 1,n
        hkmm(i) = ckmm(i) + gkmm(i)
    end do
    hrmm = fst(hkmm,- 1)
    !print *,"iteration times is ",times
end subroutine evolution_mm
!}}}
!subroutine evolution_mm_NR{{{
! iteration method to solve OZ equation
subroutine evolution_mm_NR()
    real(8)             :: deltax
    print *,"get hkmm"
    times = 0
    ckmm = 1.0
    gkmm = 1.0
    call cpu_time(t1)
    do 
        test    = ckmm
        test_cr = crmm
        gkmm(1) = 0.0
        do i = 2,n
            ctmp = dk(i) - rhom*ckmm(i)
            xtmp = rhom*ckmm(i)**2.0/ctmp
            xtmp = gkmm(i) - xtmp
            ytmp = rhom*ckmm(i)*(dk(i) + ctmp)/&
                   ctmp**2.0
            deltax  = xtmp/ytmp
            gkmm(i) = gkmm(i) - deltax
        end do
        ! inverse fft to calculate grmm
        grmm = fst(gkmm,- 1)
        ! calculate crmm with PY closure
        do i = 1,n
            crmm(i) = (dr(i) + grmm(i))*maymm(i)
        end do
        ! fft to calculate ckmm
        ckmm = fst(crmm,1)
        ! judge the convergence
        ! get a mean value of ckmm  (golden setion)
         times = times + 1

        lambda = judge(test,ckmm)
        !if(times > int(1E4))exit
        if(mod(times,10) == 0)then
            print *,"lambda = ",lambda
        endif
        if(lambda < error)exit
        !print *,"rhom = ",rhom
        !pause
    end do
    do i = 1,n
        hkmm(i) = ckmm(i) + gkmm(i)
    end do
    hrmm = fst(hkmm,- 1)
    !print *,"iteration times is ",times
end subroutine evolution_mm_NR
!}}}
!subroutine check_gr {{{
! check if the result is the solution to 
! the OZ equation
subroutine check_gr()

    filename = "check.txt"
    open(50,file = filename)
!{{{Read data
   ! filename = "gr.txt"
   ! open(60,file = filename,status = "old",iostat = ierror)

    ! read c(r) from file
   ! grmm(1) = 0.0
   ! do i = 2,n
   !     read(60,*,iostat = ierror)tmp,ctmp
   !     grmm(i) = ctmp*dr(i) 
   ! end do
!}}}
    ! calculate crmm with PY closure 
    do i = 1,n
       crmm(i) = (dr(i) + grmm(i))*maymm(i)
    end do
    ! calculate ckmm with fft
    ckmm = fst(crmm, 1)
    ! calculate gkmm with OZ equation
    gkmm(1) = 0.0
    do i = 2,n
        gkmm(i) = rhom*ckmm(i)**2.0/(dk(i) - rhom*ckmm(i))   
    end do
    ! calculate a new grmm(represented by test)
    test = fst(gkmm,- 1)
    write(50,*)"test           ","            grmm"
    do i = 2,n
        write(50,"(3f18.10)")dr(i),test(i)/dr(i),grmm(i)/dr(i)
    end do

    close(50)
   ! close(60)
end subroutine check_gr
!}}}
!subroutine check_ck {{{
! check if the result is the solution to 
! the OZ equation
subroutine check_cr()

    filename = "check.txt"
    open(50,file = filename)
!Read data{{{
   ! filename = "gr.txt"
   ! open(60,file = filename,status = "old",iostat = ierror)

    ! read c(r) from file
   ! grmm(1) = 0.0
   ! do i = 2,n
   !     read(60,*,iostat = ierror)tmp,ctmp
   !     grmm(i) = ctmp*dr(i) 
   ! end do
!}}}
    ! calculate gkmm with OZ equation
    gkmm(1) = 0.0
    do i = 2,n
        gkmm(i) = rhom*ckmm(i)**2.0/(dk(i) - rhom*ckmm(i))   
    end do
    ! calculate a new grmm(represented by test)
    grmm = fst(gkmm, - 1)
    ! calculate crmm with PY closure 
    do i = 1,n
       crmm(i) = (dr(i) + grmm(i))*maymm(i)
    end do
    ! calculate ckmm with fft
    test = fst(crmm, 1)

    write(50,*)"test           ","grmm"
    do i = 2,n
        write(50,"(3f18.10)")dr(i),test(i)/dr(i),ckmm(i)/dr(i)
    end do


    close(50)
   ! close(60)
end subroutine check_cr
!}}}
end module module_evolution
!before {{{
! iteration to get hkmm
    !times = 0
    !ckmm = 1.0
    !call cpu_time(t1)
    !do 
    !    test    = ckmm
    !    test_cr = crmm
    !    gkmm(1) = 0.0
    !    do i = 2,n
    !        gkmm(i) = rhom*ckmm(i)**2.0/(dk(i) - rhom*ckmm(i))
    !    end do
    !    ! inverse fft to calculate grmm
    !    grmm = fst(gkmm,- 1)
    !    ! calculate crmm with PY closure
    !    do i = 1,n
    !        crmm(i) = (dr(i) + grmm(i))*maymm(i)
    !    end do
    !    ! fft to calculate ckmm
    !    ckmm = fst(crmm,1)
    !    ! judge the convergence
    !    ! get a mean value of ckmm  (golden setion)
    !     times = times + 1

    !    lambda = judge(test,ckmm)
    !    !if(times > int(1E4))exit
    !    if(mod(times,1000) == 0)then
    !        print *,"lambda = ",lambda
    !    endif
    !    if(lambda < error)exit
    !    !print *,"rhom = ",rhom
    !    !pause
    !end do
!}}}
