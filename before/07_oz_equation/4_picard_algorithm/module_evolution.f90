module module_evolution
      use module_fst
      contains
!subroutine evolution_mm{{{
! iteration method to solve OZ equation
subroutine evolution_mm()
    print *,"get hkmm"
    times = 0
    call cpu_time(t1)
    !********************************
    ckmm = 1.0
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

        lambda = conver(test,ckmm)
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
    print *,"lambda = ",lambda
    !print *,"iteration times is ",times
end subroutine evolution_mm
!}}}
!subroutine evolution{{{
subroutine evolution()
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
    rate = 0.9
    !***********************************
    do  
        test  = ckfm 
        test1 = ckff
        gkfm(1) = 0.0
        gkff(1) = 0.0
        do   jj = 2,n 
           chik     = dk(jj) + rhom*hkmm(jj)
           ckffc(jj)= ckff(jj) - ckffb(jj)
           gkfm(jj) = ckfm(jj)*chik/(dk(jj) - rhof*ckffc(jj))&
                      - ckfm(jj)
           gkff(jj) = (dk(jj)**2.0*ckff(jj) + rhom*ckfm(jj)**2.0*chik -&
                      dk(jj)*rhof*ckffc(jj)**2.0)/(dk(jj) - &
                      rhof*ckffc(jj))**2.0 - ckff(jj)
           gkffb(jj)= (dk(jj)**2.0*ckff(jj) + rhom*ckfm(jj)**2.0*chik)/&
                      (dk(jj) - rhof*ckffc(jj))**2.0 - ckffb(jj)
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
        ! judge the convergence and exit the cycle
        !lambda  = conver(test ,ckfm)
        !lambda1 = conver(test1,ckff)
        !if(mod(times,6*fre) == 0 .and. rate > 0.05)then
        !    rate = rate - 0.01
        !endif
        !if(times/fre == 6 )then
        !    rate = 1.01
        !elseif(times/fre == 10)then
        !    rate = .01
        !endif
        lambda  = setion_rate(test ,ckfm,rate)
        lambda1 = setion_rate(test1,ckff,rate)
        !lambda2 = conver(test,ckffb)
        lambda = max(lambda,lambda1)
        times = times + 1
        !print *,lambda
        if(mod(times,fre) == 0)then
            print *,"lambda = ",lambda
        endif
        !if(times == 3000)exit
        if(lambda < error)exit
        !times = times + 1
        !if (mod(times,1000) == 0)then
        !    print *,"lambda = ",lambda
        !    pause
        !endif
    end do
    do jj = 1,n
    hkfm(jj) =  gkfm(jj)  + ckfm(jj) 
    hkff(jj) =  gkff(jj)  + ckff(jj) 
    hkffb(jj)=  gkffb(jj) + ckffb(jj)
    hkffc(jj)=  gkffc(jj) + ckffc(jj)
    end do
    print *,"lambda = ",lambda
end subroutine evolution
!}}}
!subroutine evolution_particle{{{
subroutine evolution_particle()
   print *,"single particle evolution" 
   !!!!!!!!***********************************
    !  get crsfb
    rate  = 0.1
    times = 0
    do i  = 1,n
        if(dr(i) < dsf)then
            crsfb(i) = dr(i) + crmm(i) + grmm(i)
        else
            crsfb(i) = crmm(i)
        endif
    end do
    cksfb = fst(crsfb, 1)
    !  get crssb
    do i = 1,n
        if(dr(i) < dss)then
            crssb(i) = dr(i) + crmm(i) + grmm(i)
        else
            crssb(i) = crmm(i)
        endif
    end do
    ckssb = fst(crssb, 1)
    !***********************************
    cksm  = 1.0
    cksf  = 1.0
    ckss  = 1.0
    do  
        test  = cksm 
        test1 = cksf
        test2 = ckss
        ! get gk with OZ equation
!OZ{{{
        gksm(1) = 0.0
        gksf(1) = 0.0
        gkss(1) = 0.0
        do   jj = 2,n 
            cksfc(jj) = cksf(jj) - cksfb(jj)
            ckssc(jj) = ckss(jj) - ckssb(jj)
            ! get hk
            hksm(jj)  = (dk(jj)*cksm(jj) + rhom*cksm(jj)*hkmm(jj)&
                        + rhof*cksfc(jj)*hkfm(jj))/dk(jj)
            hksfb(jj) = (dk(jj)*cksfb(jj) + rhom*cksm(jj)*hkfm(jj)&
                        + rhof*cksfb(jj)*hkffc(jj) + rhof*&
                        cksfc(jj)*hkffb(jj))/dk(jj)
            hksfc(jj) = (dk(jj)*cksfc(jj) + rhof*cksfc(jj)*hkffc(jj))&
                        /dk(jj)
            hkssb(jj) = (dk(jj)*ckssb(jj) + rhom*cksm(jj)*hksm(jj)&
                        + rhof*cksfc(jj)*hksfb(jj) + rhof*&
                        cksfb(jj)*hksfc(jj))/dk(jj)
            hkssc(jj) = (dk(jj)*ckssc(jj) + rhof*cksfc(jj)*hksfc(jj))&
                        /dk(jj)
            ! get gamma_k 
            gksfb(jj) = hksfb(jj) - cksfb(jj)
            gkssb(jj) = hkssb(jj) - ckssb(jj)
            gksfc(jj) = hksfc(jj) - cksfc(jj)
            gkssc(jj) = hkssc(jj) - ckssc(jj)
            !****************************************
            gksm(jj)  = hksm(jj)  - cksm(jj)
            gksf(jj)  = gksfc(jj) + gksfb(jj)
            gkss(jj)  = gkssc(jj) + gkssb(jj)
        end do
!}}}
        ! inverse fft to calculate gr
        grsm  = fst(gksm, - 1)
        grsf  = fst(gksf, - 1)
        grss  = fst(gkss, - 1)
        !grffb = fst(gkffb,- 1)
        !  calculate cr with PY closure
        do i = 1,n
            crsm(i) = (dr(i) + grsm(i))*maysm(i)
            crsf(i) = (dr(i) + grsf(i))*maysf(i)
            crss(i) = (dr(i) + grss(i))*mayss(i)
        end do

        !   fft to calculate ckfm and ckff
        cksm  = fst(crsm,1)
        cksf  = fst(crsf,1)
        ckss  = fst(crss,1)
        ! judge the convergence and exit the cycle
        !lambda  = conver(test ,cksm)
        !lambda1 = conver(test1,cksf)
        !lambda2 = conver(test2,ckss)

        lambda  = setion_rate(test ,cksm,rate)
        lambda1 = setion_rate(test1,cksf,rate)
        lambda2 = setion_rate(test2,ckss,rate)
        lambda = max(lambda,lambda1,lambda2)
        ! change rate
        !if(mod(times,10*fre) == 0)then
        !    rate = rate - 0.01
        !endif
        times = times + 1
        if(mod(times,fre) == 0)then
            print *,"lambda = ",lambda
        endif
        !if(times == 3000)exit
        if(lambda < error)exit
    end do
    print *,"lambda = ",lambda
end subroutine evolution_particle
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
!function may{{{
function may(d)
    real(8),intent(in)         :: d
    real(8)                    :: may(n)
    integer                    :: ii
    do ii = 1,n
        if(dr(ii) < d)then
             may(ii) = - 1
         else
             may(ii) = 0
         endif
    end do
end function may
!}}}
! judge the convergence
!function judge{{{
! compare the difference between two arrays
! judge the convergence
function judge(a,b)
    real(8),intent(in)      :: a(n)
    real(8),intent(in)      :: b(n)
    real(8)                 :: judge
    integer                 :: ii
    
    judge = 0
    do ii = 1,n 
        judge = judge + abs(a(ii) - b(ii))
    end do
    !judge = judge/dble(n)
    
end function judge
!}}}
!function bisetion{{{
! compare the difference between two arrays
! judge the convergence
function bisetion(a,b)
    real(8),intent(in)      :: a(n)
    real(8),intent(inout)   :: b(n)
    real(8)                 :: bisetion
    integer                 :: ii
    
    bisetion = judge(a,b)
    do ii = 1,n 
        b(ii) = (a(ii) + b(ii))/2.0
    end do
    !judge = judge/dble(n)
    
end function bisetion
!}}}
!function setion_rate{{{
! compare the difference between two arrays
! judge the convergence
function setion_rate(a,b,rate)
    real(8),intent(in)      :: a(n)
    real   ,intent(in)      :: rate
    real(8),intent(inout)   :: b(n)
    real(8)                 :: setion_rate
    integer                 :: ii
    
    setion_rate = judge(a,b)
    do ii = 1,n 
        b(ii) = a(ii)*(1.0 - rate) + b(ii)*rate
    end do
    !judge = judge/dble(n)
    
end function setion_rate
!}}}
!function conver{{{
!  judge convergence with golden setion 
function conver(a,b)
    real(8),intent(in)      :: a(n)
    real(8),intent(inout)   :: b(n)
    real(8)                 :: conver
    integer                 :: ii

    do ii = 1,n
        test1(ii) = a(ii)*gold + b(ii)*(1.0 - gold)
        test2(ii) = a(ii)*(1.0 - gold) + b(ii)*gold 
    end do

    lambda1 = judge(a,test1)
    lambda2 = judge(a,test2)
    if(lambda1 < lambda2)then
        b      = test1
        conver = lambda1
    else
        b      = test2
        conver = lambda2
    endif
    
    
end function conver
!}}}
!subroutine check_ck{{{
! check if the result is the solution to 
! the OZ equation
subroutine check_ck()

    filename = "check.txt"
    open(50,file = filename)

    gkfm(1) = 0.0
    gkff(1) = 0.0
    do   jj = 2,n 
        ckffc(jj) = ckff(jj) - ckffb(jj)
        hkffc(jj) = dk(jj)*ckffc(jj)/(dk(jj) - rhom*ckffc(jj))
        hkfm(jj)  = (dk(jj)*ckfm(jj) + rhom*ckfm(jj)*hkmm(jj))/&
                    (dk(jj) - rhof*ckffc(jj))
        hkffb(jj) = (dk(jj)*ckffb(jj) + rhom*ckfm(jj)*hkfm(jj) + &
                    rhof*ckffb(jj)*hkffc(jj))/(dk(jj) - rhof*ckffc(jj))
        gkffc(jj) = hkffc(jj) - ckffc(jj)
        gkffb(jj) = hkffb(jj) - ckffb(jj)
        gkfm(jj)  = hkfm(jj)  - ckfm(jj)
        gkff(jj)  = gkffc(jj) + gkffb(jj)
    end do
    ! inverse fft to calculate grfm and grff
    grfm  = fst(gkfm, - 1)
    grff  = fst(gkff, - 1)
    grffb = fst(gkffb,- 1)
    !grffb = fst(gkffb,- 1)
    !  calculate crfm and crff with PY closure
    do i = 1,n
        crfm(i) = (dr(i) + grfm(i))*mayfm(i)
        crff(i) = (dr(i) + grff(i))*mayff(i)
    end do

    !   fft to calculate ckfm and ckff

    test  = fst(crfm,1)
    test1 = fst(crff,1)
    do i = 2,n
        write(50,"(3f18.10)")dr(i),test(i)/dr(i),ckfm(i)/dr(i)
    end do

    close(50)
   ! close(60)
end subroutine check_ck
!}}}
end module module_evolution
