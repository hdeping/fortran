module module_evolution
      use module_fst
      implicit none
      contains
!subroutine evolution{{{
subroutine evolution()
    integer                :: jj

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
    filename = "rate.txt"
    open(90,file = filename)
    ckfm = 1.0
    ckff = 1.0
    !***********************************
    do  
        test  = ckfm 
        test1 = ckff
        gkfm(1) = 0.0
        gkff(1) = 0.0
        do   jj = 2,n 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
        !   new closure approximation
        do i = 1,n
            !if(dr(i) < dfm)then
            !else
            crfm(i) = (dr(i) + grfm(i))*mayfm(i)
            crff(i) = (dr(i) + grff(i))*mayff(i)
        end do

        !   fft to calculate ckfm and ckff
        ckfm  = fst(crfm,1)
        ckff  = fst(crff,1)
        ! judge the convergence and exit the cycle

        if(times == 100)then
        do jj = 1,100
            rate = jj*0.01
            lambda  = setion_rate(test ,ckfm,rate)
            lambda1 = setion_rate(test1,ckff,rate)
            write(90,*)rate,lambda,lambda1
        end do
        endif
        lambda  = setion_rate(test ,ckfm,rate)
        lambda1 = setion_rate(test1,ckff,rate)

        !lambda  = conver(test ,ckfm)
        !lambda1 = conver(test1,ckff)

        !lambda  = bisetion(test ,ckfm)
        !lambda1 = bisetion(test1,ckff)

        !lambda2 = conver(test,ckffb)
        lambda = max(lambda,lambda1)
        times = times + 1
        !print *,lambda
        if(mod(times,fre) == 0)then
            print *,times,"lambda = ",lambda
        endif
        !if(times == 3000)exit
        if(lambda < error)exit
        !times = times + 1
        !if (mod(times,1000) == 0)then
        !    print *,"lambda = ",lambda
        !    pause
        !endif
    end do
    close(90)
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
!subroutine evolution_mm{{{
! iteration method to solve OZ equation
subroutine evolution_mm()
    integer                     :: jj
    print *,"get hkmm"
    times = 0
    ! initial ckmm
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

        lambda = setion_rate(test,ckmm,0.99)
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
    !print *,"iteration times is ",times
end subroutine evolution_mm
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
!***********************************************
!  different methods for judging convergence 
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
function setion_rate(a,b,setion)
    real,intent(in)         :: setion
    real(8),intent(in)      :: a(n)
    real(8),intent(inout)   :: b(n)
    real(8)                 :: setion_rate
    integer                 :: ii
    
    do ii = 1,n 
        b(ii) = b(ii)*(1.0 - setion) + a(ii)*setion
    end do
    setion_rate = judge(a,b)
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
    integer                :: jj

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
