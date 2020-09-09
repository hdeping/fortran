module module_evolution
      use module_fst
      contains
!subroutine evolution{{{
subroutine evolution()
   print *,"other solutions" 
    do  
        test = ckfm 
        gkfm(1) = 0.0
        gkff(1) = 0.0
        do   jj = 2,n 
            chik      = dk(jj) + rhom*hkmm(jj)
            gkfm(jj)  = ckfm(jj)*chik/(dk(jj) - rhof*&
                        ckff(jj) + rhof*ckffb(jj)) -ckfm(jj)
            gkff(jj)  = (dk(jj)**2*ckff(jj) + rhom*&
                        ckfm(jj)**2 - dk(jj)*rhof *(ckff(jj) &
                        - ckffb(jj))**2)/(dk(jj) - rhof*ckff(jj)&
                        + rhof*ckffb(jj)) - ckff(jj)
            gkffb(jj) = (dk(jj)**2*ckff(jj) + rhom*ckfm(jj)**2)/&
                        (dk(jj) - rhof*ckff(jj) +rhof*ckffb(jj)) - ckffb(jj)
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
            !***************************** Mention *********
            crffb(i) = (dr(i) + grmm(i))*mayff(i)
            !***************************** Mention *********
        end do
        !   fft to calculate ckfm and ckff
        ckfm  = fst(crfm,1)
        ckff  = fst(crff,1)
        gkffb = fst(grffb,- 1)
        ! judge the convergence and exit the cycle
        lambda = judge(test,ckfm)
        times = times + 1
        print *,lambda
        if(mod(times,fre) == 0)then
            print *,"lambda = ",lambda
        endif
        if(times == 1400)exit
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
