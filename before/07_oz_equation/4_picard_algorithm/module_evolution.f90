module module_evolution
      use module_fst
      contains
!subroutine evolution{{{
! iteration method to solve OZ equation
subroutine evolution()
    real(8)               :: rate
    times = 0
    ckmm = 1.0
    call cpu_time(t1)
    rate  = 0.6
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

        !lambda = conver(test,ckmm)
        lambda  = ratesection(test,ckmm,rate)
        !if(times > int(1E4))exit
        if(mod(times,1000) == 0)then
            print *,"lambda = ",lambda
        endif
        if(lambda < error)exit
        !print *,"rhom = ",rhom
        !pause
    end do
    print *,"iteration times is ",times
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
!function may{{{
function may(d)
    real(8),intent(in)         :: d
    real(8)                    :: may(n)
    real(8)                    :: potential
    real(8)                    :: x
    integer                    :: ii
    do ii = 1,n
        if(ii == 1)then
            may(ii) = - 1D0
        else
            x         = d/dr(ii)
            potential = 4.0*lj*(x**12.0 - x**6.0) + lj
            may(ii)   = exp(- beta*potential) - 1D0
        endif
    end do
end function may
!}}}
!function may_hard{{{
function may_hard(d)
    real(8),intent(in)         :: d
    real(8)                    :: may_hard(n)
    real(8)                    :: potential
    real(8)                    :: x
    integer                    :: ii
    do ii = 1,n
        if(dr(ii) <= d) then
            may_hard(ii) = - 1.0
        else
            may_hard(ii) = 0.0
        endif
    end do
end function may_hard
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
!function ratesection{{{
! compare the difference between two arrays
! judge the convergence
function ratesection(a,b,rate)
    real(8),intent(in)      :: a(n)
    real(8),intent(in)      :: rate
    real(8),intent(inout)   :: b(n)
    real(8)                 :: ratesection
    integer                 :: ii
    
    ratesection = 0.0
    do ii = 1,n
        ratesection = ratesection + &
            abs(a(ii) - b(ii))
    end do
    do ii = 1,n 
        b(ii) = a(ii)*rate + b(ii)*(1.0 - rate)
    end do
    !judge = judge/dble(n)
    
end function ratesection
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
!subroutine check {{{
! check if the result is the solution to 
! the OZ equation
subroutine check()

    filename = "gr.txt"
    open(60,file = filename,status = "old",iostat = ierror)

    ! read c(r) from file
    grmm(1) = 0.0
    do i = 2,n
        read(60,*,iostat = ierror)tmp,ctmp
        grmm(i) = ctmp*dr(i) 
    end do
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

    close(60)
end subroutine check
!}}}
end module module_evolution
