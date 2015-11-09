module module_evolution
      use module_fst
      use module_algebra
      contains
!subroutine evolution{{{
subroutine evolution()
    integer              :: ii
    integer              :: inum
    real(8)              :: pre_lambda 
    wp = getwp()
    ! give an initial value of ck
    times = 0
    ! get a suitable rate
    ckpp  = 1.0
    ckwp  = 1.0
    ckww  = 1.0
    !rate  = getrate(0D0,1D0)
    rate  = 0.9
    inum  = 0
    do 
        ! get the solution of OZ 
        lambda = getoz()
        if(lambda < error)exit
        times = times + 1
        if(mod(times,fre) == 0)then
            print *,"lambda = ",lambda
        endif
        !if(mod(times,fre) == 0)then
        !    print *,rate, "lambda = ",lambda
        !    if(times/fre == 1 .or. inum == 1)then
        !        pre_lambda = lambda 
        !        inum = 0
        !    else
        !        if(pre_lambda < lambda)then
        !            rate = getrate(rate - 1D-1,rate + 1D-1)
        !            inum = 1
        !        endif
        !    endif
        !endif
    end do


    do ii = 1,n
        hkpp(ii) = gamma_kpp(ii) +  ckpp(ii)
        hkwp(ii) = gamma_kwp(ii) +  ckwp(ii)
        hkww(ii) = gamma_kww(ii) +  ckww(ii)
    end do
end subroutine evolution
!}}}
!subroutine evolution_nano{{{
subroutine evolution_nano()
    integer              :: ii
    print *,"nano particle part"
    ! give an initial value of ck
    cknp  = 1.0
    cknw  = 1.0
    times = 0
    do 
        test1  = cknp 
        test2  = cknw 
        ! get gamma_k with OZ equation
        gamma_knp(1) = 0.0
        gamma_knw(1) = 0.0
        do ii = 2,n
            hknp(ii) = (dk(ii)*wp(ii)*cknp(ii) + rhop*cknp(ii)*&
                       hkpp(ii) + rhow*hkwp(ii)*cknw(ii))/dk(ii)
            hknw(ii) = (dk(ii)*cknw(ii) + rhop*cknp(ii)*hkwp(ii)&
                       + rhow*hkww(ii)*cknw(ii))/dk(ii)
            ! calculate gamma_k
            gamma_knp(ii) = hknp(ii) - cknp(ii)
            gamma_knw(ii) = hknw(ii) - cknw(ii)
        end do
        ! calculate gamma_r with inverse fft
        gamma_rnp = fst(gamma_knp,- 1)
        gamma_rnw = fst(gamma_knw,- 1)
        ! calculate cr with PY
        do ii = 1,n
            crnp(ii) = (dr(ii) + gamma_rnp(ii))*maynp(ii)
            crnw(ii) = (dr(ii) + gamma_rnw(ii))*maynw(ii)
        end do
        ! calculate ck with fft
        cknp = fst(crnp,1)
        cknw = fst(crnw,1)
        ! judge the convergence
        !lambda  = conver(test,ckpp)
        !lambda1 = conver(test1,ckwp)
        !lambda2 = conver(test2,ckww)

        lambda1 = setion_rate(test1,cknp,0.9)
        lambda2 = setion_rate(test2,cknw,0.9)
        !lambda  = judge(test,ckpp)
        !lambda1 = judge(test1,ckwp)
        !lambda2 = judge(test2,ckww)
        lambda_final = max(lambda1,lambda2)
        !print *,times,"lambda = ",lambda_final
        if(lambda_final < error)exit
        times = times + 1
        if(mod(times,fre/10) == 0)then
            print *,"lambda = ",lambda_final
        endif
    end do
end subroutine evolution_nano
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
!function getwp{{{
function getwp()
    real(8)             :: getwp(n)
    real(8)             :: f
    integer             :: ii
    
    do ii = 1,n
        f         = sin(dk(ii)*b0)/(dk(ii)*b0)
        getwp(ii) = (nnum*(1 - f**2.0) - 2*f + 2*&
                    f**(dble(nnum) + 1.0))/(nnum*(1 - f)**2.0)
    end do


end function getwp
!}}}
!subroutine check{{{
subroutine check()
    integer              :: ii
    filename = "check.txt"
    open(10,file = filename)
        ! get gamma_k with OZ equation
!OZ{{{
        gamma_kpp(1) = 0.0
        gamma_kwp(1) = 0.0
        gamma_kww(1) = 0.0
        do ii = 2,n
            ! get the value of array a and b 
            a(1,1) = dk(ii) - rhop*wp(ii)*ckpp(ii) 
            a(1,2) = - rhow*wp(ii)*ckwp(ii)
            a(1,3) = 0.0
            a(2,1) = - rhop*ckwp(ii)
            a(2,2) = dk(ii) - rhow*wp(ii)*ckww(ii) 
            a(2,3) = 0.0
            a(3,1) = 0.0
            a(3,2) = - rhop*ckwp(ii)
            a(3,3) = dk(ii) - rhow*wp(ii)*ckww(ii) 

            b(1)   = dk(ii)*wp(ii)**2.0*ckpp(ii)
            b(2)   = dk(ii)*wp(ii)*ckwp(ii)
            b(3)   = dk(ii)*ckww(ii)
             
            x        = sol_equ(a,b,ntmp)
            !hkpp(ii) = x(1)
            !hkwp(ii) = x(2)
            !hkww(ii) = x(3)

            ! calculate gamma_k
            gamma_kpp(ii) = x(1) - ckpp(ii)
            gamma_kwp(ii) = x(2) - ckwp(ii)
            gamma_kww(ii) = x(3) - ckww(ii)
        end do
!}}}
        ! calculate gamma_r with inverse fft
        gamma_rpp = fst(gamma_kpp,- 1)
        gamma_rwp = fst(gamma_kwp,- 1)
        gamma_rww = fst(gamma_kww,- 1)
        ! calculate cr with PY
        do ii = 1,n
            crpp(ii) = (dr(ii) + gamma_rpp(ii))*maypp(ii)
            crwp(ii) = (dr(ii) + gamma_rwp(ii))*maywp(ii)
            crww(ii) = (dr(ii) + gamma_rww(ii))*mayww(ii)
        end do
        ! calculate ck with fft
        test  = fst(crpp,1)
        test1 = fst(crwp,1)
        test2 = fst(crww,1)
        do ii = 1,n
            write(10,"(6f10.5)")test(ii),ckpp(ii),test1(ii),ckwp(ii),test2(ii),ckww(ii)
        end do
        close(10)
end subroutine check 
!}}}
!subroutine getoz{{{
function getoz()
    real(8)              :: getoz
    integer              :: ii
    test   = ckpp 
    test1  = ckwp 
    test2  = ckww 
    ! get gamma_k with OZ equation
    gamma_kpp(1) = 0.0
    gamma_kwp(1) = 0.0
    gamma_kww(1) = 0.0
    do ii = 2,n
        ! get the value of array a and b 
        a(1,1) = dk(ii) - rhop*wp(ii)*ckpp(ii) 
        a(1,2) = - rhow*wp(ii)*ckwp(ii)
        a(1,3) = 0.0
        a(2,1) = - rhop*ckwp(ii)
        a(2,2) = dk(ii) - rhow*wp(ii)*ckww(ii) 
        a(2,3) = 0.0
        a(3,1) = 0.0
        a(3,2) = - rhop*ckwp(ii)
        a(3,3) = dk(ii) - rhow*wp(ii)*ckww(ii) 

        b(1)   = dk(ii)*wp(ii)**2.0*ckpp(ii)
        b(2)   = dk(ii)*wp(ii)*ckwp(ii)
        b(3)   = dk(ii)*ckww(ii)
         
        x      = sol_equ(a,b,ntmp)
        ! calculate gamma_k
        gamma_kpp(ii) = x(1) - ckpp(ii)
        gamma_kwp(ii) = x(2) - ckwp(ii)
        gamma_kww(ii) = x(3) - ckww(ii)
    end do
    ! calculate gamma_r with inverse fft
    gamma_rpp = fst(gamma_kpp,- 1)
    gamma_rwp = fst(gamma_kwp,- 1)
    gamma_rww = fst(gamma_kww,- 1)
    ! calculate cr with PY
    do ii = 1,n
        crpp(ii) = (dr(ii) + gamma_rpp(ii))*maypp(ii)
        crwp(ii) = (dr(ii) + gamma_rwp(ii))*maywp(ii)
        crww(ii) = (dr(ii) + gamma_rww(ii))*mayww(ii)
    end do
    ! calculate ck with fft
    ckpp = fst(crpp,1)
    ckwp = fst(crwp,1)
    ckww = fst(crww,1)
    lambda  = setion_rate(test, ckpp,rate)
    lambda1 = setion_rate(test1,ckwp,rate)
    lambda2 = setion_rate(test2,ckww,rate)
    getoz   = max(lambda,lambda1,lambda2)
end function getoz
!}}}
!function getrate{{{
function getrate(a,b)
    integer,parameter   :: total = 20
    real(8),intent(in)  :: a
    real(8),intent(in)  :: b
    real(8)             :: getrate
    real(8)             :: lambda_tmp ! for the judge
    integer             :: ii
    integer             :: jj
    real(8)             :: tmpww(n)
    real(8)             :: tmpwp(n)
    real(8)             :: tmppp(n)
    real(8)             :: tmpww_1(n)
    real(8)             :: tmpwp_1(n)
    real(8)             :: tmppp_1(n)

    
    filename =  "rate.txt"
    open(100,file = filename,status = "old")
    tmppp_1  =  ckpp
    tmpwp_1  =  ckwp
    tmpww_1  =  ckww
    call random_seed()
    do ii = 1,total
        call random_number(rate)
        rate = rate*(b - a) + a
        ckpp = tmppp_1
        ckwp = tmpwp_1
        ckww = tmpww_1
        do jj = 1,fre
            lambda = getoz()
        end do
        if(ii == 1)then
           lambda_tmp = lambda 
           getrate    = rate
           tmpww      = ckww
           tmpwp      = ckwp
           tmppp      = ckpp
        else
            if(lambda < lambda_tmp)then
               lambda_tmp = lambda
               getrate    = rate
               tmpww      = ckww
               tmpwp      = ckwp
               tmppp      = ckpp
            endif
        endif
        write(100,"('rate = ',f12.6,'lambda = ',f12.6)")rate,lambda
    end do
    ckww = tmpww
    ckwp = tmpwp
    ckpp = tmppp
    close(100)

end function getrate
!}}}
!****** judge the convergence **************
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
    real(8),intent(inout)   :: b(n)
    real,intent(in)         :: rate
    real(8)                 :: setion_rate
    integer                 :: ii
    
    setion_rate = judge(a,b)
    do ii = 1,n 
        b(ii) = a(ii)*(1.0 - rate) + b(ii)*rate
    end do
    
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
end module module_evolution
