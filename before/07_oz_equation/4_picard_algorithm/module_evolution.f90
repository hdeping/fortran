module module_evolution
      use module_fst
      contains
!subroutine evolution{{{
! iteration method to solve OZ equation
subroutine evolution()
    times = 0
    ckmm  = 1.0
    rate  = 0.1
    call cpu_time(t1)
    call solve_OZ(ckmm,lambda,ckmm)
    !  begin the evolution
    do 
        in_ck1(:)  = rate*test(:) + (1.0 - rate)*ckmm(:)
        in_ck2(:)  = rate*ckmm(:) + (1.0 - rate)*test(:)
        call solve_OZ(out_ck1,lambda1,in_ck1)
        call solve_OZ(out_ck2,lambda2,in_ck2)
        if(lambda1 < lambda2)then
            ckmm   = out_ck1
            test   = in_ck1
            lambda = lambda1
        else
            ckmm   = out_ck2
            test   = in_ck2
            lambda = lambda2
        endif
        if(lambda < error)exit
        times = times + 1
        if(mod(times,fre) == 0 )then
            print *,"lambda = ",lambda
        endif
    end do
    call cpu_time(t2)
    print "('time cost is ',1f10.3,' seconds')",t2 - t1
end subroutine evolution
!}}}
!subroutine solve_OZ{{{
! iteration method to solve OZ equation
subroutine solve_OZ(out_ck,lambda,in_ck)
    real(8),intent(in)            :: in_ck(n)
    real(8),intent(out)           :: out_ck(n)
    real(8),intent(out)           :: lambda
    test    = in_ck
    gkmm(1) = 0.0
    do i = 2,n
        gkmm(i) = rhom*in_ck(i)**2.0/(dk(i) - rhom*in_ck(i))
    end do
    ! inverse fft to calculate grmm
    grmm = fst(gkmm,- 1)
    ! calculate crmm with PY closure
    do i = 1,n
        crmm(i) = (dr(i) + grmm(i))*maymm(i)
    end do
    ! fft to calculate in_ck
    out_ck = fst(crmm,1)
    lambda = judge(test,out_ck)
end subroutine solve_OZ
!}}}
!subroutine evolution_old{{{
! iteration method to solve OZ equation
subroutine evolution_old()
    times = 0
    ckmm  = 1.0
    rate  = 0.5
    call cpu_time(t1)
    call solve_OZ(ckmm,lambda,ckmm)
    !  begin the evolution_old
    do 
        call solve_OZ(ckmm,lambda,ckmm)
        ckmm = ckmm*(1.0 - rate) + test*rate
        if(lambda < error)exit
        times = times + 1
        if(mod(times,fre) == 0 )then
            print *,"lambda = ",lambda
        endif
    end do
    call cpu_time(t2)
    print "('time cost is ',1f10.3,' seconds')",t2 - t1
end subroutine evolution_old
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
!************ judge the convergence **********
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
    print *,lambda1, lambda2,lambda1/lambda2
    pause
    if(lambda1 < lambda2)then
        b      = test1
        conver = lambda1
    else
        b      = test2
        conver = lambda2
    endif
    
    
end function conver
!}}}
!************ check the results **********
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
!subroutine new_check {{{
! check if the result is the solution to 
! the OZ equation
!  test hrmm and give a new one
!  compare the two 
subroutine new_check()
    include     "omp_lib.h"

    integer                     :: ii
    integer                     :: jj
    integer                     :: kk
    real(8)                     :: tmp1
    real(8)                     :: tmp2

    !$omp parallel do
    do ii = 2,n
        tmp1 = 0.0
        !$omp parallel do
        do jj = 2,n
            tmp2 = 0.0
            !$omp parallel do
            do kk = abs(ii - jj),ii + jj
                if(kk == 0)cycle
                if(kk > n)exit
                tmp2 = tmp2 + deltar*dr(kk)*crmm(kk)
            end do
            !$omp end parallel do
            tmp1 = tmp1 + deltar*dr(jj)*hrmm(jj)*tmp2
        end do
        !$omp end parallel do
        test(ii) = crmm(ii) + 2*pi*rhom*tmp1/dr(ii)
    end do
    !$omp end parallel do
    lambda  = judge(test,hrmm)
    lambda1 = 0.0
    do ii = 102,n
        lambda1 = lambda1 + abs(hrmm(ii) - test(ii))
    end do
end subroutine new_check
!}}}

end module module_evolution
