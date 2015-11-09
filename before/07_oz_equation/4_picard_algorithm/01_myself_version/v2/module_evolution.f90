module module_evolution
      use module_fst
      contains
!subroutine evolution{{{
! iteration method to solve OZ equation
subroutine evolution()
    ckmm = 1.0
    call cpu_time(t1)
    do jj = 1,30
        rho = 1.06 + jj*0.01
        times = 0
        print *,"rho = ",rho
        write(filename,"('data',i3.3,'.txt')")106 + jj
        open(100,file = filename)
        do 
            test    = ckmm
            test_cr = crmm
            gkmm(1) = 0.0
            do i = 2,n
                gkmm(i) = rho*ckmm(i)**2.0/(dk(i) - rho*ckmm(i))
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
        end do   ! 
        print *,"times = ",times,"lambda = ",lambda
        do i = 2,n
            write(100,*)dr(i),(grmm(i) + crmm(i))/dr(i) 
        end do
        close(100)
    end do  ! jj
end subroutine evolution
!}}}
!subroutine evolution_new{{{
! iteration method to solve OZ equation
subroutine evolution_new()
    ckmm = 1.0
    call cpu_time(t1)

    filename = "g_r_summary.txt"
    open(100,file = filename)
    do jj = 1,100
        rho = 0.06 + jj*0.01
        times = 0
        print *,"rho = ",rho
        do 
            test    = ckmm
            test_cr = crmm
            gkmm(1) = 0.0
            do i = 2,n
                gkmm(i) = rho*ckmm(i)**2.0/(dk(i) - rho*ckmm(i))
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
        end do   ! 
        print *,"times = ",times,"lambda = ",lambda
        i = 102
        tmp = (grmm(i) + crmm(i))/dr(i) + 1.0
            write(100,*)rho,tmp
    end do  ! jj
        close(100)
end subroutine evolution_new
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
