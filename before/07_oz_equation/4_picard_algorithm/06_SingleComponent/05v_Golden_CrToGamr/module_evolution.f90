module module_evolution
    use  module_fst
    contains
!subroutine evolution{{{
subroutine evolution()
    times = 0
    grmm = 1.0
    call cpu_time(t1)
    do 
        test = grmm
        ! calculate crmm with PY closure
        do i = 1,n
            crmm(i) = (dr(i) + grmm(i))*maymm(i)
        end do
        ! fft to calculate ckmm
        ckmm = fst(crmm,1)
        do i = 2,n
            gkmm(i) = rhom*ckmm(i)**2.0/(dk(i) - rhom*ckmm(i))
        end do
        ! inverse fft to calculate grmm
        grmm = fst(gkmm,- 1)
        ! judge the convergence
        ! get a mean value of grmm
                
        
        times = times + 1
        lambda = conver(test,grmm)
        !if(times > int(1E4))exit
        if(mod(times,1000) == 0)then
            print *,"lambda = ",lambda
        endif
        if(lambda < error)exit
    end do
end subroutine evolution
!}}}
!subroutine check gr{{{
subroutine check_gr()

    filename = "gr.txt"
    open(60,file = filename,status = "old",iostat = ierror)

    ! read c(r) from file
    do i = 1,n
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
end subroutine check_gr
!}}}
!subroutine check ck{{{
subroutine check_ck()

    filename = "ck.txt"
    open(60,file = filename,status = "old",iostat = ierror)

    ! read ck from file
    ckmm(1) = 0.0
    do i = 2,n
        read(60,*,iostat = ierror)tmp,ctmp
        ckmm(i) = ctmp*dk(i) 
    end do
    ! calculate gkmm with OZ equation
    gkmm(1) = 0.0
    do i = 2,n
        gkmm(i) = rhom*ckmm(i)**2.0/(dk(i) - rhom*ckmm(i))   
    end do
    ! calculate grmm(i) with fft
    grmm = fst(gkmm,- 1)
    ! calculate crmm with PY closure 
    do i = 1,n
       crmm(i) = (dr(i) + grmm(i))*maymm(i)
    end do
    ! calculate a new ckmm(represented by test)
    test = fst(crmm, 1)

    close(60)
end subroutine check_ck
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
function conver(a,b)
    real(8),intent(in)      :: a(n)
    real(8),intent(inout)   :: b(n)
    real(8)                 :: conver
    integer                 :: ii
    
    conver = 0
    do ii = 1,n 
        conver = conver + abs(a(ii) - b(ii))
    end do
    do i = 1,n
        b(i) = a(i)*gold + b(i)*(1.0 - gold)
    end do

    
end function conver
!}}}
end module module_evolution
