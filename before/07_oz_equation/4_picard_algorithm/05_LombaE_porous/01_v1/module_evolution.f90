module module_evolution
    use module_fst
    use module_common
    integer       ::  ii
    integer       ::  jj
    contains
!function may{{{
function may(d)
    real(8)       ::  d
    real(8)       ::  may(n)

    do ii = 1,n
        xtmp = dr(ii)/d
        if(xtmp < 1)then
            may(ii) = - 1
        else
            may(ii) = 0
        endif
    end do
    
end function may
!}}}
!subroutine evolution{{{
subroutine evolution()
    
    do  
        test = ckfm 
        gkfm(1) = 0.0
        gkff(1) = 0.0
        do  jj = 2,n 
            chik = dk(jj) + rhom*hkmm(jj)
            gkfm(jj) = ckfm(jj)*chik/(dk(jj) - rhof*&
              ckff(jj) + rhof*ckffb(jj)) -ckfm(jj)
            gkff(jj) = (dk(jj)**2*ckff(jj) + rhom*&
              ckfm(jj)**2 - dk(jj)*rhof *(ckff(jj) &
              - ckffb(jj))**2)/(dk(jj) - rhof*ckff(jj)&
              + rhof*ckffb(jj)) - ckff(jj)
         !   gkffb(jj) = (dk(jj)**2*ckff(jj) + rhom*ckfm(jj)**2)/&
         !      (dk(jj) - rhof*ckff(jj) +rhof*ckffb(jj)) - ckffb(jj)
        end do
        ! inverse fft to calculate grfm and grff
        grfm = fst(gkfm,- 1)
        grff = fst(gkff,- 1)
        !grffb = fst(gkffb,- 1)
        !  calculate crfm and crff with PY closure
        do i = 1,n
            crfm(i) = (dr(i) + grfm(i))*mayfm(i)
            crff(i) = (dr(i) + grff(i))*mayff(i)
        end do
        !   fft to calculate ckfm and ckff
        ckfm = fst(crfm,1)
        ckff = fst(crff,1)
        ! judge the convergence and exit the cycle
        lambda = judge(test,ckfm)
        print *,"lambda = ",lambda
        print *,"times = ",times
        if(lambda < error)exit
        !times = times + 1
        !if (mod(times,1000) == 0)then
        !    print *,"lambda = ",lambda
        !    pause
        !endif
    end do
end subroutine evolution
!}}}
!function judge{{{
function judge(a,b)
    real(8),intent(in)   :: a(n)
    real(8),intent(in)   :: b(n)
    real(8)              :: judge
    
    judge = 0.0
    do ii = 1,n
        judge = judge + abs(a(ii) - b(ii))
    end do
    judge = judge/dble(n)
end function judge
!}}}
end module module_evolution
