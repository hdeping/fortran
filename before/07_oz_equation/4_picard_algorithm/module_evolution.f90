module module_evolution
      use module_fst
      use module_algebra
      use module_judge
      implicit none
      contains
!subroutine evolution{{{
subroutine evolution()
    integer              :: ii
    real(8)              :: tmp_array(n)
    real(8)              :: a(3,3)
    real(8)              :: b(3)
    real(8)              :: hk(3)
    ! give an initial value of ck
    ck11  = 1.0
    ck12  = 1.0
    ck22  = 1.0
    rate  = 0.1
    times = 0
!solve{{{
    do 
        ! get previous ck
        test1  = ck11  
        test2  = ck12  
        test3  = ck22  
        ! get gamma_k with OZ equation
        gamma_k11 = 1.0
        gamma_k12 = 1.0
        gamma_k22 = 1.0
        do ii = 2,n
            ! get the value of array a and b 
            a(1,1) = dk(ii) - rho1*ck11(ii)
            a(1,2) = - rho2*ck12(ii)
            a(1,3) = 0.0 
            a(2,1) = - rho1*ck12(ii)
            a(2,2) = dk(ii) - rho2*ck22(ii)
            a(2,3) = 0.0
            a(3,1) = 0.0
            a(3,2) = - rho1*ck12(ii)
            a(3,3) = dk(ii) - rho2*ck22(ii)
            b(1)   = dk(ii)*ck11(ii)
            b(2)   = dk(ii)*ck12(ii)
            b(3)   = dk(ii)*ck22(ii)
            ! get hk
            hk   = sol_equ(a,b,3)
            hk11 = hk(1)
            hk12 = hk(2)
            hk22 = hk(3)
            ! calculate gamma_k
            gamma_k11(ii) = hk11(ii) - ck11(ii)
            gamma_k12(ii) = hk12(ii) - ck12(ii)
            gamma_k22(ii) = hk22(ii) - ck22(ii)

        end do   ! to get the new gamma_k
        !******************************************************
        ! iteration and judge convergence
        !******************************************************
        ! calculate gamma_r with inverse fft
        gamma_r11 = fst(gamma_k11,- 1)
        gamma_r12 = fst(gamma_k12,- 1)
        gamma_r22 = fst(gamma_k22,- 1)
        ! calculate cr with PY
        do ii = 1,n
            cr11(ii) = (gamma_r11(ii) + dr(ii))*may11(ii)
            cr12(ii) = (gamma_r12(ii) + dr(ii))*may12(ii)
            cr22(ii) = (gamma_r22(ii) + dr(ii))*may22(ii)
        end do
        ! calculate ck with fft
        ck11 = fst(cr11,1)
        ck12 = fst(cr12,1)
        ck22 = fst(cr22,1)
        ! judge convergence
        lambda = 0.0
        lambda1 = setion_rate(test1,ck11,rate)
        lambda = lambda + lambda1
        lambda1 = setion_rate(test2,ck12,rate)
        lambda = lambda + lambda1
        lambda1 = setion_rate(test3,ck22,rate)
        lambda = lambda + lambda1
        !print *,"times = ",times,"lambda = ",lambda
        !print *,times,lambda
        if(mod(times,fre) == 0)then
            print *,rate, "lambda = ",lambda
        endif
        if(lambda < error)exit
        times = times + 1
    end do  ! end solving
!}}}
    print *,"lambda = ",lambda
end subroutine evolution
!}}}
!function may{{{
function may(diameter)
    real(8),intent(in)         :: diameter
    real(8)                    :: may(n)
    integer                    :: ii
    do ii = 1,n
        if(dr(ii) < diameter)then
             may(ii) = - 1
         else
             may(ii) = 0
         endif
    end do
end function may
!}}}
end module module_evolution
