module module_evolution
      use module_fst
      use module_algebra
      use module_judge
      implicit none
      contains
!subroutine evolution{{{
subroutine evolution()
    integer              :: ii
    integer              :: jj
    real(8)              :: tmp_array(n)
    ! give an initial value of ck
    do i = 1,m
        do j = i,m
                ck(i,j,:) = 1.0
        end do
        do j = 1,m - 1
                ck(i,j,:) = ck(j,i,:)
        end do
    end do
    
    !  get delta function
    do i = 1,m
        do j = i,m
            if(i == j)then
                deltafun(i,j) = 1.0
            else
                deltafun(i,j) = 0.0
            endif
        end do
    end do
! OZ equation
    times = 0
    rate  = 0.9
!solve{{{
    do 
        times = times + 1
        if(times > 100)exit
        do i = 1,m
            do j = i,m
                test(i,j,:) = ck(i,j,:)
            end do
        end do
        ! get gamma_k with OZ equation
        do i = 1,m
            do j = i,m
                    gamma_k(i,j,1) = 0.0
            end do
        end do
        do ii = 2,n
            ! get the value of array a and b 
            do i = 1,m
                do j = 1,m 
                    a(i,j) = deltafun(i,j)*dk(ii) - ck(i,j,ii)*rho(j)
                    b(i,j) = dk(ii)*ck(i,j,ii)
                end do
            end do
            ! get hk
            hk(:,:,ii)     = sol_mat(a,b,m)
            ! calculate gamma_k
            do i = 1,m
                do j = i,m
                    gamma_k(i,j,ii) = hk(i,j,ii) - ck(i,j,ii)
                end do
            end do
        end do   ! to get the new gamma_k
        !******************************************************
        ! iteration and judge convergence
        ! calculate gamma_r with inverse fft
        do i = 1,m
            do j = i,m
                do jj = 1,n
                    tmp_array(jj) = gamma_k(i,j,jj)
                end do
                gamma_r(i,j,:) = fst(tmp_array,- 1)
            end do
        end do
        ! calculate cr with PY
        do i = 1,m
            do j = i,m
                cr(i,j,:) = (dr(:) + gamma_r(i,j,:))*mayfun(i,j,:)
            end do
        end do
        ! calculate ck with fft
        do i = 1,m
            do j = i,m
                ck(i,j,:) = fst(cr(i,j,:),1)
            end do
            do j = 1,i-1
                ck(i,j,:) = ck(j,i,:)
            end do
        end do
    !    ! judge convergence
    !    lambda = 0.0
    !    do i = 1,m
    !        do j = i,m
    !            lambda1 = setion_rate(test(i,j,:),ck(i,j,:),rate)
    !             lambda = lambda + lambda1 
    !        end do
    !    end do
    !    print *,"times = ",times,"lambda = ",lambda
    !    !print *,times,lambda
    !    if(mod(times,fre) == 0)then
    !        print *,rate, "lambda = ",lambda
    !    endif
    !    if(lambda < error)exit
    !    times = times + 1
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
