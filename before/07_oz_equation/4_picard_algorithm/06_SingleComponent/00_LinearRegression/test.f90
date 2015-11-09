program main
    use module_regre
    use module_common
    implicit none
    real(8)        atmp(n,n)
    real(8)        btmp(n)
    real(8)        btmp_new(n)
    real(8)        res(n)
    real(8)        value_tmp
    
    call random_seed()
    do i = 1,n
        do j = 1,n
            call random_number(x1)
            atmp(i,j) = (i*1.0 + j*1.0)*x1
        end do
        btmp(i) = i**2.0
    end do
!!test the function getbar{{{
!    do i = 1,n
!        value_tmp = getbar(atmp(i,:),n)
!        print *,i,value_tmp
!    end do
!!}}}
!!test the function getmat{{{
!    do i = 1,n
!        do j = 1,n
!           value_tmp = getmat(atmp(i,:),atmp(j,:),n)
!           print *,i,j,value_tmp
!        end do
!    end do
!!}}}
!!test algebra{{{
!    res      = sol_equ(atmp,btmp,n)
!    btmp_new = multi_vec(atmp,res,n,n)
!    do i = 1,n
!        print *,i,res(i),btmp_new(i)
!    end do
!!}}}
!  test the multi_regre
!multi_regre{{{
    b = 3.890
    a = (/5.0, 6.0, 13.3, 18.2/)
    do j = 1,n
        y(j) = b 
        do i = 1,m
            call random_number(x1)
            x(j,i) = dble(j)*x1
            y(j) =  y(j) +  (a(i) + 3.0*x1)*x(j,i)
        end do
    end do

    call multi_regre(a_new,b_new,r,x,y,n,m) 
    print "(18a2)","a_new = ","a"
    do i = 1,m
        print *, a_new(i),a(i)
    end do
    print *,"b_new = ",b_new
    print *,"b = ",b
    print *,"r = ",r
!}}}

    
    
    



end program main
