module module_judge
    use module_common
    contains
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
    real(8),intent(in)      :: rate
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
end module module_judge
