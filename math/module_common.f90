module module_common
    implicit none
    integer,parameter            :: n = 100
    integer                      :: i
    integer                      :: j
    integer                      :: k
    real(8)                      :: x(n)
    real(8)                      :: y(n)
    real(8)                      :: val
    character(10)                :: filename 
!*******************************************
    contains
!function ave_sum{{{
! get the average of an array
function ave_sum(x)
    real(8),intent(in)           :: x(n)
    real(8)                      :: ave_sum
    integer                      :: ii
    real(8)                      :: tmp

    do ii = 1,n
        if(ii == 1)then
            tmp = x(ii)
        else
            tmp = (tmp*dble((ii - 1)) + x(ii))/dble(ii)
    end do
    ave_sum = tmp

end function ave_sum
!}}}
!function ave_mul{{{
! get the average of an array
function ave_mul(x)
    real(8),intent(in)           :: x(n)
    real(8)                      :: ave_mul
    real(8)                      :: y(n)
    integer                      :: ii
    real(8)                      :: tmp

    do ii = 1,n
        y(ii) = log(x(ii))
    end do
    ave_mul = ave_sum(y)
    ave_mul = exp(ave_mul)
end function ave_mul
!}}}
end module module_common
