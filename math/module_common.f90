module module_common
    implicit none
    integer,parameter            :: n = 100
    integer                      :: i
    integer                      :: j
    integer                      :: k
    character(10)                :: filename 
    contains
!function multi_mat{{{
function multi_mat(a,b,n)
    integer,intent(in)       :: n 
    real(8),intent(in)       :: a(n,n)
    real(8),intent(in)       :: b(n,n)
    real(8)                  :: multi_mat(n,n)
    real(8)                  :: tmp
    integer                  :: ii
    integer                  :: jj
    integer                  :: kk

    do ii = 1,n
        do jj = 1,n
            tmp = 0.0
            do kk = 1,n
                tmp = tmp + a(ii,kk)*b(kk,jj)
            end do
            multi_mat(ii,jj) = tmp
        end do
    end do


end function multi_mat
!}}}
!function power_mat{{{
function power_mat(a,n,m)
    integer,intent(in)       :: n 
    integer,intent(in)       :: m
    real(8),intent(in)       :: a(n,n)
    real(8)                  :: power_mat(n,n)
    real(8)                  :: tmp
    integer                  :: ii

    power_mat = a
    do ii = 1,m - 1
        power_mat = multi_mat(power_mat,a,n) 
    end do


end function power_mat
!}}}
end module module_common
