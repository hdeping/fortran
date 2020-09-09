module module_common
    implicit none
    integer,parameter            :: n = 100
    integer,parameter            :: m = 100
    integer                      :: i
    integer                      :: j
    integer                      :: k
    real(8)                      :: x(m)
    real(8)                      :: y(n,m)
    character(20)                :: filename 

    contains
!function legendre{{{
function legendre(x)
    real(8),intent(in)    :: x(m)
    real(8)               :: legendre(n,m)
    integer               :: ii
    integer               :: jj
    integer               :: kk

    do ii = 1,m
        legendre(1,ii) = 1.0
        legendre(2,ii) = x(ii)
    enddo !cycle ends

    do ii = 2,n - 1
        do jj = 1,m
            legendre(ii + 1,jj) =  legendre(ii,jj)*dble((2*jj+1))&
                *x(jj)/dble((jj+1)) - legendre(ii - 1,jj)*dble(jj)&
                /dble((jj+1))
        enddo !cycle ends
    enddo !cycle ends
     
     
end function legendre
!}}}
!function laguerre{{{
function laguerre(x)
    real(8),intent(in)    :: x(m)
    real(8)               :: laguerre(n,m)
    integer               :: ii
    integer               :: jj
    integer               :: kk

    do ii = 1,m
        laguerre(1,ii) = 1.0
        laguerre(2,ii) = 1.0 - x(ii)
    enddo !cycle ends

    do ii = 2,n - 1
        do jj = 1,m
            laguerre(ii + 1,jj) =  laguerre(ii,jj)*dble((2*jj+1 - x(jj)))&
                /dble((jj+1)) - laguerre(ii - 1,jj)*dble(jj)&
                /dble((jj+1))
        enddo !cycle ends
    enddo !cycle ends
     
     
end function laguerre
!}}}
end module module_common
