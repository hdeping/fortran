module module_common
    implicit none
    integer,parameter            :: n  = 300
    integer,parameter            :: m  = 5 ! rank of polynomial
    real(8),parameter            :: pi = 3.141592653
    real(8),parameter            :: dt = 0.03
    integer                      :: i
    integer                      :: j
    integer                      :: k
    character(10)                :: filename 
    real(8)                      :: x
    real(8)                      :: y
    real(8)                      :: r
    real(8)                      :: theta
    real(8)                      :: tmp
    real(8)                      :: a(m+1) ! coefficient
    real(8)                      :: t1
    real(8)                      :: t2

    contains
!subroutine getr_theta{{{
subroutine getr_theta(r,theta,x,y)
    real(8),intent(in)          :: x
    real(8),intent(in)          :: y
    real(8),intent(out)         :: r
    real(8),intent(out)         :: theta
    integer                     :: ii

    r = sqrt(x**2.0 + y**2.0)
    !  first judge about theta
    if(x == 0.0)then
        if(y > 0.0)then
            theta = 0.0
        else
            theta = pi
        endif
    else
        theta  = atan(y/x) 
    endif
    
    !  second judge about theta
    if(x < 0.0)then
        if(y > 0.0)then
            theta = theta + pi/2.0
        else
            theta = theta - pi/2.0
        endif
    endif
end subroutine getr_theta
!}}}
!function getNorm{{{
function getNorm(x,y,a)
    real(8),intent(in)          :: a(m+1)
    real(8),intent(in)          :: x
    real(8),intent(in)          :: y
    real(8)                     :: getNorm
    real(8)                     :: tmpx
    real(8)                     :: tmpy
    integer                     :: ii
    call getr_theta(r,theta,x,y)
    tmpx    = getg(r,theta,a)
    tmpy    = geth(r,theta,a)
    getNorm = sqrt(tmpx**2.0 + tmpy**2.0)

end function getNorm
!}}}
!function getPoly{{{
function getPoly(x,a,dima)
    real(8),intent(in)          :: x
    integer,intent(in)          :: dima
    real(8),intent(in)          :: a(dima)
    real(8)                     :: getPoly
    integer                     :: ii

    getPoly = 0.0
    do ii = dima,1,- 1
        getPoly = getPoly*x + a(ii)
    end do
end function getPoly
!}}}
!function getg{{{
function getg(r,theta,a)
    real(8),intent(in)          :: r
    real(8),intent(in)          :: theta
    real(8),intent(in)          :: a(m+1)
    real(8)                     :: b(m+1)
    real(8)                     :: getg
    integer                     :: ii

    ! get the coefficients of g(r,theta)
    do ii = 1,m+1
         b(ii) = a(ii)*cos((ii - 1)*theta) 
         !print *,theta,b(ii)
         !pause
    end do
    getg = getPoly(r,b,m+1)

end function getg
!}}}
!function geth{{{
function geth(r,theta,a)
    real(8),intent(in)          :: r
    real(8),intent(in)          :: theta
    real(8),intent(in)          :: a(m+1)
    real(8)                     :: b(m+1)
    real(8)                     :: geth
    integer                     :: ii

    ! get the coefficients of g(r,theta)
    do ii = 1,m+1
         b(ii) = a(ii)*sin((ii - 1)*theta) 
    end do
    geth = getPoly(r,b,m+1)

end function geth
!}}}

end module module_common
