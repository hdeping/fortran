module module_common
    implicit none
    integer,parameter            :: n     = 100
    integer,parameter            :: total = 2000
    real(8),parameter            :: g     = 6.61E-11
    real(8),parameter            :: pi    = 3.141592653
    integer                      :: i
    integer                      :: j
    integer                      :: k
    character(10)                :: filename 
    contains
!function getintegral{{{
function getintegral(a,b)
    real(8),intent(in)              :: a(3)
    real(8),intent(in)              :: b(3)
    real(8)                         :: getintegral(3)
    real(8)                         :: d
    real(8)                         :: dr
    real(8)                         :: dphi
    real(8)                         :: dtheta
    real(8)                         :: r
    real(8)                         :: phi
    real(8)                         :: theta
    real(8)                         :: x
    real(8)                         :: y
    real(8)                         :: z
    real(8)                         :: tmp
    integer                         :: ii
    integer                         :: jj
    integer                         :: kk
    integer                         :: itmp

    d       = 2.0
    dr      = d/dble(total)
    dtheta  = pi/dble(total)
    dphi    = 2*pi/dble(total)
    getintegral = 0.0
    do ii = 1,total  ! r
        do jj = 1,total  ! theta
            do kk = 1,total   ! phi
                r     = dr*dble(ii)
                theta = dtheta*dble(jj)
                phi   = dphi*dble(kk)
                x     = r*sin(theta)*cos(phi) + a(1) - b(1)
                y     = r*sin(theta)*sin(phi) + a(2) - b(2)
                z     = r*cos(theta) + a(3) - b(3)
                tmp   = r**2*sin(theta)*dr*dtheta*dphi/&
                        sqrt(x**2.0 + y**2.0 + z**2.0)**1.5
                getintegral(1) = getintegral(1) - tmp*x
                getintegral(2) = getintegral(2) - tmp*y
                getintegral(3) = getintegral(3) - tmp*z
                 
            enddo !cycle ends
        enddo !cycle ends
    enddo !cycle ends
    getintegral = getintegral/(4.0*pi*d**3.0/3.0)
end function getintegral
!}}}
end module module_common
