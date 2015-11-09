module module_common
    implicit none
    integer,parameter            :: n     = 200
    integer,parameter            :: total = int(1E5)
    real(8),parameter            :: pi    = 3.141592653
    real(8),parameter            :: delta = pi/2/dble(total)
    integer                      :: i
    integer                      :: j
    integer                      :: k
    real(8)                      :: t1
    real(8)                      :: t2
    character(10)                :: filename 

    contains
!function getlength{{{
function getlength(a,b)
    real(8),intent(in)        :: a
    real(8),intent(in)        :: b
    real(8)                   :: getlength
    real(8)                   :: x
    integer                   :: ii

    getlength = 0D0
    do ii = 1,total
        x    = dble(ii)*delta
        getlength = getlength + delta*sqrt(a**2.0*cos(x)**2.0&
                    +  b**2.0*sin(x)**2.0)
    enddo !cycle ends


end function getlength
!}}}
end module module_common
