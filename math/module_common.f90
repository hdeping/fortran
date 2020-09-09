module module_common
    implicit none
!variables{{{
    integer,parameter            :: s     = 3
    integer,parameter            :: n     = 3*10**s
    integer,parameter            :: m     = 3*10**s
    integer,parameter            :: nnum  = 10**(s - 2)
    integer,parameter            :: mnum  = 10**(s - 2)
    integer,parameter            :: nfre  = n/nnum
    integer,parameter            :: mfre  = m/mnum
    real(8),parameter            :: dt    = 3.0/dble(n)
    real(8)                      :: f(nfre,mfre)
    real(8)                      :: a1(n)
    real(8)                      :: a(n)
    real(8)                      :: b(m)
    real(8)                      :: tmp
    !real(8)                      :: fx(n,m)
    !real(8)                      :: fy(n,m)
    real(8)                      :: x(n)
    real(8)                      :: y(m)
    real(8)                      :: t1
    real(8)                      :: t2
    integer                      :: i
    integer                      :: j
    integer                      :: itmp 
    integer                      :: jtmp
    integer                      :: k
    character(10)                :: filename 
!}}}

    contains
!function obj{{{
function obj(x,y)
    real(8),intent(in)         :: x
    real(8),intent(in)         :: y
    real(8)                    :: obj

    obj = x*sin(y) + 5.0*y*cos(x)
end function obj
!}}}
end module module_common
