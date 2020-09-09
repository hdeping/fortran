program main
    implicit none
    integer,parameter   :: n    = 499
    integer,parameter   :: m    = 1200
    real(8),parameter   :: dt   = 1D-3
    real(8),parameter   :: dt1  = 1D-2
    integer             :: i
    integer             :: j
    real(8)             :: x
    real(8)             :: a(n)
    real(8)             :: b(m)
    character(10)       :: filename 

    filename = "data.txt"
    open(10,file = filename)

    do i = 1,n
        x    = dble(i)*dt
        a(i) = getgamma(x)
        !b(i) = getgamma2(x)
        !write(10,"(3f10.5)")x,a(i),b(i)
        write(10,"(2f10.5)")x,a(i)
    end do
    close(10)

    filename = "data1.txt"
    open(10,file = filename)
    do i = 1,m
        x    = dble(i)*dt1
        b(i) = getgamma2(x)
        write(10,"(2f10.5)")x,b(i)
    end do
    close(10)
    contains
!function getgamma{{{
function getgamma(x)
    real(8),intent(in)       :: x
    real(8)                  ::getgamma
    getgamma = gamma(1 - x)**2.0/gamma(1 - 2.0*x)
end function getgamma
!}}}
!function getgamma2{{{
function getgamma2(x)
    real(8),intent(in)       :: x
    real(8)                  ::getgamma2
    getgamma2 = gamma(1 + x)**2.0/gamma(1 + 2.0*x)
end function getgamma2
!}}}
end program main
