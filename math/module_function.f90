module module_function
    use module_common
    use module_algebra

    contains
!subroutine iter_NR{{{
subroutine iter_NR()
    real(8)               :: x1
    real(8)               :: y1
    real(8)               :: x2
    real(8)               :: y2
        do 
            x1 = x(1)
            y1 = x(2)
            x2 = x(3)
            y2 = x(4)

            a(1,1) = 2.0*x1
            a(1,2) = 4.0*y1
            a(1,3) = 0.0 
            a(1,4) = 0.0 
            a(2,1) = 0.0 
            a(2,2) = 0.0 
            a(2,3) = 6.0*x2 - 2.0*y2 - 20.0
            a(2,4) = -2.0*x2 - 6.0*y2 - 20.0
            a(3,1) = 4.0
            a(3,2) = - 1.0
            a(3,3) = - 4.0
            a(3,4) = 1.0
            a(4,1) = x2 + 3.0*y2 + 10.0 
            a(4,2) = 3.0*x2 - y2 - 10.0
            a(4,3) = -2.0*x2 - (3.0*y2 + 10.0)&
                     + 3.0*(y1 - y2)
            a(4,4) = 3.0*(x1 - x2) + 2.0*y2 &
                     + 10.0 - 3.0*x2 - y1
            b(1)   = fun1(x)
            b(2)   = fun2(x)
            b(3)   = fun3(x)
            b(4)   = fun4(x)
            deltax = sol_equ(a,b,n)
            lambda = asum(deltax)
            if(lambda < error)exit
            x = x - deltax
            times = times + 1
            if(times == fre)exit
        end do
end subroutine iter_NR
!}}}
!function fun1{{{
function fun1(x)
    real(8),intent(in)         :: x(n)
    real(8)                    :: fun1
    real(8)                    :: y1
    real(8)                    :: y2
    real(8)                    :: x1
    real(8)                    :: x2

    x1 = x(1)
    y1 = x(2)
    x2 = x(3)
    y2 = x(4)

    fun1 = x1**2.0 + 4.0*y1**2.0 - 4.0
end function fun1
!}}}
!function fun2{{{
function fun2(x)
    real(8),intent(in)         :: x(n)
    real(8)                    :: fun2
    real(8)                    :: y1
    real(8)                    :: y2
    real(8)                    :: x1
    real(8)                    :: x2

    x1 = x(1)
    y1 = x(2)
    x2 = x(3)
    y2 = x(4)

    fun2 = 3.0*x2**2.0 - 2.0*x2*y2 - 3.0*y2**2.0 &
           - 20.0*x2 - 20.0*y2 + 99.0
end function fun2
!}}}
!function fun3{{{
function fun3(x)
    real(8),intent(in)         :: x(n)
    real(8)                    :: fun3
    real(8)                    :: y1
    real(8)                    :: y2
    real(8)                    :: x1
    real(8)                    :: x2

    x1 = x(1)
    y1 = x(2)
    x2 = x(3)
    y2 = x(4)

    fun3 = 4.0*x1 - y1 - 4.0*x2 + y2
end function fun3
!}}}
!function fun4{{{
function fun4(x)
    real(8),intent(in)         :: x(n)
    real(8)                    :: fun4
    real(8)                    :: y1
    real(8)                    :: y2
    real(8)                    :: x1
    real(8)                    :: x2

    x1 = x(1)
    y1 = x(2)
    x2 = x(3)
    y2 = x(4)

    fun4 = (x1 - x2)*(x2 + 3.0*y2 + 10.0)&
           + (y1 - y2)*(3.0*x2 - y2 - 10.0)
end function fun4
!}}}
!function asum{{{
function asum(x)
    real(8),intent(in)     :: x(n)
    real(8)                :: asum
    integer                :: ii

    asum = 0.0
    do ii = 1,n
        asum = asum + abs(x(ii))
    end do
end function asum
!}}}
!subroutine power{{{
subroutine power(xnew,ynew,x,y,z)
    real,intent(in)          :: x   
    real,intent(in)          :: y   
    real,intent(in)          :: z   
    real,intent(out)         :: xnew   
    real,intent(out)         :: ynew
    real                     :: r
    real                     :: theta

    r = sqrt(x**2.0 + y**2.0)
    theta = atan(y/x)
    xnew = r**z*cos(theta*z)
    ynew = r**z*sin(theta*z)
end subroutine power
!}}}
end module module_function
