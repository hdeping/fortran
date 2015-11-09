program main
    use module_algebra
    implicit none
    integer,parameter        :: n     = 2
    integer,parameter        :: m     = int(1E8)
    integer,parameter        :: fre   = 100
    real(8),parameter        :: delta = 1E-2
    real(8),parameter        :: error = 1E-6
    integer                  :: i
    !integer                  :: j
    !integer                  :: k
    integer                  :: times
    real(8)                  :: x1
    real(8)                  :: x2
    real(8)                  :: lambda
    real(8)                  :: x(n)
    real(8)                  :: deltax(n)
    real(8)                  :: a(n,n)
    real(8)                  :: b(n)
    character(20)            :: filename

    filename = "solution.txt"
    open(10,file = filename)

    call random_seed()
    do i = 1,100
        call random_number(x1) 
        call random_number(x2) 
        x1 = x1*1000 - 500
        x2 = x2*2000 - 500
        x  = (/x1,x2/)
        times = 0
        do 
            a(1,1) = 2.0*x(1)
            a(1,2) = 6.0*x(2) 
            a(2,1) = 1.0 + 5.0*x(2) 
            a(2,2) = - 2.0 + 5.0*x(1) 
            b(1)   = fun1(x)
            b(2)   = fun2(x)
            deltax = sol_equ(a,b,n)
            lambda = asum(deltax)
            if(lambda < error)exit
            x = x - deltax
            times = times + 1
            if(mod(times,fre) == 0)then
                print *,x
                print *,"lambda = ",lambda
            endif
        end do
        write(10,"(2f12.6)")x(:)
    end do

    close(10)

    contains
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
!function fun1{{{
function fun1(x)
    real(8),intent(in)        :: x(n)
    real(8)                   :: fun1

    fun1 = x(1)**2.0 + 3.0*x(2)**2.0 - 10.0
end function fun1
!}}}
!function fun2{{{
function fun2(x)
    real(8),intent(in)        :: x(n)
    real(8)                   :: fun2

    fun2 = x(1) - 2.0*x(2) + &
           5.0*x(1)*x(2) - 9.0
end function fun2
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



end program main
