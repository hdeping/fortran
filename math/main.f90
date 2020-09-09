program main
    implicit none
    integer,parameter   :: n     = 4999
    real(8),parameter   :: dt    = 1D-4
    real(8),parameter   :: delta = 1D-7
    real(8),parameter   :: bg    = 0.0
    real(8),parameter   :: ed    = 5.0
    integer             :: i
    integer             :: j
    real(8)             :: x
    real(8)             :: a(n)
    real(8)             :: b(n)
    character(10)       :: filename 

    filename = "data.txt"
    open(10,file = filename)

    do i = 1,n
         x   = dble(i)*dt
        b(i) = sol_b(x)
        !write(10,"(3f10.5)")x,a(i),b(i)
        write(10,"(4f12.5)")x,b(i),1.0/b(i),log(b(i))
        !print *,getgamma(x) - getgamma2(b(i))
        !pause
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
!function fun{{{
function fun(x,y)
    real(8),intent(in)       :: x
    real(8),intent(in)       :: y
    real(8)                  ::fun
    fun = getgamma2(x) - y
end function fun
!}}}
!function sol_b{{{
function sol_b(x)
    real(8),intent(in)   :: x
    real(8)              :: sol_b
    real(8)              :: y
    real(8)              :: ra
    real(8)              :: rb
    real(8)              :: rc
    integer              :: times 

    ra    = bg
    rb    = ed
    !if(x > 0.499)then
    !    ra = ed
    !endif
    times = 0
    y     = getgamma(x)

    do 
        rc = (ra + rb)/2.0 
        if(fun(rc,y) > 0.0)then
            ra = rc
        else
            rb = rc
        endif
        !  convergence
        if(abs(ra - rb) < delta)exit
    end do
    sol_b = rc


end function sol_b
!}}}
end program main
