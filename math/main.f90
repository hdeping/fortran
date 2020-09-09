program main
    implicit none
    integer,parameter   :: n     = 49
    integer,parameter   :: num   = 14
    real(8),parameter   :: delta = 1D-9
    real(8),parameter   :: fn    = 5D-1
    real(8)             :: bg   
    real(8)             :: ed   
    real(8)             :: dt   
    real(8)             :: reta ! restart
    integer             :: i
    integer             :: j
    real(8)             :: x
    real(8)             :: a(n)
    real(8)             :: b(n)
    character(10)       :: filename 

    filename = "data.txt"
    open(10,file = filename)

    bg   = 0.0
    ed   = 5.0
    dt   = 1D-2
    do i = 1,n
         x   = dble(i)*dt
        b(i) = sol_b(x)
        write(10,"(2f18.9)")x,b(i)
    end do
    reta = dble(n)*dt
    do j = 1,num
        dt   = dt/10.0
        do i = 1,9
             x   = reta + dble(i)*dt
            b(i) = sol_b(x)
             x   = - log(fn - x)
            write(10,"(2f18.9)")x,b(i)
        end do
        reta  = x
    end do

    close(10)
    contains
!function getgamma{{{
function getgamma(x)
    real(8),intent(in)       :: x
    real(8)                  ::getgamma
    getgamma = log(1 - 2.0*x) + 2.0*x
end function getgamma
!}}}
!function getgamma2{{{
function getgamma2(x)
    real(8),intent(in)       :: x
    real(8)                  ::getgamma2
    getgamma2 = log(1 + 2.0*x) - 2.0*x
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
!function par_f{{{
function par_f(x)
    real(8),intent(in)       :: x
    real(8)                  ::par_f
    par_f = 2.0/(1.0 + 2.0*x) - 2.0*x
end function par_f
!}}}
!function sol_b{{{
function sol_b(x)
    real(8),intent(in)   :: x
    real(8)              :: sol_b
    real(8)              :: y
    real(8)              :: ra
    real(8)              :: rb
    real(8)              :: rc

    y  = getgamma(x)
    ra = bg
    !print *,"x = ",x
    do 
        if(fun(ed,y) < 0.0)exit
        ed = 2.0*ed  
    end do
    rb = ed

    do 
        rc  = (ra + rb)/2.0
        if(fun(rc,y) > 0.0)then
            ra = rc
        else
            rb = rc
        endif
        if(abs(ra - rb) < delta)exit
    end do
    sol_b = rc
    !pause
    !print *,"times = ",times


end function sol_b
!}}}

end program main
