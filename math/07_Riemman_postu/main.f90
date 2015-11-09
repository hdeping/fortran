program main
    use module_common
    implicit none
    integer             :: i 
    integer             :: j 
    character(10)       :: filename 
    
    filename = "data.txt"
    open(10,file = filename)

    call random_seed()
    ! test the functions
    !a   = 1.0
    !b   = 1.0
    !call cpu_time(t1)
    !x   = getzero(a,b)
    !call cpu_time(t2)
    !print *,"time cost is ==> ",t2 - t1
    !print *,"x = ",x(:)
    !print *,"f(a,b) = ",fab(x(1),x(2))
    !print *,"g(a,b) = ",gab(x(1),x(2))
    
    !a = 1.0
    !b = 1.0
    !x = getzero(a,b)
    !print *,"x = ",x(:)
    do i = 1,100000
       call random_number(x1)
       call random_number(x2)
       x1 = 1000*x1 + 1
       x2 = 2*pi*x2
       a  = x1*cos(x2)
       b  = x1*sin(x2)
       x  = getzero(a,b)
       if(lambda > 100)cycle
       write(10,"(2f12.4)")x(:)
    end do

    close(10)
end program main
