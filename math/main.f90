program main
    use module_regre
    implicit none
    integer,parameter   :: n     = 998
    integer,parameter   :: m     = n/2
    real(8),parameter   :: dt1   = 1D-3
    real(8),parameter   :: dt    = 1D-2
    real(8),parameter   :: delta = 1D-7
    real(8),parameter   :: bg    = 0.0
    real(8),parameter   :: ed    = 5.0
    integer             :: i
    integer             :: j
    real(8)             :: x
    real(8)             :: y
    real(8)             :: r
    real(8)             :: tmp
    real(8)             :: a(n)
    real(8)             :: b(n)
    character(10)       :: filename 

    filename = "data.txt"
    open(10,file = filename)

    do i = m,1,- 1
        x    = - dble(i)*dt1
        a(i) = x
        b(i) = getgamma2(x)
    end do
    do i = 1,m
        x        = dble(i)*dt
        a(i+m) = x
        b(i + m) = getgamma2(x)
    end do

    tmp  = 0.0
    do i = 1,n
        a(i) = log(2.0*a(i) + 1) - 2.0*a(i)
        b(i) = log(b(i))
        tmp  = (tmp*(i - 1) + b(i)/a(i))/dble(i)
        write(10,"(3f12.6)")a(i),b(i),tmp
    end do
    print *,"tmp = ",tmp
    call regre(x,y,r,a(1:m),b(1:m),m)
    print *,"x = ",x
    print *,"y = ",y
    print *,"r = ",r
    call regre(x,y,r,a(m+1:n),b(m+1:n),m)
    print *,"x = ",x
    print *,"y = ",y
    print *,"r = ",r


    close(10)
end program main
