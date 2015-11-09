program main
    implicit none
    real(8),parameter  :: pre   = 2.0/sqrt(3.141592653)
    integer,parameter  :: total = int(1E4)
    real(8),parameter  :: pi    = 3.141592653
    real(8),parameter  :: delta = pi/dble(total)
    real(8),parameter  :: dt    = 0.01
    real(8)            :: x         ! the up-limit
    integer            :: m   ! for the summation cycle
    integer            :: i   ! for the main cycle
    integer            :: j   ! for the main cycle
    real               :: t1
    real               :: t2
    character(20)      :: filename

    !do i = 0,10
    !    write(filename,"('bessel',i2.2,'.txt')")i
    !    open(10,file = filename)
    !    call cpu_time(t1)
    !    do j = 1,1000
    !        x = j*dt
    !        write(10,"(2f12.6)")x,besjn(x,i)
    !    end do
    !    call cpu_time(t2)
    !    close(10)
    !    print *,"time cost is ",t2 - t1
    !end do
    filename = "rdata.txt"
    open(10,file = filename)
    do i = 1,100
        do j = 1,100
            x = sqrt(dble(i**2.0 + j**2.0))*dt
            write(10,"(3f12.6)")i*dt,j*dt,besjn(x,0)
        end do
    end do
    close(10)
    contains
!{{{
function besjn(x,n)
    integer,intent(in)            :: n
    real(8),intent(in)            :: x 
    real(8)                       :: besjn
    real(8)                       :: tau
    real(8)                       :: tmp
    integer                       :: ii

    
    tmp = 0.0
    do ii = 1,total
        tau = ii*delta
        tmp = tmp + delta*cos(n*tau - x*sin(tau))
    end do
    besjn = tmp/pi
end function besjn
!}}}
!{{{
!}}}
end program main
