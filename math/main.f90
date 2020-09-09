program main
    implicit none
    real(8),parameter  :: pre   = 2.0/sqrt(3.141592653)
    integer,parameter  :: total = int(1E6)
    real(8)            :: jalpha    ! Bessel function
    real(8)            :: alpha     ! the rank of jalpha
    real(8)            :: x         ! the up-limit
    integer            :: m   ! for the summation cycle
    integer            :: i   ! for the main cycle
    real(8)            :: tmp
    real               :: t1
    real               :: t2
    character(20)      :: filename

    filename = "erf.txt"
    open(10,file = filename)
    do i = 1,300
        x = - 1.5 + i/100.0
        write(10,*)x,erf(x)
    end do
    close(10)
    x = 2.0
    tmp = gamma(35.0)
    print "(E20.4)",tmp
    print *,"j0",j0(1.0)
    !print *,erf(2.1)
    !print *,erf(3.0/sqrt(2.0))
    !x = 3.0/sqrt(2.0)
    !print *,"my function",error(x)
    !print *,erf(sqrt(2.0))
    !print *,gamma(0.5)

    !do i = 1,10
    !    if(i == 2)cycle
    !    alpha = dble(i)
    !    write(filename,"('data',i2.2,'.txt')")i 
    !    open(10,file = filename)
    !    do m = 1,1000
    !        x   = m/100.0
    !        tmp = p_gamma(x,alpha)
    !        write(10,*)x,tmp
    !    end do
    !    close(10)
    !    call cpu_time(t2) 
    !    print *,alpha,"time cost is",t2 - t1
    !    t1 = t2
    !end do

    contains
!function igamma{{{
function igamma(x)
    real(8),intent(in)   :: x
    real(8)              :: igamma
    integer              :: ii

    igamma = 0.0
end function igamma
!}}}
! function p_gamma{{{
! incomplete function gamma
function p_gamma(x,alpha)
    real(8),intent(in)   :: x
    real(8),intent(in)   :: alpha
    real(8)              :: p_gamma
    real(8)              :: tmp
    real(8)              :: deltat
    real(8)              :: t ! independent variable
    integer              :: ii

    tmp = 0.0
    deltat = x/dble(total)
    do ii = 1,total
         t = ii*deltat
       tmp = tmp + deltat*(t**(alpha - 1.0)*exp(- t)) 
    end do
    p_gamma = tmp
end function p_gamma
!}}}
!function error{{{
!gauss error function
function error(x)
    real(8),intent(in)   :: x
    real(8)              :: error
    real(8)              :: tmp  
    real(8)              :: eta  ! independent variable 
    real(8)              :: deltat
    integer              :: ii

    tmp    = 0.0
    deltat = x/dble(total)
    do  ii = 1,total
        eta = deltat*ii
        tmp = tmp + deltat*exp(- eta**2.0) 
    end do
    error = tmp*pre 
end function error
!}}}
!{{{
!}}}
!{{{
!}}}
end program main
