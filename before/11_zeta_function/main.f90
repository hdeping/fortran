program main
    implicit none
    integer,parameter     :: n = int(1E6)
    integer               :: i
    integer               :: j
    real                  :: a
    real                  :: b
    real                  :: t1
    real                  :: t2
    real(8)               :: tmp1
    real(8)               :: tmp2
    real(8)               :: sa(5)
    real(8)               :: sb(5)
    character(12)         :: filename
    !  tmp1 for the real    part
    !  tmp2 for the imagine part
    !tmp1 = zeta_re(1.0,2.0)
    !tmp2 = zeta_im(1.0,2.0)
    !print *,tmp1,tmp2
    sa = (/2.463161, 1.286496, 2.307570 &
         , 1.382763, 0.964685/) 
    sb = (/23.298320, 31.708250, 38.489983,&
         42.290964, 48.847159/)
    call cpu_time(t1)
    do i = 1,5
        a = sa(i)
        b = sb(i)
        tmp1 = zeta_re(a,b)
        tmp2 = zeta_im(a,b)
        print "('a,b,re,im ==> ',4f18.9)",a,b,tmp1,tmp2
    end do
    print *,"gamma ==> ",gamma(0.6)
    !do i = 1,10
    !    a = dble(i)/10.0
    !    write(filename,"('data',i2.2,'.txt')")i
    !    open(10,file = filename)
    !    do j = 1,10000
    !        b = 1.0 + dble(j)/100.0
    !        tmp1 = zeta_re(a,b)
    !        tmp2 = zeta_im(a,b)
    !        write(10,"(3f12.6)")b,tmp1,tmp2
    !    end do
    !    close(10)
    !end do
    !call cpu_time(t2)
    !filename = "time.txt"
    !open(10,file = filename)
    !write(10,*)"time cost is ==>",t2 - t1
    !close(10)
    contains
!function zeta_re{{{
    function zeta_re(a,b)
        real                  :: a
        real                  :: b
        real(8)               :: zeta_re
        real(8)               :: tmp
        integer               :: k

        tmp = 0
        do k = 1,n
            tmp = tmp + exp(- a*log(dble(k)))*&
                   cos(b*log(dble(k))) 
        end do
        zeta_re = tmp
        
    end function zeta_re
!}}}
!function zeta_im{{{
    function zeta_im(a,b)
        real                  :: a
        real                  :: b
        real(8)               :: zeta_im
        real(8)               :: tmp
        integer               :: k

        tmp = 0
        do k = 1,n
            tmp = tmp - exp(- a*log(dble(k)))*&
                   sin(b*log(dble(k))) 
        end do
        zeta_im = tmp
        
    end function zeta_im
!}}}
end program main
