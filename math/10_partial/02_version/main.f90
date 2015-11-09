program main
    use module_common
    implicit none

    filename = "data.txt"
    open(10,file = filename)
    ! get x and y
    do i = 1,n
        x(i) = i*dt
        y(i) = i*dt
    end do
    ! initial f
    do i = 1,n
        a(i) = exp(x(i))
        a1(i) = cos(y(i))
    end do

    call cpu_time(t1)
    do i = 2,n
        b(1) = a1(i)
        do j = 2,m
           b(j) = (obj(x(i),y(j))*dt - x(i)*a(j)&
                  -y(j)*b(j - 1))/(x(i) + y(j) + 5.0*dt) 
           if((mod(i,nnum) == 0).and.(mod(j,mnum) == 0) )then
               itmp = i/nnum 
               jtmp = j/nnum 
               f(itmp,jtmp) = b(j)
           endif
        end do
        a = b
    end do
    call cpu_time(t2)
    print *,"time cost is ==> ",t2 - t1

    call cpu_time(t1)
    do i = 1,nfre
        do j = 1,mfre
            itmp  =  nnum*i
            jtmp  =  mnum*j
            write(10,"(3f12.6)")x(itmp),y(jtmp),f(i,j)
        end do
        if(i == nfre/2)then
            write(10,*)"    "
        endif
    end do
    call cpu_time(t2)
    print *,"time cost is ==> ",t2 - t1

    close(10)

end program main
