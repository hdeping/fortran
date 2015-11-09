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
        f(i,1) = exp(x(i))
        f(1,i) = cos(y(i))
    end do
    !print *,f(n,1)
    !pause

    call cpu_time(t1)
    do i = 2,n
        do j = 2,m
           tmp    = 5.0*dt*(2.0*f(i - 1,j - 1) + f(i - 1,j) &
                    + f(i, j - 1) )/4.0
           f(i,j) = (obj(x(i),y(j))*dt - tmp &
                    - x(i)*f(i - 1,j) -y(j)*f(i,j - 1))/&
                    (x(i) + y(j)) 
           !f(i,j) = (obj(x(i),y(j))*dt - tmp &
           !         - x(i)*f(i - 1,j) -y(j)*f(i,j - 1))/&
           !         (x(i) + y(j)) 
           !f(j,i) = (obj(x(j),y(i))*dt - x(j)*f(j - 1,i)&
           !         -y(i)*f(j,i - 1))/(x(j) + y(i)) 
        end do
    end do
    call cpu_time(t2)
    print *,"time cost is ==> ",t2 - t1

    call cpu_time(t1)
    do i = 1,nfre
        do j = 1,mfre
            itmp  =  nnum*i
            jtmp  =  mnum*j
            write(10,"(3f12.6)")x(itmp),y(jtmp),f(itmp,jtmp)
        end do
        if(i == nfre/2)then
            write(10,*)"     "
        endif
    end do
    call cpu_time(t2)
    print *,"time cost is ==> ",t2 - t1

    close(10)

end program main
