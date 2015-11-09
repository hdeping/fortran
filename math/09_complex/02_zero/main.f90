program main
    use module_common
    implicit none
    real(8)               :: tmpnew
    
    filename = "data.txt"
    !open(10,file = filename,status = "old",iostat = ierror)
    open(10,file = filename)

    !a(1:10) = (/5.0,2.0,- 4.0,10.0,- 1.0,- 1.0,- 9.0,9.8,3.9,10.9/)
    !a(11:20) = (/5.0,2.4,- 4.0,10.0,- 0.0,- 1.0,- 9.0,9.8,3.9,10.9/)
    !a(21:30) = (/9.0,1.0,- 4.0,10.0,- 1.0,- 1.0,- 9.0,9.0,3.9,10.9/)
    !a(31:40) = (/0.0,1.3,- 4.0,00.0,- 0.0,- 1.0,- 9.0,9.8,3.9,10.9/)
    !a(41:51) = (/1.0,2.0,- 4.0,10.0,- 1.0,- 1.0,- 9.0,9.8,444.0,10.9,90.9/)

    call random_seed()

    !  get coefficients 
    do i = 1,m+1
        call random_number(a(i))
        a(i) = a(i)*500.0 - 250.0
    end do
    ! get the roots
    roots   = getall()
        

    print *,"cycle times is ",cyTimes
    print *,"time cost is ",t2 - t1
    !do i = 1,m+1
    !    write(10,"(f12.6)")a(i)
    !end do
    do i = 1,m
        write(10,"(2f12.6)")roots(i,:)
    end do
    
    !  check
    call check()
    
   



    close(10)
    
end program main
!test{{{
!    do i = - 10,10
!        do j = - 10,10
!           x = i*dt
!           y = j*dt
!           call getr_theta(r,theta,x,y)
!           print "('x = ',f18.9,'y = ',f18.9,'r = ',f18.9,'theta =&
!           ',f18.9)",x,y,r,theta
!           pause
!        end do
!    end do
    !x = 2.0
    !y = getPoly(x,a)
    !print *,"y = ",y

    !r     = 2.0
    !theta = pi/2.0
    !print *,getg(r,theta,a)

    !r     = 2.0
    !theta = pi/2.0
    !print *,geth(r,theta,a)

    !do i = - n,n
    !    do j = - n,n
    !       x   = i*dt
    !       y   = j*dt
    !       !call cpu_time(t1)
    !       write(10,"(3f18.6)")x,y,getNorm(x,y,a)
    !       !call cpu_time(t2)
    !       !print *,"time cost is ",t2 - t1
    !       !pause
    !    end do
    !end do

    !tmp = 3.0
    !do 
    !    tmpnew = - getPoly(tmp,a,m+1)/getPoly(tmp,c,m)
    !    judge  = abs(tmpnew)
    !    if(judge < error)exit
    !    tmp    = tmp + tmpnew
    !    print *,"judge = ",judge
    !    pause
    !end do
    !print *,"f(x) = ",getPoly(tmp,a,m+1)
    !print *,"x = ",tmp
    !r     = tmp
    !theta = 0.0
    !print *,"g = ",getg(r,theta,a)," h = ",geth(r,theta,a)
    
    !theta = 1000
    !print *,gettheta(theta)
    !pause

    !do i = - n,n
    !    do j = - n,n
    !        x = i*dt
    !        y = j*dt
    !        call getr_theta(r,theta,x,y)
    !        write(10,*)x,y,getroot(r,theta)
    !    end do
    !end do

    !filename = "data.txt"
    !open(10,file = filename,status = "old",iostat = ierror)
    !do 
    !    read(10,*,iostat = ierror)r,theta
    !    if(ierror /= 0)exit
    !    print *,"norm is => ",getNorm(r,theta,a)
    !    print *,"norm is => ",getNorm(r,theta,a)
    !end do

    !do i = 1,m
    !    c(i) = i*a(i+1)
    !end do



!}}}
