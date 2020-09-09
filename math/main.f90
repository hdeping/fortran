program main
    use module_common
    implicit none
    integer                :: times
    
    filename = "data.txt"
    open(10,file = filename)

    a = (/5.0,2.0,- 4.0,10.0,- 1.0,- 1.0/)
   

    times = 0
    do i = - n,n
        do j = - n,n
           !times = times + 1
           !if(mod(times,40) /= 0)cycle
           x     = i*dt
           y     = j*dt
           !call cpu_time(t1)
           write(10,"(3f18.6)")x,y,getNorm(x,y,a)
           !call cpu_time(t2)
           !print *,"time cost is ",t2 - t1
           !pause
        end do
    end do


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
!}}}
