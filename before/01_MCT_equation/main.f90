program main
    !use module_fst 
    use module_mct
    integer               :: itmp
    integer               :: jtmp
    call getmatrixes()
! get MCT
dt = 5.0E-6
!do i = 1,48
    dt = dt*2
    call cpu_time(t1)
    ! get MCT equations
    lambda = getmct() 
    call cpu_time(t2)
    print *,"time cost is ==> ",t2 - t1
    do i = 1,n
        write(filename,"('q',i3.3,'.txt')")i
        print *,"wirte ",i
        open(10,file = filename)
        do j = 1,tmnum/2
            write(10,"(4f15.6)")j*dt,f(1,1:2,i,j),f(2,2,i,j)
            !print *,j
        end do
        close(10)
    end do


end program main
!code before{{{
    !do q = 1,300
    !        
    !    wri
    !    open(10,file = filename)
    !    do i = 1,n
    !        write(10,"(4f12.5)")(i - 1)*h,f(1,1,i,3),f(1,2,i,3),f(2,2,i,3)
    !    end do
    !    close(10)
    !end do


!}}}
