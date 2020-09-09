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
        ! write the data into files
        call cpu_time(t1)

        filename = "mct.txt"

        !open(10,file = filename,form = "binary")
        !do i = 1,ncut
        !    do j = 1,tmnum/2
        !        write(10)j*dt,f(1,1:2,i,j),f(2,2,i,j)
        !    end do
        !end do
        !close(10)
        do i = 1,tmnum/2
            write(filename,"('t',i3.3,'.txt')")i
            open(10,file = filename)
            do j = 1,ncut
                write(10,"(4D18.6)")dk(j),f(1,1:2,j,i),f(2,2,j,i)
                !if(j == 1)then
                !print *,j*dt,f(1,1:2,i,j),f(2,2,i,j)
                !endif
            end do
            close(10)
        end do


        call cpu_time(t2)
    !enddo
    print *,"time cost of writing is ==> ",t2 - t1


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

    !filename = "data.dat"
    !open(10,file = filename)
    !do i = 1,ncut
    !    tmp = (i*deltar)**2.0
    !    write(10,*)i*deltar,tmp
    !end do
    !close(10)



!}}}
