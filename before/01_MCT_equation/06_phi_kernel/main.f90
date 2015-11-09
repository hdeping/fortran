program main
    !use module_fst 
    use module_mct
    integer               :: itmp
    integer               :: jtmp
    integer               :: sunit(ncut)

    !  prepare some variables and parameters
    call getmatrixes()
    !  open files
    do i  = 1,ncut
        sunit(i) = 100 + i
        write(filename,"('q',i3.3,'.txt')")i
        open(sunit(i),file = filename)
    end do
    ! get MCT
    dt = 1.0E-5

    call cpu_time(t1)
    call getmct()
    call cpu_time(t2)
    print *,"time cost is ==> ",t2 - t1
    ! print first tmnum data
    do itmp = 1,ncut
        do jtmp = 1,tmnum
           write(sunit(itmp),"(2f18.6)")dt*jtmp,f(itmp,jtmp)
        end do
    end do
    !  get the data of the rest time
    do i = 1,14
        call cpu_time(t1)
        dt = dt*2.0
        ! get a solution cycle
        call getmct2()
        do itmp = 1,ncut
            do jtmp = tmnum/2 + 1,tmnum
               write(sunit(itmp),"(2f18.6)")dt*jtmp,f(itmp,jtmp)
            end do
        end do
        call cpu_time(t2)
        print *,"time cost is ==> ",t2 - t1
    enddo
     
     do i  = 1,ncut
         close(sunit(i))
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

    !filename = "data.dat"
    !open(10,file = filename)
    !do i = 1,ncut
    !    tmp = (i*deltar)**2.0
    !    write(10,*)i*deltar,tmp
    !end do
    !close(10)

    !open(10,file = filename,form = "binary")
    !do i = 1,ncut
    !    do j = 1,tmnum/2
    !        write(10)j*dt,f(1,1:2,i,j),f(2,2,i,j)
    !    end do
    !end do
    !close(10)
!}}}
