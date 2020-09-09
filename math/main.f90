program main
    use module_common
    implicit none
    
    do i = 1,m
        x(i) = 1D-2*dble(i)
    enddo !cycle ends
    y = laguerre(x) 
    do i = 1,m
        write(filename,"('laguerre',i3.3,'.txt')")i
        !write(filename,"('legendre',i3.3,'.txt')")i
        open(10,file = filename)
        do j = 1,n
            write(10,*)j,y(j,i)
        enddo !cycle ends
        close(10) 
    enddo !cycle ends
     
end program main
