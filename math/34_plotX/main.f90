program main
    use module_common
    implicit none

    real(8)   x,y

    
    filename = "data.txt"
    open(10, file = filename)
    

    x = - 2.99
    do i = 1,600
        x = x + 0.01
        if ( x < 0.0 )then
            y = - x*(x - 2.0)
        else
            y = x*(x - 2.0)
        endif ! if ends
        write(10,*)x,y
        
    enddo !cycle ends
     

    close(10)
end program main
