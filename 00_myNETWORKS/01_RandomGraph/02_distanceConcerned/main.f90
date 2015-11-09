program main
    use module_common
    implicit none

    real                :: degree
    real                :: beta
    integer             :: ii

    filename = "data.txt"
    open(10, file = filename)
    

    !$omp parallel
       !$omp do 
       
    ii = 1
    do  
        beta = ii*1.0
        degree = getMeanDegree(beta)
        write(10,*)beta,degree
        if ( ii < 1000 )then
            ii = ii + 10
        else
            ii = ii + 200
        endif ! if ends
        
        if ( ii > 10000 )then
            exit
        endif ! if ends
        
    enddo !cycle ends
       !$omp end do
    !$omp end parallel
    

    close(10)

end program main
