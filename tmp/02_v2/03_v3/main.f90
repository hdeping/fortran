program main
    implicit none
    integer,parameter         :: n = 10
    integer   i,j,k

    do i = 1,n
        do j = 1,n 
            do k = 1,n
                print *,i+j+k
            enddo !cycle ends
             
        enddo !cycle ends
         
    enddo !cycle ends
     
end program main
