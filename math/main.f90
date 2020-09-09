program main
    use module_common
    implicit none
    integer         :: times = 0

    filename = "data.txt"
    open(10,file = filename)

    write(10,*)"times"

    do i = 1,n - 4
        do j = i + 1,n - 3
            do k = j + 1,n - 2
                do kk = k+1,n - 1
                    do jj = kk + 1,n
                        if (  i+j+k+kk+jj < 100)then
                            cycle
                        endif ! if ends
                        if ( tail(i,j,k,kk,jj) /= 0 )then
                            times = times + 1
                            write(10,"(8i6)")times,i,j,k,kk,jj
                        endif ! if ends
                    enddo !cycle ends
                enddo !cycle ends
            enddo !cycle ends
        enddo !cycle ends
    enddo !cycle ends
    
     
end program main
