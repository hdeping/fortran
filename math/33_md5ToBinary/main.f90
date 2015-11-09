program main
    use module_common
    implicit none
    integer  ierror
    character(len = 1)     :: ch(m),ch_input(m)
    character(len = 4)     :: ch4(m)
    character(len = 128)   :: ch_out = ""


    ! initial the array
    ch4=(/"0000","0001","0010","0011",  &
          "0100","0101","0110","0111",  &
          "1000","1001","1010","1011",  &
          "1100","1101","1110","1111"/)
    ch=(/"0","1","2","3","4","5","6","7", &
         "8","9","a","b","c","d","e","f"/)
    
    filename = "data.txt"
    open(10, file = filename,status = "old",iostat = ierror)

    do 
        read(10,"(16a)",iostat = ierror)ch_input(:)
        ch_out = ""
        write(*,*)ch_input
        do i = 1,m
            do j = 1,m
                if ( ch_input(i) == ch(j) )then
                    ch_out = trim(ch_out)//ch4(j)
                    exit
                endif ! if ends
            enddo !cycle ends
        enddo !cycle ends
        write(*,*)ch_out
        if ( ierror /= 0 )then
            exit
        endif ! if ends
        
    enddo !cycle ends
    
    

    close(10)
    




end program main
