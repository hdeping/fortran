program test
    use module_common
    implicit none
    real(8)           :: x1 = 10000
    angle = x1
    angle = ch_range(angle,num,2.0*pi)
    do i = 1,num
        print *,i,angle(i)
    enddo !cycle ends
     
    
end program test 
