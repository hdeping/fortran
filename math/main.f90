program main
    use module_common
    !it is a problem that there is a sequence {1..9},
    !how to insert "+","-" or "" into any position between
    !the numbers such as 1 - 2 + 345678 + 9, and make the 
    !final answer is 100? Try to write a program to 
    !solve it

    integer times,kk
    integer symbol
    integer asum,bsum  ! used in the calculation
    character(len = 1)   :: ch(3) = (/' ','+','-'/)
    character(len = 1)   :: ch1(n) 

    do i = 1,n
        a(i) = i
    enddo !cycle ends
     
    b(1) = 1
    do i = 0,total
        b(2:n) = getTernary(i)
        asum = 0
        bsum = 0
        do j = 1,n
            if ( b(j) == 1 )then
                symbol = 1
                asum = asum + bsum
                bsum = symbol*a(j)
            elseif ( b(j) == 2 )then
                symbol = - 1
                asum = asum + bsum
                bsum = symbol*a(j)
            elseif ( b(j) == 0 )then
                bsum = 10 * bsum + symbol*a(j)
            endif ! if ends
        enddo !cycle ends
        asum = asum + bsum
        if ( asum == 100 )then
            print 100,b
        endif ! if ends
    enddo !cycle ends

     

100 format (10i1)

end program main
