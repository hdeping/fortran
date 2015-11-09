program main
    use module_common
    implicit none
    ! record the status of the link
    ! 1 for connect , 0 for disconnect
    integer    :: a(n,n)
    integer    :: b(n),c(m)
    real       :: x1,x2,x3



    a = 0
    call random_seed() 
    ! random network evolution
    do k = 1,total
        call random_number( x1 )
        i = int(n*x1) + 1
        call random_number( x1 )
        j = int(n*x1) + 1
        call random_number (x1)
        if ( x1 < 0.5 )then
            a(i,j) = 1 - a(i,j) 
            a(j,i) = a(i,j)
        endif ! if ends

    enddo !cycle ends

    ! record the number
    do i = 1,n
        b(i) = sum(a(i,:))
    enddo !cycle ends
 

    ! print the b(n) 
    filename = "data.txt"
    open(10, file = filename)
    do i = 1,n
        write(10,*)b(i)
    enddo !cycle ends

    ! get the frequency  c(m)

    c = 0
    do i  = 1,n
       j = b(i)/10 + 1
       if ( j == m + 1 )then
           c(m) = c(m) + 1
       else
           c(j) = c(j) + 1
       endif ! if ends
    enddo !cycle ends
 

    ! print the c(m)
    filename = "file.txt"
    open(100, file = filename)
    do i = 1,m
        write(100,"(2i6)")i,c(i)
    enddo !cycle ends

    close(10)


    ! get the sigma of c(m)

    write(*,*)"sigma = ",getSigma(c,m)

    ! print the average link number

    write(*,*)"average number is ",sum(b(:))/ dble(n)

end program main

integer,allocatable :: seed(:)
integer             :: l
! spicial for gfortran
call random_seed( size=l )
allocate(seed(l))
call system_clock(count=clock)
do i = 1 , l 
    seed(i)=clock+37*(i-1)   
enddo
call random_seed(put=seed)   
deallocate(seed) 
!..........................................  
call random_number(s)  
