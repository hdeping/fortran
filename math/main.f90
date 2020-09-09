program main
    use module_common
    implicit none

    integer               :: n

    ! initial the status

    do n = 10,50
         allocate(a(n))
         allocate(b(n))
         a = 0
         c = getFinal(n,m)
         print *,n,c
         deallocate(a)
         deallocate(b)
    enddo !cycle ends
    ! game begins
     

end program main
