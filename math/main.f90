program main
    use module_common
    implicit none
    real(8)        :: a(n,n)
    real(8)        :: c(n,n)


    call random_seed()
    a(1,:) = (/1.0,3.0,4.0/)
    a(2,2:3) = (/2.0,5.0/)
    a(3,3) = 1.0
    do i = 1,n
        do j = 1,n - i 
            a(i,j) = a(j,i)
        enddo !cycle ends
    enddo !cycle ends
    
    c = mat_sqare_root(a,n)

end program main
