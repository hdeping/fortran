program main
    use module_common
    implicit none
    real(8)                :: a(n,n)
    real(8)                :: b(n,n)
    real(8)                :: x1
    integer                :: m


    call random_seed()
    do i = 1,n
        do j = 1,n
            call random_number(x1)
            a(i,j) = 10*x1
        end do
    end do
    m = 10
    b = power_mat(a,n,m)

end program main
