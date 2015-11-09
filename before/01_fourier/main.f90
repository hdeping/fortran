program main
    use blas95
    use mkl_dfti
    use lapack95
    implicit none
    integer,parameter   :: n = int(1E6)
    integer             :: i
    integer             :: j
    integer             :: k
    real(8)             :: x1
    real(8)             :: y1
    real(8)             :: x2
    real(8)             :: y2
    real(8)             :: x
    real(8)             :: y
    real(8)             :: a(n)

    do i = 1,10
        print "(i3)",i
    end do
    
    do i = 1,n
        a(i) = i
    end do
    x1 = asum(a)
    print *,x1
end program main
