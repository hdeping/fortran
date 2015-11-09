program main
    use blas95
    implicit none
    integer,parameter       :: m = 2
    integer,parameter       :: n = 3
    real                    :: alpha = 3
    real                    :: beta  = 2
    real                    :: a(m,n)
    real                    :: x(n)
    real                    :: y(m)
    real                    :: res(m)
    real                    :: tmp 
    integer                 :: kl    ! the number of the super-diagonals    
    integer                 :: ku    ! the number of the sub-diagonals
    integer                 :: i
    integer                 :: j
    integer                 :: k
    integer                 :: u
    character               :: trans
     
    a(1,:) = (/2,3,1/)
    a(2,:) = (/3,1,4/)
    x(:)   = (/2,3,4/)
    y(:)   = (/1,3/)

    tmp = asum(y)
    print *,tmp

    do i = 1,m
       tmp = 0 
       do j = 1,n
           tmp = tmp + a(i,j)*x(j)
       end do
       res(i) = alpha*tmp + beta*y(i)
    end do
    print *,"the normal one"
    do i = 1,m
        print *,res(i)
    end do
   

    trans = 'n'
    kl = 1
!    call gbmv(a,x,y,kl,m,,alpha,beta,trans)
    print *,"matrix a"
    do i = 1,m
        print *,a(i,:)
    end do
    print *,"array x"
    print *,x(:)
    print *,"array y"
    print *,y(:)






end program main
