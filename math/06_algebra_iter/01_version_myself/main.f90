program main
    use module_algebra
    implicit none
    integer,parameter      :: n = 4
    real(8)                :: a(n,n)
    real(8)                :: b(n)
    real(8)                :: x(n)
    real(8)                :: s(n)   ! abs(s(i)) < 1
    real(8)                :: t1
    real(8)                :: t2
    real(8)                :: x1
    integer                :: i
    integer                :: j
    integer                :: ii
    integer                :: jj

    call random_seed()
    do ii = 1,n
        do jj = 1,n
            call random_number(x1)
            a(ii,jj) = x1*10000.0 - 5000.0
            if(ii == jj)a(ii,jj) = 0.0
        end do
        call random_number(x1)
        b(ii) = 10000*x1 - 5000.0
        call random_number(s(ii))
    end do
    a = sol_mat(a,multi_mat(a,b,n,n,n),n)
!program before{{{
   ! a(1,:) = (/0.0,0.0,2.0,3.0/)
   ! a(2,:) = (/5.0,4.0,0.0,0.0/)
   ! a(3,:) = (/3.0,0.0,3.0,0.0/)
   ! a(4,:) = (/0.0,3.0,0.0,9.0/)
!}}}
    call cpu_time(t1)
    x = iter_sol_equ(a,b,n) 
    !x = iter_sol_new(a,b,n) 
    !x1 = det(a,n)
    print *,"det = ",x1
    call cpu_time(t2)
    print *,"time cost is ==>  ",t2 - t1
end program main
