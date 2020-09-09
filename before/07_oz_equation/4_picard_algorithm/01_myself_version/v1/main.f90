program main
    use module_fst
    integer,parameter      ::  l = 5
    integer,parameter      ::  n = 2**l
    integer                ::  itmp 

    real(8)     a(n)
    real(8)     b(n)
    real(8)     c(n)

    call cpu_time(t1)
    do itmp = 1,n
        a(itmp) = itmp - 1 
    end do
    !   subroutine
    !    call fst(a,b,  1)
    !    call fst(b,c,- 1)
    !   function 
        b = fst(a,n,l,  1)
        c = fst(b,n,l,- 1)
    do itmp = 1,n
        print *,a(itmp),b(itmp),c(itmp)
    end do
    call cpu_time(t2)
    print *,"the time is ==>",t2 - t1

end program main
