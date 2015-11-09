program main
    use module_function

    implicit none

    filename = "solution.txt"
    open(10,file = filename)

    call random_seed()
    !do i = 1,100000
    !    call random_number(t1) 
    !    call random_number(t2) 
    !    call random_number(t3) 
    !    call random_number(t4) 
    !    t1 = t1*10000 - 5000
    !    t2 = t2*10000 - 5000
    !    t3 = t3*10000 - 5000
    !    t4 = t4*10000 - 5000
    !    x  = (/t1,t2,t3,t4/)
    !    times = 0
    !    call iter_NR()
    !    if(times /= fre)write(10,"(<n>f12.6)")x(:)
    !    if(mod(i,fre) == 0)print *,x
    !end do
    x = (/1.974018,0.160667,2.473845,2.159975/)
    print *,fun1(x),fun2(x),fun3(x),fun4(x)

    close(10)

end program main
