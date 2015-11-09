program main
    use module_common
    use module_fft
    real(8)     a(n)
    real(8)     b(n)
    real(8)     c(n)
    real        t1
    real        t2

    call cpu_time(t1)
    do i = 1,n
        a(i) = i
    end do
    call fft(b,a,n,l,1)
    ! print 
    call fft(c,b,n,l,-1)
    call cpu_time(t1)

    print *,"time cost is ",t2 - t1

    !write(*,*)"a     b      c"
    !do i = 1,n
    !    write(*,"(3f10.5)")a(i),b(i),c(i)
    !end do
end program main
