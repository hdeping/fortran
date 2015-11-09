program main
    use module_common
    use module_fst
    real(8)     a(n)
    real(8)     b(n)
    real(8)     c(n)

    do i = 1,n
        a(i) = i
    end do
    call fst(b,a,n,l,1)
    ! print 
    call fst(c,b,n,l,-1)
    !write(*,*)"a     b      c"
    !do i = 1,n
    !    write(*,"(3f10.5)")a(i),b(i),c(i)
    !end do
end program main
