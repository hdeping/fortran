program main
    use module_fst 
    use module_common
    !use fst_module
    integer  itmp
    real(8)          a(n)
    real(8)          b(n)
    real(8)          c(n)
    
    do itmp = 1,n
        !a(itmp) = sin(deltar*dble(itmp))
       a(itmp) = 0.5
    end do
    b = fst(a, 1)
    c = fst(b,- 1)
    do i = 1,n
        write(*,"(i5,3f12.5)")i,a(i),b(i),c(i)
    end do
    print *,"n = ",n


end program main
