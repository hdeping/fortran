program main
    use module_common
    use module_fst
    integer       ii

    do ii = 1,n
        a(ii) = ii
    end do

    test  = a
    do ii = 1,1000
        b = fst(test,1)
        c = fst(b,- 1)
        test = c
    end do
    do ii =  1,n
        print *,a(ii),b(ii),c(ii)
    end do
    lambda = judge(a,c)
    print *,"lambda is ==>",lambda
    
    

end program main
