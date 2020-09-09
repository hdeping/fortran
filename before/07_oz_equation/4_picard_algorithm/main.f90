program main
    use module_common
    use module_fst

    filename = "result.txt"
    open(20, file = filename)
    ! initial dr and dk
    do i = 1,n
        dr(i) = (i - 1)*deltar
        dk(i) = (i - 1)*deltak
    end do
    ! calculate mayer function of m-m
    maymm = may(dmm)
    
    !*********  check the results ************

    call check()

    !******************************************

    !  print crmm and test to make a comparision
    write(20,*)"      r       ","    new_grmm    ","    old_grmm   "
    do i = 2,n
        write(20,"(3f15.7)")dr(i), test(i)/dr(i), grmm(i)/dr(i)
    end do
    close(20)
end program main
