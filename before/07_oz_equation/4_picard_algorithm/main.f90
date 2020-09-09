program main
    use module_fst
    implicit none

    filename = "result.txt"
    open(20, file = filename)
    ! initial dr and dk
    do i = 1,n
        dr(i) = (i - 1)*deltar
        dk(i) = (i - 1)*deltak
    end do
    
    !*******check the result *************

    call check()

    !************************************

    !  print crmm and test to make a comparision
    write(20,*)"  i   ","    new_crmm    ","    old_crmm   "
    do i = 2,n
        write(20,"(3f15.7)")dr(i), test(i)/dr(i), crmm(i)/dr(i)
    end do
    close(10)
    close(20)
end program main
