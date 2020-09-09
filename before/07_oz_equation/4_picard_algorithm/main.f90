program main
    !use module_fst 
    use module_common
    use module_fst
    integer  itmp

    real(8)    :: a(n)
    real(8)    :: b(n)
    real(8)    :: c(n)
    

    filename = "gr.txt"
    open(10,file = filename)
    filename = "cr.txt"
    open(20,file = filename)
    filename = "sk.txt"
    open(30,file = filename)
    !  initial the k and r
    do i = 1,n
        dk(i) = (i - 1)*deltak
        dr(i) = (i - 1)*deltar
    end do

    maymm = may(dmm)
!********************************************************** 

    !  solve the equation
    call cpu_time(t1)
    call evolution()
    call cpu_time(t2)

    print *,"lambda is ",lambda
    !  write ck to the file
    g_rmm(1) = 0.0

    do i = 2,n 
       g_rmm(i) = (crmm(i)  +  grmm(i) )/dr(i)  + 1
    end do   !  i
    ! calculate crmm
    crmm(1) = 0.0
    do i = 2,n
        crmm(i) = crmm(i)/dr(i)
    end do
    ! print the results
    do i = 1,n
        write(10,*)dr(i),g_rmm(i) 
        write(20,*)dr(i),crmm(i) 
    end do
    print *,"time cost is ",t2 - t1

    close(10)
    !print *,hkmm(10)
end program main
