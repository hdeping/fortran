program main
    !use module_fst 
    use module_common
    use module_evolution
    call random_seed()

!get the solution of OZ{{{
    filename = "gr.txt"
    open(20,file = filename)
    filename = "sk.txt"
    open(30,file = filename)

    !  initial the k and r
    do i = 1,n
        dk(i) = (i - 1)*deltak
        dr(i) = (i - 1)*deltar
    end do
    !  calculate the mayer function
    maypp = may(dpp)
    maywp = may(dwp)
    mayww = may(dww)
    maynp = may(dnp)
    maynw = may(dnw)

    call cpu_time(t1)
!*********** solve the equation **********
    call evolution()
!*********** nano-particle **********
    call evolution_nano()
!*****************************************
    !call check()
    !lambda = judge(test,ckpp)
    !lambda1 = judge(test1,ckwp)
    !lambda2 = judge(test2,ckww)
    !print *,lambda,lambda1,lambda2

    call cpu_time(t2)

    print *,"lambda is ",lambda
    !  write ck to the file
     grpp(1) = 0.0 
     grwp(1) = 0.0 
     grpp(1) = 0.0 
    do i = 2,n
       grpp(i) = (crpp(i)  +  gamma_rpp(i) )/dr(i)  + 1
       grwp(i) = (crwp(i)  +  gamma_rwp(i) )/dr(i)  + 1
       grww(i) = (crww(i)  +  gamma_rww(i) )/dr(i)  + 1
       grnp(i) = (crnp(i)  +  gamma_rnp(i) )/dr(i)  + 1
       grnw(i) = (crnw(i)  +  gamma_rnw(i) )/dr(i)  + 1
    end do   !  i
    do i = 1,n
        write(20,"(6f10.5)")dr(i),grpp(i),grwp(i),grww(i),grnp(i),grnw(i)
    end do

    ! get sk

    close(30)
    close(20)
!}}}
end program main
