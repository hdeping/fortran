program main
    !use module_fst 
    use module_common
    use module_evolution
    integer               :: itmp
    call random_seed()

!get the solution of OZ{{{
    filename = "gr.txt"
    open(10,file = filename)

    !  initial the k and r
    do i = 1,n
        dk(i) = (i - 1)*deltak
        dr(i) = (i - 1)*deltar
    end do
    !  calculate the mayer function
    may11 = may(d11)
    may12 = may(d12)
    may22 = may(d22)

    call cpu_time(t1)
!*********** solve the equation **********
    call evolution()
!*****************************************
    call cpu_time(t2)
    print *,"time cost is ",t2 - t1

    !  write gr to the file
            do itmp = 2,n
                gr11(itmp) = (cr11(itmp) + gamma_r11(itmp))/dr(itmp) + 1.0
                gr12(itmp) = (cr12(itmp) + gamma_r12(itmp))/dr(itmp) + 1.0
                gr22(itmp) = (cr22(itmp) + gamma_r22(itmp))/dr(itmp) + 1.0
            end do
    do ii = 2,n
        write(10,"(4f12.5)")dr(ii),gr11(ii),gr12(ii),gr22(ii)
    end do
    close(10)
!}}}
end program main
