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
    do i = 1,m
        do j = i,m
            mayfun(i,j,:) = may(d(i + j -1))
        end do
    end do

    call cpu_time(t1)
!*********** solve the equation **********
    call evolution()
!*****************************************
    call cpu_time(t2)
    print *,"time cost is ",t2 - t1

    !  write gr to the file
    do i = 1,m
        do j = i,m
            do itmp = 2,n
                gr(i,j,itmp) = (cr(i,j,itmp) + gamma_r(i,j,itmp))/dr(itmp) + 1.0
            end do
        end do
    end do
    do ii = 1,n
        write(10,"(4f12.5)")dr(ii),gr(1,1:m,ii),gr(2,m,ii)
    end do
    close(10)
!}}}
end program main
