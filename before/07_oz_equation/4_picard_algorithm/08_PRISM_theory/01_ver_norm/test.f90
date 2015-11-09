program test
    !use module_fst 
    use module_common
    use module_evolution
    integer,parameter         :: nnew = 1100
    real(8)                   :: a_array(nnew,nnew)
    real(8)                   :: b_array(nnew)
    real(8)                   :: res(nnew)
    real(8)                   :: b_new_array(nnew)
    call random_seed()
    !  test module_algebra
    do i = 1,nnew
        do j = 1,nnew
            call random_number(x1)
            a_array(i,j) = dble(i + j)*x1
        end do
        b_array(i) = dble(i)
    end do
    call cpu_time(t1)
    res = sol_equ(a_array,b_array,nnew)
    call cpu_time(t2)
    print *,t2 - t1
    pause
    b_new_array = multi_vec(a_array,res,nnew,nnew)
    do i = 1,nnew
        print *,i,b_new_array(i)
    end do
    

    

end program test
