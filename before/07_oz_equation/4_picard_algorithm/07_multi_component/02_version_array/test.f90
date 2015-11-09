program test
    !use module_fst 
    use module_algebra
    use module_fst
    use module_evolution
    implicit none
    integer,parameter         :: nnew = 3
    real(8)                   :: a_array(nnew,nnew)
    real(8)                   :: b_array(nnew)
    real(8)                   :: res(nnew,nnew)
    real(8)                   :: x_res(nnew)
    real(8)                   :: b_new_array(nnew)
    real(8)                   :: mayer_r(n)
    real(8)                   :: mayer_k(n)
    real(8)                   :: mayer_O(n)
    integer                   :: ii
    call random_seed()

    ! test the iteration method to
    ! solve the linear equations

    a_array(1,:) = (/4.0,- 1.0,2.0/)
    a_array(2,:) = (/- 1.0,- 5.0,1.0/)
    a_array(3,:) = (/2.0,- 1.0,6.0/)
    b_array(:)   = (/1.0,2.0,3.0/)

    x_res = iter_sol_equ(a_array,b_array,nnew)
    do i = 1,nnew
        print *,i,x_res(i)
    end do
    filename = "mayer_O.txt"
    open(10,file = filename)
    
    do ii = 1,n
        dr(ii) = (ii - 1.0)*deltar
        dk(ii) = (ii - 1.0)*deltak
    end do
    mayer_r = may(d12)
    do ii = 1,n
        mayer_r(ii) = mayer_r(ii)*dr(ii)
    end do
    mayer_k = fst(mayer_r,1)
    do ii  = 1,n
        mayer_k(ii) = mayer_k(ii)**2.0
    end do
    mayer_O = fst(mayer_k,- 1)
    do ii = 2,n
        write(10,"(2f15.7)")dr(ii),mayer_O(ii)/dr(ii)
    end do


    close(10)
    





    !  test module_algebra
    !do i = 1,nnew
    !    do j = 1,nnew
    !        call random_number(x1)
    !        a_array(i,j) = dble(i + j)*x1
    !        b_array(i,j) = dble(i+j)
    !    end do
    !end do
    !call cpu_time(t1)
    !res = sol_mat(a_array,b_array,nnew)
    !call cpu_time(t2)
    !print *,t2 - t1
    !b_new_array = multi_mat(a_array,res,nnew,nnew,nnew)
    !do i = 1,nnew
    !    print "(i4,10f12.5)",i,b_new_array(i,1:10)
    !end do
    

    

end program test
