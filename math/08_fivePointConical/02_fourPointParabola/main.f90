! In general, a conical section
! can be given by five points(for 2D)
! In a n dimmesional space, the number 
! of the points are n*(n+3)/2
program main
    use module_common
    use module_algebra
    implicit none
    real(8)             :: arraynew(m - 1)
    real(8)             :: arraytest(m - 1,m)
    real(8)             :: judge_delta
    
    !  read data
    filename = "data.txt"
    open(10,file = filename,status = "old",iostat = ierror)
    do i = 1,m - 1
        read(10,*,iostat = ierror)x(i,:)
        if(ierror /= 0)exit
    end do
    close(10)

    !  get a and b
    b = 1.0

    do i = 1,m - 1
        num = 0
        do j = 1,d
            do k = j,d
                num = num + 1
                a(i,num) = x(i,j)*x(i,k)
            end do
        end do
        do j = num + 1, m
            a(i,j) = x(i,j - num)
        end do
    end do
    ! save the current array a
    arraytest = a
    !  a,c,d,e be expressed by a 
    ! so, we should make transform 
    ! 1,2,3,4,5 to 1,3,4,5,2
    do i = 2,4
        arraynew(:) = a(:,i)
        a(:,i)      = a(:,i+1)
        a(:,i+1)    = arraynew(:)
    end do
    ! test the transform
    do i = 1,m - 1
        print "(<m>f18.6)",a(i,:)
    end do
    print *,"original array"
    do i = 1,m - 1
        print "(<m>f18.6)",arraytest(i,:)
    end do
    x_tmp = sol_equ2(a,b,m - 1)
    !  print x_tmp
    do i = 1,m - 1
        print "(<m>f18.9)",x_tmp(i,:)
    end do
    pause
    ! according to b**2 - 4ac = 0
    coef(1) = 1.0 - 4.0*x_tmp(1,2)*x_tmp(2,2)
    coef(2) = - 4.0*(x_tmp(1,1)*x_tmp(2,2) + &
              x_tmp(1,2)*x_tmp(2,1))
    coef(3) = - 4.0*x_tmp(2,1)*x_tmp(1,1)
    judge_delta = coef(2)**2.0 - 4.0*coef(1)*coef(3)
    print *,"coef = "
    print "(3D15.5)",coef(:)
    print *,"judge_delta = ",judge_delta
    pause
    if(judge_delta > 0.0)then
        bvalue(1) = (- coef(2) + sqrt(judge_delta))/2.0/coef(1)
        bvalue(2) = (- coef(2) - sqrt(judge_delta))/2.0/coef(1)
    endif
    ! output the coefficients with two occassions
    do i = 1,2
        do j = 1,m
            if(j == 1)then
                y(j) = x_tmp(j,1) + x_tmp(j,2)*bvalue(i)
            elseif(j == 2)then
                y(j) = bvalue(i)
            else
                y(j) = x_tmp(j - 1,1) + x_tmp(j - 1,2)*bvalue(i)
            endif
        end do
        ! output a,b,c,d,e
        print "(<m>D18.9)",y(:)
    end do

    ! print the coefficients
end program main
