! In general, a conical section
! can be given by five points(for 2D)
! In a n dimmesional space, the number 
! of the points are n*(n+3)/2
program main
    use module_common
    use module_algebra
    implicit none
    
    !  read data
    filename = "data.txt"
    open(10,file = filename,status = "old",iostat = ierror)
    do i = 1,m
        read(10,*,iostat = ierror)x(i,:)
        if(ierror /= 0)exit
    end do

    !  get a and b
    b = 1.0

    do i = 1,m
        num = 0
        do j = 1,d
            do k = j,d
                num = num + 1
                a(i,num) = x(i,j)*x(i,k)
            end do
        end do
        num = num + 1
        do j = num, m
            a(i,j) = x(i,j - num)
        end do
    end do
    y = sol_equ(a,b,m)
    ! print the coefficients
    print *,"y = "
    do i = 1,m
        print *,y(i)
    end do
    print *,"x = "
    do i = 1,m
        print *,x(i,:)
    end do
    close(10)
end program main
