!  to make a fourier transform
!  y = F*x
program main
    implicit none
    integer,parameter          :: l  = 6
    integer,parameter          :: n  = int(1E3) 
    real,parameter             :: pi = 3.14159265 
    real,parameter             :: dx = 1E-2   ! the unit of x
    real                       :: f_re(n,n)   ! the real part of the matrix
    real                       :: f_im(n,n)   ! the imagine part of the matrix
    real                       :: x_re(n)        ! the input array
    real                       :: x_im(n)        ! the input array
    real                       :: y_re(n)        ! the output array
    real                       :: y_im(n)        ! the output array
    real                       :: tmp
    real                       :: t1
    real                       :: t2
    integer                    :: i 
    integer                    :: j 
    integer                    :: si             ! si = i -1             
    integer                    :: sj             ! sj = j -1
    integer                    :: k 
    character(10)              :: filename

    filename = "data.txt"
    open(10, file = filename)
    ! initial the input matrix
    do i = 1,n
        x_re(i) =  func(i)
    end do
    ! calculate the matrix F 
    call cpu_time(t1)
    do i = 1,n
        do j = i,n
            si = i - 1
            sj = j - 1
            f_re(i,j) = cos(2*pi*si*sj/n)
            f_im(i,j) = sin(2*pi*si*sj/n)
        end do
    end do
    do i = 2,n
        do j = 1, i - 1
            f_re(i,j) = f_re(j,i)
        end do
    end do
    call cpu_time(t2)
    print *,"time cost is ",t2 - t1
    ! calculate the forward DFT
    ! real part
    do i = 1,n
        tmp = 0 
        do j = 1,n 
            tmp = tmp + f_re(i,j)*x_re(j) + f_im(i,j)*x_im(j)
        end do
        y_re(i) = tmp
    end do
    ! imaginary part
    do i = 1,n
        tmp = 0 
        do j = 1,n 
            tmp = tmp + f_re(i,j)*x_im(j) + f_im(i,j)*x_re(j)
        end do
        y_im(i) = tmp
    end do
    ! print the result
    write(10,*)"the forward DFT of x"
    do i = 1,n
        write(10,"(3f12.5)")i*dx, y_re(i), y_im(i)
    end do
!****************************************************
    ! calculate the backward DFT
    ! real part
    do i = 1,n
        tmp = 0 
        do j = 1,n 
            tmp = tmp + f_re(i,j)*x_re(j) - f_im(i,j)*x_im(j)
        end do
        y_re(i) = tmp/n
    end do
    ! imaginary part
    do i = 1,n
        tmp = 0 
        do j = 1,n 
            tmp = tmp + f_re(i,j)*x_im(j) - f_im(i,j)*x_re(j)
        end do
        y_im(i) = tmp/n
    end do
    ! print the result
    write(10,*)"the forward DFT of x"
    do i = 1,n
        write(10,"(3f12.5)")i*dx, x_re(i), x_im(i),func(i)
    end do
    stop

    close(10)


    contains
!function func{{{
      real function func(i)
          integer,intent(in)    :: i
          real                  :: x
          x = i*dx
          func = sqrt(x)
      end function func
!}}}
end program main
