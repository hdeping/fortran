module module_new
   use ogpf
   use module_common
  
   
   
    contains
!subroutine data_line{{{
   subroutine data_process(x,y,output,n)
       ! draw the data line with * representing the point
        integer,intent(in)       :: n
        real(8),intent(in)       :: x(n),y(n)
        character(10),intent(in) :: output
        type(gpf)                :: gp
        character(50)            :: setout

        ! sample i. plot on the screen
        !annotation, set title, xlabel, ylabel
        call gp%title('example 1. a simple xy plot')
        call gp%xlabel('my x axis x')
        call gp%ylabel('my y axis y')
        call gp%options('set style data linespoints')
        setout = 'set output "'//trim(output)//'.pdf"'
        call gp%options(setout)

        !call plot to draw a vector against a vector of data
        call gp%plot(x,y)
    end subroutine data_process
!}}}
!{{{
function besj0(x)
    real(8),intent(in)  :: x
    real(8)             :: tau
    real(8)             :: besj0
    integer             :: ii
    real                :: t1
    real                :: t2

    !call cpu_time(t1)
    besj0 = 0.0
    do ii = 1,total
        tau   = dble(ii)*dx
        besj0 = besj0 + dx/pi*cos(x*sin(tau))
    end do
    !call cpu_time(t2)
    !print *,"time cost is ==>",t2 - t1

    
end function besj0
!}}}
end module module_new
