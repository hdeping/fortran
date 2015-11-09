!interpretation{{{
    !-------------------------------------------------------------------------------
    !    gnuplot interface
    !-------------------------------------------------------------------------------
    !    purpose:   object based interface to gnuplot from fortran (ogpf)
    !    platform:  windows xp/vista/7
    !               (it should work on other platforms, see the write2gnuplot subroutine below)
    !    language:  fortran 2003 and 2008
    !    requires:  1. fortran 2003 compiler (e.g gfortran 4.7, ivf 12.1, ...)
    !               2. gnuplot 4.5 and higher (other previous version can be used
    !    author:    mohammad rahmani
    !               chem eng dep., amirkabir uni. of tech
    !               tehran, ir
    !               url: aut.ac.ir/m.rahmani
    !               email: m[dot]rahmani[at]aut[dot]ac[dot]ir


! this file demonstrate the capability of ogpf module
! an object based fortran interface to gnuplot


    ! version:  0.12
    ! date:     feb 9th, 2012
    ! new examples for semilogx, semilogy, loglog
    ! new set options method

    ! version:  0.11
    ! date:     feb 9th, 2012

!}}}
module demo
   use ogpf
  
   
   
    contains
!}}}
!subroutine choose_exmp{{{
    !  select a number between  1 and 20
    subroutine choose_exmp(i)
     select case(i)
 case(1)
 call exmp02()
 case(2)
 call exmp03()
 case(3)
 call exmp04()
 case(4)
 call exmp05()
 case(5)
 call exmp06()
 case(6)
 call exmp07()
 case(7)
 call exmp08()
 case(8)
 call exmp09()
 case(9)
 call exmp10()
 case(10)
 call exmp11()
 case(11)
 call exmp12()
 case(12)
 call exmp13()
 case(13)
 call exmp14()
 case(14)
 call exmp15()
 case(15)
 call exmp16()
 case(16)
 call exmp17()
 case(17)
 call exmp18()
 case(19)
 call exmp20()
 case(20)
 call exmp21()
 case default
 print *,"no number, please choose another one"
 end select


    end subroutine choose_exmp
!}}}
!subroutine data_line{{{
   subroutine data_line(x,y,n)
       ! draw the data line with * representing the point
        type(gpf):: gp
        integer  n
        real(8):: x(n),y(n)

        ! sample i. plot on the screen
        !annotation, set title, xlabel, ylabel
        call gp%title('example 1. a simple xy plot')
        call gp%xlabel('my x axis ...')
        call gp%ylabel('my y axis ...')
        call gp%options('set style data linespoints')

        !call plot to draw a vector against a vector of data
        call gp%plot(x,y)
    end subroutine data_line
!}}}
!subroutine exmp02{{{
 subroutine exmp02()
        !...............................................................................
        !example 2: set line specification and legends
        !...............................................................................
        type(gpf):: gp
        integer, parameter:: n=17
        real(8):: x(n)
        real(8):: y(n)
        ! input data
        x=dble([- 8,- 7,- 6,- 5,- 4,- 3,- 2,- 1,0,1,2,3,4,5,6,7,8])
        y=dble([66,51,38,27,18,11,6,3,2,3,6,11,18,27,38,51,66])

        ! sample i. plot on the screen
        !annotation, set title, xlabel, ylabel
        call gp%title('example 2. a simple xy plot')
        call gp%xlabel('my x axis ...')
        call gp%ylabel('my y axis ...')

        !call plot to draw a vector against a vector of data
        !the last argument defines the line specification
        call gp%plot(x,y,'with linespoints lt 2 pt 4')
 end subroutine exmp02
!}}}
!{{{subroutine exmp03
    subroutine exmp03()
        !...............................................................................
        ! example 3: plot several data set at the same time
        !...............................................................................
        type(gpf):: g
        integer, parameter:: n=50
        integer, parameter:: m=65
        real(8):: x(n)
        real(8):: y(n)
        real(8):: xv(m)
        real(8):: yv(m)
        real(8), parameter :: pi=4.d0*atan(1.d0)
        ! input data
        x=linspace(-pi,pi,n)  !linspace is a utility function in ogpf module
        y=sin(x)              !linspace(a,b,n) create a linear vector in [a,b] with n elements

        xv=linspace(0.d0, 2.d0*pi,m)
        yv=cos(2.d0*xv)

        ! annotation, set title, xlabel, ylabel
        call g%title('example 3. plot two data series using gnuplot')
        call g%xlabel(' x axis ...')
        call g%ylabel(' y axis ...')

        ! sample 1: plot to draw two set of data
        call g%plot(x,y,'title "sin"',xv,yv,'title "cos"') !linespec is empty

        !sample 2: use keyword arguments to plot the same example as above
        call g%title('example 3. another plot using keyword arguments...')
        call g%plot(x1=x,y1=y,ls1='',x2=xv,y2=yv,ls2='')

        ! sample 3: an xy plot with line specification no legends
        call g%title('example 3. another plot using keyword arguments...')
        call g%plot(x,y,'title "sin(x)" lt 6', xv,yv,'title "cos(2x)" with points pt 7 lc rgb "#993300"')

    end subroutine exmp03
!}}}
!subroutine exmp04{{{
      subroutine exmp04()
        !...............................................................................
        ! example 4: plot four data series, the maximum number of data series can be plotted!
        !...............................................................................
        type(gpf):: gp
        integer, parameter:: n=50
        integer, parameter:: m=65
        real(8):: x(n)
        real(8):: y(n)
        real(8)::  xv(m)
        real(8):: yv(m)
        real(8), parameter :: pi=4.d0*atan(1.d0)
        ! input data
        x=linspace(-pi,pi,n)  !linspace is a utility function from module utils
        y=sin(x)

        xv=linspace(0.d0, 2.d0*pi,m)
        yv=cos(2.d0*xv)
        !           this is the maximum number of plot can be drawn at the same time
        !           if you have more data see, you can plot can be used with matrices!
        call gp%title('example 4. plot four data sets using gnuplot')
        call gp%plot(x,y, 'title "sin(x)"', &
        xv,yv, 'with lp lt 6 title "cos(2x)"', &
        xv, 2.d0*yv, 'title "2cos(2x)" lt 7', &
        xv, 0.5d0*yv, 'title "0.5cos(2x)" with points pt 8')

        ! another example with keyboard arguments
        call gp%filename('exmp04.plt') ! the plot commands are saved into exmp04.plt for further use
        call gp%plot(x1=x,y1=y,x2=xv,y2=yv)

    end subroutine exmp04
!}}}
!subroutine exmp05{{{
    subroutine exmp05()
        !...............................................................................
        ! example 5: use line style and legends
        !...............................................................................
        type(gpf):: gplot
        integer, parameter :: n=50
        real(8), parameter :: pi=4.d0*atan(1.d0)
        real(8)            :: x(n)
        real(8)            :: ys(n)
        real(8)            :: yc(n)
        real(8)            :: ysc(n)
        ! input data
        x=linspace(-2.d0*pi,2.d0*pi,n)  !linspace is a utility function from module utils
        ys=sin(x)
        yc=exp(-0.1d0*x)*cos(x)
        ysc=sin(x/2.d0)*cos(2.d0*x)

        ! annotation, set title, xlabel, ylabel
        call gplot%title('example 5. a sample with style and legends')
        call gplot%xlabel('x, rad')
        call gplot%ylabel('y, dimensionless')

        ! plot to draw three set of data
        call gplot%plot(x,ys,'title "sin" with lines lt 5 lc rgb "#0008b0"', &
                        x,yc,'title "cos" with points lt 6 lc rgb "#ff1100"', &
                        x,ysc,'title "sin(x/2)cos(2x)" with lp lt 7 lc rgb "#00aa04"' )

    end subroutine exmp05
!}}}
!subroutine exmp06{{{
   subroutine exmp06()
        !...............................................................................
        ! example 6: plot a single point along with a series of data
        !...............................................................................
        type(gpf):: gplot
        integer, parameter:: n=125

        real(8):: x(n)
        real(8):: y(n)

        real(8), parameter :: pi=4.d0*atan(1.d0)
        ! input data
        x=linspace(0.d0,pi*2.d0,n)  !linspace is a utility function from module utils
        y=sin(6.d0*x)*exp(-x)


        ! annotation, set title, xlabel, ylabel
        call gplot%title('example 6. a sample shows sin(x) and its zero on the plot')
        call gplot%xlabel('x, rad')
        call gplot%ylabel('sin(x), dimensionless')
        call gplot%options('set grid')

        ! plot to draw three set of data
        call gplot%plot(x,y,'title "sin(x)" with lines lt 3', &
                        [pi],[0.d0],'title "zero" with points pt 7 ps 2 lc rgb "#ff0000"')
    end subroutine exmp06
!}}}
!subroutine exmp07{{{
    subroutine exmp07()
        !plot a matrix against a vector
        type(gpf):: matplot
        integer, parameter:: n=25
        real(8), parameter :: pi=4.d0*atan(1.d0)
        real(8):: x(n)
        real(8):: y(n,6)

        !create data
        x=linspace(-pi,pi,n)
        y(:,1)=sin(x)
        y(:,2)=cos(x)
        y(:,3)=cos(0.5d0*x)
        y(:,4)=sin(0.5d0*x)
        y(:,5)=sin(x)*cos(x)
        y(:,6)=sin(x)*exp(-x**2)


        !draw the matrix y againest vector x
        call matplot%title('example 7. plotting a matrix against a vector')
        call matplot%xlabel ('my x axis')
        call matplot%ylabel ('my y axis')
    !    call matplot%plot(x, y)    2015-05-21 17:35:20    
    end subroutine exmp07
!}}}
!subroutine exmp08{{{
    subroutine exmp08()
        !plot a matrix against a vector
        type(gpf):: matplot
        integer, parameter:: n=25
        real(8):: tf
        real(8):: vo
        real(8):: g
        real(8):: t(n)
        real(8):: y(n,6)

        !create data
        tf=10.d0
        g=32.d0;
        t=linspace(0.d0,tf,n)
        vo=150.d0;
        y(:,1)=vo*t-0.5d0*g*t**2
        vo=125.d0;
        y(:,2)=vo*t-0.5d0*g*t**2
        vo=100.d0;
        y(:,3)=vo*t-0.5d0*g*t**2
        vo=75.d0;
        y(:,4)=vo*t-0.5d0*g*t**2
        vo=50.d0;
        y(:,5)=vo*t-0.5d0*g*t**2
        vo=25.d0;
        y(:,6)=vo*t-0.5d0*g*t**2

        !draw the matrix y againest vector x
        call matplot%title('example 8. plotting a matrix against a vector')
        call matplot%xlabel ('t, sec')
        call matplot%ylabel ('y, feet')
        call matplot%options('set xrange[0:10];set yrange [0:400];')
        !  call matplot%plot(t, y)    2015-05-21 17:36:06    
    end subroutine exmp08
!}}}
!subroutine exmp09{{{
   subroutine exmp09()
        !...............................................................................
        ! example 09: use gnuplot interactively
        !...............................................................................
        type(gpf):: xyplot
        integer, parameter:: n=17
        real(8):: x(n)
        real(8):: y(n)
        ! input data
        x=dble([-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8])
        y=dble([66,51,38,27,18,11,6,3,2,3,6,11,18,27,38,51,66])

        ! annotation, set title, xlabel, ylabel
        call xyplot%title('example 9. a simple xy plot in an interactive session')
        call xyplot%xlabel('my x axis ...')
        call xyplot%ylabel('my y axis ...')

        ! keep gnuplot open for an interactive session
        ! you can enter any valid gnuplot commands and finally
        ! type q to exit
        call xyplot%hold('on')
        ! call plot to draw a vector against a vector of data
        call xyplot%plot(x,y)

     end subroutine exmp09
!}}}
!subroutine exmp10{{{
    subroutine exmp10()
        !use gnuplot options
        type(gpf):: mp
        real(8):: x(10)
        real(8):: y(10)
        !option is a string of 320 character and can set all kind of gnuplot
        !global options

        call mp%options(&
               "set logscale xy;&
               &set xrange [0.1:100];&
               &set yrange [0.01:10000];")

        x=linspace(0.1d0,100d0,10); y=x**2;
        call mp%title("example 10. x vs. x^2")
        call mp%plot(x,y)
        call mp%reset()
        call mp%plot(x,2*y)
    end subroutine exmp10
!}}}
!subroutine exmp11{{{
!...............................................................................
! example 11: a simple polar plot
!...............................................................................
    subroutine exmp11()
        type(gpf):: gp
        integer, parameter :: n=75
        real(8):: t(n)
        real(8):: r(n)
        real(8):: pi=4.d0*atan(1.d0)
    !1. reset gplot
     call gp%reset()

    !2. set option, and set plot as polar
        call gp%options("&
            &set polar;&
            &set trange [-pi/2:pi/2]")

    ! 3. create data
        t=linspace(-pi/2.d0,pi/2.d0,n)
        r=sin(3.d0*t)

    !annotation, set title, xlabel, ylabel
    call gp%title("example 11: simple polar plot")
    call gp%xlabel("x,...")
    call gp%ylabel("y,...")

    !call plot method
    call gp%plot(t,r)

    end subroutine exmp11
!}}}
!subroutine exmp12{{{
!...............................................................................
! example 12: a simple plot with logarithmic x axis
!...............................................................................
    subroutine exmp12()
        type(gpf):: gp
        integer, parameter :: n=75
        real(8):: x(n)
        real(8):: y(n)


    ! 1. create data
        x=linspace(0.1d0,10.d0,n)
        y=5.d0*x**3+4.d0*x**2+3.d0*x+1.d0

    !annotation, set title, xlabel, ylabel
    call gp%title("example 12: a semi-log x plot")
    call gp%xlabel("x,logarithmic scale")
    call gp%ylabel("y, normal scale")

    !call plot method
    call gp%semilogx(x,y)


    end subroutine exmp12
!}}}
!subroutine exmp13{{{
!...............................................................................
! example 13: a simple plot with logarithmic y axis
!...............................................................................
    subroutine exmp13()
        type(gpf):: gp
        integer, parameter :: n=75
        real(8):: x(n)
        real(8):: y(n)


    ! 1. create data
        x=linspace(0.1d0,10.d0,n)
        y=5.d0*x**3+4.d0*x**2+3.d0*x+1.d0

    !annotation, set title, xlabel, ylabel
    call gp%title("example 13: a semi-log y plot")
    call gp%ylabel("y,logarithmic scale")
    call gp%xlabel("x, normal scale")

    !call plot method
    call gp%semilogy(x,y)

    end subroutine exmp13
!}}}
!subroutine exmp14{{{
!...............................................................................
! example 13: a simple plot with logarithmic xy axes
!...............................................................................
    subroutine exmp14()
        type(gpf):: gp
        integer, parameter :: n=75
        real(8):: x(n)
        real(8):: y(n)
        real(8):: pi=4.d0*atan(1.d0)

    ! 1. create data
        x=exp(linspace(0.d0,2.d0*pi,n))
        y=50.d0+exp( 3.d0* linspace(0.d0,2.d0*pi,n) )

    !annotation, set title, xlabel, ylabel
    call gp%title("example 14: a semi-log y plot")
    call gp%xlabel("x,logarithmic scale")
    call gp%ylabel("y,logarithmic scale")
    !set grid on
    call gp%options('set grid xtics ytics mxtics')


    !call plot method
    call gp%loglog(x,y)

    end subroutine exmp14
!}}}
!subroutine exmp15{{{
     subroutine exmp15()
        !a demo to plot a function using fplot method
        type(gpf):: gp
        real(8):: pi=4.d0*atan(1.d0)

        call gp%title("example 15. a fplot example")
        call gp%xlabel("x...")
        call gp%ylabel("y...")
        call gp%filename('example15.plt') !save the results into a file
        call gp%fplot(myfun,[0d0,15.d0*pi],150)
        print*, 'plot commands were written in example15.plt successfully'
        print*, 'open gnuplot and load this script file to plot the results!'
    end subroutine exmp15

    function myfun(x)
        !the function to be plotted
        !see example 15
        real(8), intent(in) :: x
        real(8):: myfun
        myfun=x*sin(x)
    end function myfun
!}}}
!subroutine exmp16{{{
    subroutine exmp16()
        !demo for gnuplot script
        !a script is a gnuplot set of commands to be sent through a file for gnuplot

        type(gpf):: gp

        call gp%title("example 16. script files")
        !example 16-1
        call gp%script("set xrange[-pi:pi]; set title 'script plot'; plot sin(x)")

        !example 16-2
        call gp%script("splot [-2:2][-2:2] exp(-(x**2 + y**2))*cos(x/4)*sin(y)*cos(2*(x**2+y**2))")

        !example 16-3
        call gp%script("set polar;plot t*sin(t);&
                        &set trange [-2*pi:2*pi]; set rrange [0:3];plot t*sin(t)")

    end subroutine exmp16
!}}}
!subroutine exmp17{{{
    subroutine exmp17()
        !demo for gnuplot script
        type(gpf):: gp
        call gp%title("example 17. another sample script file")
        call gp%script('splot x**2+y**2')

    end subroutine exmp17
!}}}
!subroutine exmp18{{{

    subroutine exmp18()
        !use gnuplot script
        !to send a special script file to gnuplot
        !the file is an external file here is called "simple.dem"
        type(gpf):: gp
        call gp%title("example 18. running an external script file")
        call gp%script("load 'simple.plt'")

    end subroutine exmp18
!}}}
!subroutine exmp20{{{
! 3d plots
 subroutine exmp20()
        !a simple 3d plot
        type(gpf):: gp
        real(8), allocatable:: x(:,:)
        real(8), allocatable:: y(:,:)
        real(8), allocatable:: z(:,:)
        real(8):: a=0.5d0
        real(8):: b=2.0d0
        integer:: m
        integer:: n
        call meshgrid(x,y,dble([-10,10,2]),dble([-10,10,2]))
        m=size(x,1)
        n=size(x,2)
        allocate( z(m,n) )
        z=(x**2/a - y**2/b)

        !#annotation, label colors are also set!
        call gp%filename('example20.plt')
        !here options has been called several times instead of one time with a long string
        call gp%options("set style data lines")
        call gp%options("unset key;set contour base")
        call gp%options("set xyplane relative 0.1")

        call gp%title('example 20: simple 3d plot using splot')
        call gp%xlabel('x-axis,...')
        call gp%ylabel('y-axis,...')
        call gp%zlabel('z-axis,...')

        !plot the 3d data
        call gp%surf(z,lspec= 'title "x^2/a - y^2/b"')
    end subroutine exmp20
!}}}
!subroutine exmp21{{{
    subroutine exmp21()
         !another simple 3d plot
        type(gpf):: gp
        real(8), allocatable:: x(:,:)
        real(8), allocatable:: y(:,:)
        real(8), allocatable:: z(:,:)
        integer:: m
        integer:: n
        call meshgrid(x, y, [real(8)::-15.,15.,0.75] )
        m=size(x,1)
        n=size(x,2)
        allocate( z(m,n) )
        !z=x*x*exp(-x*x) * y*y*exp(-y*y)
        !z=x*y**3-y*x**3
        ! z=-1/(x**2+y**2+0.25)
        z=sin(sqrt(x**2+y**2))/sqrt(x**2+y**2+0.025)

        call gp%title('example 21: simple 3d plot using splot')
        call gp%xlabel('x-axis,...')
        call gp%ylabel('y-axis,...')
        call gp%zlabel('z-axis,...')

        !plot the 3d data
        call gp%surf(x,y,z,lspec= 'title "sample 3d plot"')

    end subroutine exmp21
!}}}
end module demo
