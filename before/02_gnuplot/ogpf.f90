!README{{{
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

    ! revision history

    ! version:  0.12
    ! date:     feb 9th, 2012
    !   minor corrections
    !   new semilogx, semilogy, loglog methods
    !   new options method, allow to be called several times to set the gnuplot options



    ! version:  0.11
    ! date:     feb 9th, 2012
    !   minor corrections
    !   use of newuint specifier from fortran 2008
    !   added configuration parameters
    !   extra procedures have been removed
    !   temporary file is now deleted using close(...,status='delete')

    !
    ! version:  0.1
    ! date:     jan 5th, 2012
    !  first object-based version
!}}}

    module ogpf
    use ifport
    implicit none
    private
    public linspace, meshgrid
    integer,public  ::  ii,jj,kk

    ! configuration parameters
    character(len=*), parameter             :: gnuplot_terminal        ='pdf'             !output terminal
    character(len=*), parameter             :: gnuplot_font_name       ='calibri,10'     !font
    character(len=*), parameter             :: gnuplot_plot_size       ='360,240'        !plot window size
    character(len=*), parameter             :: gnuplot_datastyle       ='linespoints'    !data style
    character(len=*), parameter             :: gnuplot_output_filename ='ogpf_temp_script_file.plt'   !temporary file for output

!   can be adjusted as required
    integer, parameter                      :: len_options   = 320
    integer, parameter                      :: len_labels    = 50
    integer, parameter                      :: len_msg       = 80

    !type gpf{{{
    type, public                            :: gpf
        private
        character(len=len_labels)           :: txtplottitle   = ""
        character(len=len_labels)           :: txtxlabel      = ""
        character(len=len_labels)           :: txtylabel      = ""
        character(len=len_labels)           :: txtzlabel      = ""
        character(len=len_options)          :: txtoptions     = ""  !four 80 characters lines
        logical                             :: displayplot    =.true.
        logical                             :: persist        = .false.
        logical                             :: hasplottitle   =.false.
        logical                             :: hasxlabel      =   .false.
        logical                             :: hasylabel      =   .false.
        logical                             :: haszlabel      =   .false.
        logical                             :: hasoptions     =  .false.
        real(8)                             :: xrange(2)
        real(8)                             :: yrange(2)
        real(8)                             :: range(2)

        character(len=len_msg)              :: msg        = ""   !message from plot procedures
        integer                             :: status     = 0 !status from plot procedures
        character(len=25)                   :: txtfilename=gnuplot_output_filename
        character(len=8)                    :: plotscale

    contains
    procedure, pass, public                 :: options  => set_options
    procedure, pass, public                 :: title    => set_plottitle
    procedure, pass, public                 :: xlabel   => set_xlabel
    procedure, pass, public                 :: ylabel   => set_ylabel
    procedure, pass, public                 :: zlabel   => set_zlabel
    procedure, pass, public                 :: filename => set_filename
    procedure, pass, public                 :: reset    => reset_to_defaults
    procedure, pass, public                 :: hold     => set_persist
    procedure, pass, public                 :: surf     => splot
    procedure, pass, public                 :: script   => gnuplotscript
    procedure, pass, public                 :: fplot
    procedure, pass, public                 :: semilogx
    procedure, pass, public                 :: semilogy
    procedure, pass, public                 :: loglog

    procedure, pass, public                 :: plot   => plot2d_vector_vs_vector
    procedure, pass, private                :: plot2d_vector_vs_vector
    procedure, pass, private                :: plot2d_matrix_vs_vector
    end type gpf
!}}}


    contains
!subroutine fplot{{{
subroutine fplot(this, func,xrange,np)
    ! fplot, plot a function in the range xrange=[xmin, xamx] with np points
    ! if np isnot sent, then np=50 is assumed!
    ! func is the name of function to be plotted

    class(gpf):: this
    interface
function func(x)
    real(8), intent(in) :: x
    real(8) :: func
    end function func
    end interface
    real(8), intent(in) :: xrange(2)
    integer, optional, intent(in):: np

    integer:: n
    integer:: i
    integer:: alloc_err
    real(8), allocatable :: x(:)
    real(8), allocatable :: y(:)

    if (present(np)) then
        n=np
    else
        n=50
    end if
    allocate(x(1:n), y(1:n), stat=alloc_err)
    if (alloc_err /=0) then
        stop "allocation error in fplot procedure..."
    end if
    !create set of xy data
    x=linspace(xrange(1),xrange(2), n)
    y=[ (func(x(i)), i=1, n) ]

    call plot2d_vector_vs_vector(this,x,y)


    if (allocated(x)) deallocate(x)
    if (allocated(y)) deallocate(y)


    end subroutine fplot
!!}}}
!subroutine semilogx{{{
    !..............................................................................
subroutine semilogx(this, x1, y1, ls1, &
    x2, y2, ls2, &
    x3, y3, ls3, &
    x4, y4, ls4  )
    !..............................................................................
    !   this procedure is the same as plotxy with logarithmic x axis
    !..............................................................................
    class(gpf)                                    :: this
    ! input vector
    real(8),  intent(in)                          :: x1(:)
    real(8),  intent(in), optional                :: y1(:)
    real(8),  intent(in), dimension(:), optional  :: x2
    real(8),  intent(in), dimension(:), optional  :: y2
    real(8),  intent(in), dimension(:), optional  :: x3
    real(8),  intent(in), dimension(:), optional  :: y3
    real(8),  intent(in), dimension(:), optional  :: x4
    real(8),  intent(in), dimension(:), optional  :: y4

    character(len=*),  intent(in), optional       ::  ls1
    character(len=*),  intent(in), optional       :: ls2
    character(len=*),  intent(in), optional       :: ls3
    character(len=*),  intent(in), optional       ::  ls4

    this%plotscale='semilogx'
    call plot2d_vector_vs_vector(this, x1, y1, ls1, x2, y2, ls2, x3, y3, ls3, x4, y4, ls4  )
    ! set the plot scale as normal. it means log scale is off
    this%plotscale='normal'

    end subroutine semilogx
!!}}}
!subroutine semilogy{{{
    !..............................................................................
subroutine semilogy(this, x1, y1, ls1, &
    x2, y2, ls2, &
    x3, y3, ls3, &
    x4, y4, ls4  )
    !..............................................................................
    !   this procedure is the same as plotxy with logarithmic y axis
    !..............................................................................
    class(gpf):: this
    ! input vector
    real(8),  intent(in)            :: x1(:)
    real(8),  intent(in), optional  :: y1(:)
    character(len=*),  intent(in), optional   ::  ls1

    real(8),  intent(in), dimension(:), optional  :: x2
    real(8),  intent(in), dimension(:), optional  :: y2
    character(len=*),  intent(in), optional       :: ls2

    real(8),  intent(in), dimension(:), optional  :: x3
    real(8),  intent(in), dimension(:), optional  :: y3
    character(len=*),  intent(in), optional       :: ls3

    real(8),  intent(in), dimension(:), optional  :: x4
    real(8),  intent(in), dimension(:), optional  :: y4
    character(len=*),  intent(in), optional   ::  ls4

    this%plotscale='semilogy'
    call plot2d_vector_vs_vector(this, x1, y1, ls1, x2, y2, ls2, x3, y3, ls3, x4, y4, ls4  )
    ! set the plot scale as normal. it means log scale is off
    this%plotscale='normal'


    end subroutine semilogy
!!}}}
!subroutine loglog{{{
    !..............................................................................
subroutine loglog(this, x1, y1, ls1, &
    x2, y2, ls2, &
    x3, y3, ls3, &
    x4, y4, ls4  )
    !..............................................................................
    !   this procedure is the same as plotxy with logarithmic x axis
    !..............................................................................
    class(gpf):: this
    ! input vector
    real(8),  intent(in)            :: x1(:)
    real(8),  intent(in), optional  :: y1(:)
    character(len=*),  intent(in), optional   ::  ls1

    real(8),  intent(in), dimension(:), optional  :: x2
    real(8),  intent(in), dimension(:), optional  :: y2
    character(len=*),  intent(in), optional       :: ls2

    real(8),  intent(in), dimension(:), optional  :: x3
    real(8),  intent(in), dimension(:), optional  :: y3
    character(len=*),  intent(in), optional       :: ls3

    real(8),  intent(in), dimension(:), optional  :: x4
    real(8),  intent(in), dimension(:), optional  :: y4
    character(len=*),  intent(in), optional   ::  ls4

    this%plotscale='loglog'
    call plot2d_vector_vs_vector(this, x1, y1, ls1, x2, y2, ls2, x3, y3, ls3, x4, y4, ls4  )
    ! set the plot scale as normal. it means log scale is off
    this%plotscale='normal'

    end subroutine loglog
!!}}}
!subroutine plot2d_vector_vs_vector{{{
    !..............................................................................
subroutine plot2d_vector_vs_vector(this, x1, y1, ls1, &
    x2, y2, ls2, &
    x3, y3, ls3, &
    x4, y4, ls4  )
    !..............................................................................
    !   this procedure plots:
    !   1. a vector against another vector (xy plot)
    !   2. a vector versus its element indices.
    !   3. can accept up to 4 data sets as x,y pairs!
    ! arguments
    ! xi, yi vectors of data series,
    ! lsi a string maximum 80 characters containing the line specification, legends, ...

    class(gpf):: this
    ! input vector
    real(8),  intent(in)            :: x1(:)
    real(8),  intent(in), optional  :: y1(:)
    character(len=*),  intent(in), optional   ::  ls1

    real(8),  intent(in), dimension(:), optional  :: x2
    real(8),  intent(in), dimension(:), optional  :: y2
    character(len=*),  intent(in), optional       :: ls2

    real(8),  intent(in), dimension(:), optional  :: x3
    real(8),  intent(in), dimension(:), optional  :: y3
    character(len=*),  intent(in), optional       :: ls3

    real(8),  intent(in), dimension(:), optional  :: x4
    real(8),  intent(in), dimension(:), optional  :: y4
    character(len=*),  intent(in), optional   ::  ls4

    !   local variables
    !----------------------------------------------------------------------

    integer:: nx1
    integer:: ny1
    integer:: nx2
    integer:: ny2
    integer:: nx3
    integer:: ny3
    integer:: nx4
    integer:: ny4
    integer:: number_of_plots
    character(len=3)::  plottype
    integer:: file_unit
    integer:: i
    character(len=80)::  pltstring(4)  !four 80 character lines

    !initialize variables
    plottype=''
    pltstring=''
    !   check the input
    nx1=size(x1)
    if ((present(y1) )) then
        ny1=size(y1)
        if (checkdim(nx1,ny1)) then
            plottype='xy1'
            number_of_plots=1
        else
            print*, 'gpf error: length of x1 and y1 doesnot match'
            return
        end if
    else !plot only x againest its element indices
        plottype='xi'
        number_of_plots=1
    end if

    !process line spec for first data set if present
    call process_linespec(1, pltstring(1),ls1)

    if (present(x2) .and. present (y2)) then
        nx2=size(x2)
        ny2=size(y2)
        if (checkdim(nx2,ny2)) then
            plottype='xy2'
            number_of_plots=2
        else
            return
        end if
        !process line spec for 2nd data set if present
        call process_linespec(2, pltstring(2),ls2)
    end if

    if (present(x3) .and. present (y3)) then
        nx3=size(x3)
        ny3=size(y3)
        if (checkdim(nx3,ny3)) then
            plottype='xy3'
            number_of_plots=3
        else
            return
        end if
        !process line spec for 3rd data set if present
        call process_linespec(3, pltstring(3),ls3)
    end if

    if (present(x4) .and. present (y4)) then
        nx4=size(x4)
        ny4=size(y4)
        if (checkdim(nx4,ny4)) then
            plottype='xy4'
            number_of_plots=4
        else
            return
        end if
        !process line spec for 4th data set if present
        call process_linespec(4, pltstring(4),ls4)
    end if

    ! open the output file
    open ( unit= file_unit, file = this%txtfilename, status = 'replace',iostat = this%status )

    if (this%status /= 0 ) then
        this%msg= "gpf error: cannot open file for output"
        print*, this%msg
        return !an error has been occurred
    end if

    ! write plot title, axis labels and other annotations
    call processcmd(this, file_unit)

    ! write plot command and line styles and legend if any
    if (number_of_plots ==1) then
        write ( file_unit, '(a)' ) trim(pltstring(1))
    else
        write ( file_unit, '(a)' ) ( trim(pltstring(i)) // ' \' , i=1, number_of_plots-1)
        write ( file_unit, '(a)' )   trim(pltstring(number_of_plots))
    end if
    ! write xy data into file
    select case (plottype)
    case ('xi')
        do i = 1, size(x1)
            !         write ( file_unit, * ) x1(i)
        end do
        !     write ( file_unit, '(a)' ) 'e'  !end of jth set of data
    case ('xy1')
        call write_xydata(file_unit,nx1,x1,y1)
    case ('xy2')
        call write_xydata(file_unit,nx1,x1,y1)
        call write_xydata(file_unit,nx2,x2,y2)
    case ('xy3')
        call write_xydata(file_unit,nx1,x1,y1)
        call write_xydata(file_unit,nx2,x2,y2)
        call write_xydata(file_unit,nx3,x3,y3)
    case ('xy4')
        call write_xydata(file_unit,nx1,x1,y1)
        call write_xydata(file_unit,nx2,x2,y2)
        call write_xydata(file_unit,nx3,x3,y3)
        call write_xydata(file_unit,nx4,x4,y4)
    end select


    if (this%persist .and. this%displayplot) then
        close ( unit = file_unit )
        write(*,*)
        write(*,*) 'this is gnuplot interactive mode!'
        write(*,*)  'type q and press enter to exit'
        write(*,*)
    else
        write ( file_unit, '(a)' ) 'pause -1 "press a key to continue..."'
        write ( file_unit, '(a)' ) 'q'
        close ( unit = file_unit )
    end if
    !   now plot the results
    if (this%displayplot) then
        call write2gnuplot(this%txtfilename,this%persist)
    else
        this%displayplot=.true. !reset display plot to its value
        this%txtfilename=gnuplot_output_filename
    end if
    !
    !: end of plot2d_vector_vs_vector
    end subroutine plot2d_vector_vs_vector
!!}}}
!subroutine splot{{{
    !'++++++++++++++++++++++++++++++++++++++++++++++++++
    !..............................................................................
subroutine splot(this, x, y, z, lspec)
    !..............................................................................
    class(gpf):: this
    ! input vector
    real(8),  intent(in)            :: x(:,:)
    real(8),  intent(in), optional  :: y(:,:)
    real(8),  intent(in), optional  :: z(:,:)
    character(len=*),  intent(in), optional   ::  lspec

    !   local variables
    !----------------------------------------------------------------------

    integer:: ncx
    integer:: nrx
    integer:: file_unit
    integer:: i
    integer:: j
    logical:: xyz_data
    character(len=80)::  pltstring

    pltstring=''
    !   check the input data
    ncx=size(x,dim=2)
    nrx=size(x,dim=1)
    if (present(y) .and. present(z)) then
        xyz_data=.true.
    elseif (present(y)) then
        print*, "gpf error: z matrix was not sent to 3d plot routine"
        return
    else
        xyz_data=.false.
    end if

    ! open the output file
    open ( unit = file_unit, file = this%txtfilename, status = 'replace',iostat = this%status )
    if (this%status /= 0 ) then
        this%msg= "gpf error: cannot open file for output"
        print*, this%msg
        return !an error has been occurred
    end if


    ! set the plot scale as normal. it means log scale is off
    this%plotscale='normal'
    ! write titles and other annotations
    call processcmd(this, file_unit)

    if ( present(lspec) ) then
        if (hastitle(lspec)) then
            pltstring='splot "-" '//trim(lspec)
        else
            pltstring='splot "-" notitle '//trim(lspec)
        end if
    else
            pltstring='splot "-" notitle '
        end if

    write ( file_unit, '(a)' ) trim(pltstring)

    ! write xy data into file
    write ( file_unit, '(a)' ) '#data x y z'


    if (xyz_data) then
        do j=1,ncx
            do i=1, nrx
                write ( file_unit, * ) x(i,j), y(i,j), z(i,j)
            enddo
            write( file_unit, '(a)' )  !put an empty line
        enddo
        write ( file_unit, '(a)' ) 'e'  !end of set of data
    else !only z has been sent (i.e. single matrix data)
        do j=1,ncx
            do i=1, nrx
                write ( file_unit, * ) i, j, x(i,j)
            enddo
            write( file_unit, '(a)' )  !put an empty line
        enddo
        write ( file_unit, '(a)' ) 'e'  !end of set of data
    end if

    if (this%persist .and. this%displayplot) then
        close ( unit = file_unit )
        write(*,*)
        write(*,*) 'this is gnuplot interactive mode!'
        write(*,*)  'type q and press enter to exit'
        write(*,*)
    else
        write ( file_unit, '(a)' ) 'pause -1 "press a key to continue..."'
        write ( file_unit, '(a)' ) 'q'
        close ( unit = file_unit )
    end if
    !   now plot the results
    if (this%displayplot) then
        call write2gnuplot(this%txtfilename,this%persist)
    else
        this%displayplot=.true. !reset display plot to its value
        this%txtfilename='ogpf_temp_script_file.plt'
    end if
    !
    !: end of splot
    end subroutine splot
!!}}}
!subroutine write2gnuplot{{{
    !..............................................................................
subroutine write2gnuplot(filename, persist)
    !..............................................................................
    !   this subroutine call gnuplot through system command
    implicit none
    character(len=*):: filename
    logical, intent(in) :: persist
    integer   ::  icolor

    !local vars
    integer :: lun, ios

    if  (persist) then
        !fortran standard recommend to use call 
        !execute_command_line to invoke another program
        !from within fortran, the old method was to use call system
        !here by default fortran standard is used, 
        !if you have a compiler does not support
        !call execute_command_line, uncomment the following call system
        !call system('/usr/bin/gnuplot -persist '//filename)              !obsolete method, use with old compilers
        icolor=system('/usr/bin/gnuplot -persist '//filename) !fortran standard
    else
        !call system ('gnuplot '//filename) !obsolete method, use with old compilers
        icolor=system('/usr/bin/gnuplot '//filename)
    end if

!   the following lines actually do the deletion of temporary file
!   this method is used to have a portable code!

    open ( unit = lun, file = filename, status = 'old',iostat = ios )
    if (ios /= 0 ) then
        print*, "gpf error: cannot open file for output"
        return !an error has been occurred
    end if
    close(lun,status='delete')  !delete file

    end subroutine write2gnuplot
!!}}}
!subroutine writestring{{{
    !..............................................................................
subroutine writestring(str, file_unit)
    !..............................................................................
    !writestring, separate a string using ";" delimiters
    !and write each statement in a separate line in the file indicated by file_unit
    character(len=*), intent(in):: str
    integer, intent(in)   :: file_unit
    character, parameter :: delimiter=';'
    integer::n
    integer:: m
    integer:: k


    k=len_trim(str) !length with removed trailing blanks
    n=scan(str,delimiter)
    if (n==0) then  !this is a single statement
        write(file_unit,'(a)') str
        return
    end if
    m=1
    do while (n/=0 .and. m<k)
        if (n/=1) then
            write(file_unit,'(a)') str(m:m+n-2)
        end if
        m=n+m
        n=scan(str(m:k),delimiter)
    end do
    if (m<k) then !write the last statement
        write(file_unit,'(a)') str(m:k)
    end if
    end subroutine writestring
!!}}}
!subroutine processcmd{{{
    !..............................................................................
subroutine processcmd(this, file_unit)
    !..............................................................................
    !   this subroutine writes all the data into plot file
    !   to be read by gnuplot

    class(gpf), intent(in)   :: this
    integer, intent(in)     :: file_unit
    ! the following lines set the gnuplot terminal
    ! the data style to lines+symbols:linespoints
    ! can be overwritten by options
    write ( file_unit, '(a)' ) 'set term ' //gnuplot_terminal// &
                               ' font ' // '"'// gnuplot_font_name // '"' // &
                               ' size '// gnuplot_plot_size          !set output terminal
    write ( file_unit, '(a)' ) 'set style data '//gnuplot_datastyle !set data style

     ! write options
    if ( this%hasoptions ) then
        call writestring(trim(this%txtoptions),file_unit)
    end if
    ! check with plot scale: i.e normal, logx, logy, or log xy
    select case (this%plotscale)
    case ('semilogx')
        write ( file_unit, '(a)' ) 'set logscale  x'
    case ('semilogy')
        write ( file_unit, '(a)' ) 'set logscale  y'
    case ('loglog')
        write ( file_unit, '(a)' ) 'set logscale  xy'
    case default !for normal xy plot or 3d plots
        !pass
    end select
    !   write the plot options to script file
    if (this%hasplottitle) then
        write ( file_unit, '(a)' ) 'set title  "' // trim(this%txtplottitle)// '"'
    end if
    if (this%hasxlabel) then
        write ( file_unit, '(a)' ) 'set xlabel "'// trim(this%txtxlabel)//'"'
    end if
    if (this%hasylabel) then
        write ( file_unit, '(a)' ) 'set ylabel "'//trim(this%txtylabel)//'"'
    end if
    if (this%haszlabel) then
        write ( file_unit, '(a)' ) 'set zlabel "'//trim(this%txtzlabel)//'"'
    end if
    !    write ( file_unit, '(a)' ) 'set xrange '//trim(this%xrange)
    !    write ( file_unit, '(a)' ) 'set yrange '//trim(this%yrange)
    !    write ( file_unit, '(a)' ) 'set logscale '//trim(this%logscale)

    end subroutine processcmd
    !!}}}
!subroutine gnuplotscript{{{
    !..............................................................................
subroutine gnuplotscript(this, strscript)
    !..............................................................................
    ! write a gnuplot script in a file and then call gnuplot to execute the script
    class(gpf):: this
    character(len=*), intent(in):: strscript
    !local variables
    integer:: file_unit

    ! open the output file
    open ( unit = file_unit, file = this%txtfilename, status = 'replace',iostat = this%status )
    if (this%status /= 0 ) then
        this%msg= "gpf error: cannot open file for output"
        print*, this%msg
        return !an error has been occurred
    end if


    ! write gnuplot script in the file
    call writestring(strscript, file_unit)

    if (this%persist .and. this%displayplot) then
        close ( unit = file_unit )
        write(*,*)
        write(*,*) 'this is gnuplot interactive mode!'
        write(*,*)  'type q and press enter to exit'
        write(*,*)
    else
        write ( file_unit, '(a)' ) 'pause -1 "press a key to continue..."'
        write ( file_unit, '(a)' ) 'q'
        close ( unit = file_unit )
    end if
    !   now plot the results
    if (this%displayplot) then
        call write2gnuplot(this%txtfilename,this%persist)
    else
        this%displayplot=.true. !reset display plot to its value
        this%txtfilename=gnuplot_output_filename
    end if

    end subroutine gnuplotscript
!!}}}
!subroutine reset_to_defaults{{{
    !..............................................................................
subroutine reset_to_defaults(this)
    !..............................................................................
    !reset all params to their default values
    class(gpf):: this
    this%txtplottitle=""
    this%txtxlabel=""
    this%txtylabel=""
    this%txtoptions=""
    this%txtfilename=gnuplot_output_filename
    this%displayplot=.true.
    this%hasoptions= .false.
    this%hasplottitle= .false.
    this%hasxlabel= .false.
    this%hasylabel= .false.
    this%haszlabel= .false.

    this%persist= .false.
    end subroutine reset_to_defaults
!!}}}
!subroutine set_filename{{{
    !..............................................................................
subroutine set_filename(this,string)
    !..............................................................................
    !set a file name for plot command output
    !this file can be used later by gnuplot as an script to reproduce the plot
    class(gpf):: this
    character(len=*), intent(in) :: string
    this%txtfilename=trim(string)
    this%displayplot=.false.  ! set dispalyplot to false to only write plot commands into file
    end subroutine set_filename
!!}}}
!subroutine set_options{{{
  !..............................................................................
subroutine set_options(this,string)
    !..............................................................................
    !set the plot title
    class(gpf):: this
    character(len=*), intent(in) :: string
    integer, save :: strlength=0
    strlength=strlength+len_trim(string)
    if (strlength < len_options) then
        this%txtoptions=trim(this%txtoptions)//';'//trim(string)
    else
        print*, 'ogpf warning: the length of options exceeds than the set value :', len_options
        print*, 'options is truncated to ', len_options
        this%txtoptions=trim(this%txtoptions)//';'//trim(string)
    end if
    this%hasoptions=.true.
    end subroutine set_options
!!}}}
!subroutine set_plottitle{{{
    !..............................................................................
subroutine set_plottitle(this,string)
    !..............................................................................
    !set the plot title
    class(gpf):: this
    character(len=*), intent(in) :: string
    this%txtplottitle=trim(string)
    this%hasplottitle=.true.
    end subroutine set_plottitle
!!}}}
!subroutine set_xlabel{{{
    !..............................................................................
subroutine set_xlabel(this,string)
    !..............................................................................
    !set the xlabel
    class(gpf):: this
    character(len=*), intent(in) :: string
    this%txtxlabel=trim(string)
    this%hasxlabel=.true.
    end subroutine set_xlabel
!!}}}
!subroutine set_ylabel{{{
    !..............................................................................
subroutine set_ylabel(this,string)
    !..............................................................................
    !set the ylabel
    class(gpf):: this
    character(len=*), intent(in) :: string
    this%txtylabel=trim(string)
    this%hasylabel=.true.
    end subroutine set_ylabel
!!}}}
!subroutine set_zlabel{{{
    !..............................................................................
subroutine set_zlabel(this,string)
    !..............................................................................
    !set the z label
    class(gpf):: this
    character(len=*), intent(in) :: string
    this%txtzlabel=trim(string)
    this%haszlabel=.true.
    end subroutine set_zlabel
!!}}}
!subroutine set_persist{{{
    !..............................................................................
subroutine set_persist(this,string)
    !..............................................................................
    !set persist for gnuplot
    ! -on  keeps gnuplot in interactive mode
    ! -off closess gnuplot after drawing the plot
    class(gpf):: this
    character(len=*), intent(in) :: string
    select case(lcase(string))
    case ('on')
        this%persist=.true.
    case ('off')
        this%persist=.false.
    case default
        this%persist=.false.
    end select
    end subroutine set_persist
!!}}}
!function linestyle{{{
function linestyle(string)
    
             !this function accepts both abbreviated and full text line style and
             !return the correct linestyle value for gnuplot
             !if use pass wrong value, 'lines' is passed as the line style
            character(len=*), intent(in) :: string
            character(len=11 ):: linestyle
    
            select case( lcase(string) )
                case ('linespoints','lp')
                    linestyle='linespoints'
                case ('points','p')
                    linestyle='points'
                case ('lines','l')
                    linestyle='lines'
                case default  ! in the case of error use lines
                    linestyle='lines'
            end select
        end function linestyle
!!}}}
!subroutine errhandler{{{
    !           ***

    !..............................................................................
subroutine errhandler(msg)
    !..............................................................................
    character(len=*), intent(in)    :: msg
    write(6,*)    msg // ' failed'
    write(6,*)    'error in gpf '
    write(6,*)    'see gpf gui module....   '
    stop
    end subroutine errhandler
!!}}}
!subroutine  plot2d_matrix_vs_vector{{{
    !..............................................................................
subroutine  plot2d_matrix_vs_vector(this, xv,ymat, lspec)
    !..............................................................................
    !plot2d_matrix_vs_vector accepts a vector xv and a matrix ymat and plots columns of ymat against xv
    !lspec is an optional array array defines the line specification for each data series
    !if a single element array is sent for lspec then all series are plotted using the same
    !linespec

    implicit none
    class(gpf):: this
    ! input arrays
    real(8),  intent(in)    :: xv(:)
    real(8),  intent(in)    :: ymat(:,:)
    character(len=*),  intent(in), optional   :: lspec(:)
    !----------------------------------------------------------------------
    !       local variables
    integer:: nx
    integer:: ny
    integer:: file_unit
    integer:: ns
    integer:: number_of_plots
    integer::  i
    integer:: j
    integer:: ierr
    character(len=80), allocatable ::  pltstring(:)
    !
    !   check the input
    nx=size(xv)
    ny=size(ymat,dim=1)
    if (.not. checkdim(nx,ny)) then
        print*, 'gpf error: the length of arrays does not match'
        return
    end if

    ! open the output file
    open ( unit = file_unit, file = this%txtfilename, status = 'replace',iostat = this%status )
    if (this%status /= 0 ) then
        this%msg= "gpf error: cannot open file for output"
        print*, this%msg
        return !an error has been occurred
    end if

    ! write titles and other annotations
    call processcmd(this, file_unit)
    !   process legends and style
    number_of_plots=size(ymat,dim=2)
    allocate(pltstring(number_of_plots), stat=ierr)
    if (ierr /=0) then
        print*, 'allocation error'
        return
    end if

    if ( present(lspec) ) then
        call process_linespec(1,pltstring(1),lspec(1))
        ns=size(lspec)
        if (ns>1) then
            do i=2, number_of_plots
                call process_linespec(i,pltstring(i),lspec(i))
            end do
        else ! use the same setting for all plots
            do i=2, number_of_plots
                call process_linespec(i,pltstring(i),lspec(1)) !use the same lspec(1)
            end do
        end if
    else !no lspec is available
        pltstring(1)=' plot "-" notitle,'
        pltstring(2:number_of_plots-1)='"-" notitle,'
        pltstring(number_of_plots)='"-" notitle'
    end if

     ! write plot command and line styles and legend if any
    write ( file_unit, '(a)' ) ( trim(pltstring(i)) // ' \' , i=1, number_of_plots-1)
    write ( file_unit, '(a)' )   trim(pltstring(number_of_plots))



    !   write data into script file
    do j=1, number_of_plots
        do i = 1, nx
            write ( file_unit, * ) xv(i),ymat(i,j)
        end do
        write ( file_unit, '(a)' ) 'e'  !end of jth set of data
    end do

    if (this%persist .and. this%displayplot) then
        close ( unit = file_unit )
        write(*,*)
        write(*,*) 'this is gnuplot interactive mode!'
        write(*,*)  'type q and press enter to exit'
        write(*,*)
    else
        write ( file_unit, '(a)' ) 'pause -1 "press a key to continue..."'
        write ( file_unit, '(a)' ) 'q'
        close ( unit = file_unit )
    end if
    !   now plot the results
    if (this%displayplot) then
        call write2gnuplot(this%txtfilename,this%persist)
    else
        this%displayplot=.true. !reset display plot to its value
        this%txtfilename='ogpf_temp_script_file.plt'
    end if
    !release memory
    if (allocated(pltstring)) then
        deallocate(pltstring)
    end if
    !: end of plot2d_matrix_vs_vector
    end subroutine  plot2d_matrix_vs_vector
!!}}}
!subroutine process_linespec{{{    
subroutine process_linespec(order, lsstring, lspec)
    !accepts the line specification and interpret it into a format
    !to be sent to gnuplot

    integer, intent(in) :: order !1 for the first data series
    character(len=*), intent(out) :: lsstring
    character(len=*), intent(in), optional :: lspec

    select case(order)
    case(1)
        if ( present(lspec) ) then
            if (hastitle(lspec)) then
                lsstring='plot "-" '//trim(lspec)
            else
                lsstring='plot "-" notitle '//trim(lspec)
            end if
        else
            lsstring='plot "-" notitle '
        end if
    case  default !e.g. 2, 3, 4, ...
        if (present(lspec)) then
            if (hastitle(lspec)) then
                lsstring=', "-" '// trim(lspec)
            else
                lsstring=', "-" notitle '// trim(lspec)
            end if
        else
            lsstring=', "-" notitle '
        end if
    end select
    end subroutine process_linespec
!!}}}
!function hastitle{{{
    !..............................................................................
function hastitle(string)
    !..............................................................................
    character(len=*), intent(in) :: string
    logical:: hastitle
    integer:: idx1
    integer:: idx2

    idx1=index( lcase(string),'title')     !check if title is passed
    idx2=index(' ' // lcase(string),' t ') !check if the abbreviated title 't' is passed. extra space is added
                                           ! at the beginning of string to find starting 't'
    if (idx1 /=0 .or. idx2 /=0 ) then
        hastitle=.true.
    else
        hastitle=.false.
    end if

    end function hastitle
!!}}}
!function lcase{{{
    !!!!other utility

    !..............................................................................
function lcase(string)
    !..............................................................................
    character(len=*), intent(in) :: string
    character(len=len(string)):: lcase
    integer:: i
    integer:: n
    character(1):: chr

    do i=1, len(string)
        chr=string(i:i)
        n=ichar(chr)
        if (n >=65 .and. n <= 90) then
            lcase(i:i)=char(n+32)
        else
            lcase(i:i)=chr
        end if
    end do
    end function lcase
!!}}}
!function checkdim{{{
    !..............................................................................
function checkdim(nx,ny)
    !..............................................................................
    integer, intent(in):: nx
    integer, intent(in):: ny
    logical:: checkdim
    if (nx/=ny) then
        print*,  ' gpf - fatal error!'//char(13)// &
        ' the length of x, and y vectors must be equal.'
        checkdim=.false.
    else
        checkdim=.true.
    end if

    end function checkdim
!!}}}
!subroutine write_xydata{{{
    !..............................................................................
subroutine write_xydata(file_unit,ndata,x,y)
    !..............................................................................
    ! writes set of xy data into a file
    integer, intent(in) :: file_unit
    integer, intent(in) :: ndata
    real(8), intent(in) :: x(:)
    real(8), intent(in) :: y(:)

    integer:: i
    do i = 1, ndata
        write ( file_unit, * ) x(i), y(i)
    end do
    write ( file_unit, '(a)' ) 'e'  !end of set of data
    end subroutine write_xydata
!!}}}
!function linspace{{{
function linspace(a,b,n_elements)
    !
    !   returns a linearly spaced vector with n points in [a, b]
    !   if n is omitted, 100 points will be considered
    !
    real(8), intent(in) ::a
    real(8), intent(in) :: b
    integer, intent(in), optional ::  n_elements
    real(8), allocatable:: linspace(:)
    !   local vars
    real(8):: dx
    integer:: i
    integer:: n
    integer:: ierr
    if (present(n_elements)) then
        n=n_elements
    else
        n=100
    end if

    allocate(linspace(n), stat=ierr)
    if (ierr /= 0) then
        print*, "fatal error, allocation failed in linspace function"
        return
    end if
    dx=(b-a)/dble(n-1)
    linspace=[(i*dx+a, i=0,n-1)]
    end function linspace
!!}}}
!subroutine meshgrid{{{
subroutine meshgrid(x,y,xgv,ygv, ierr)
    !meshgrid generate mesh grid over a rectangular domain of [xmin xmax, ymin, max]
    !     xgv, ygv are vector in form of [start, stop, step] or [start, stop]
    !     when step has not been given, nx=m and ny=n are used to calculate steps
    !     and generate nx and ny data points respectively.
    !     now these values are set to 25.
    !     x and y are matrix each of size [ny by nx] contains the grid data.
    !     the coordinates of point (i,j) is [x(i,j), y(i,j)]
    !     """
    !     # example
    !     # call meshgrid(x, y, [0,3,1],[5,8,1])
    !     # x
    !     # [[0.0, 1.0, 2.0, 3.0],
    !     # [0.0, 1.0, 2.0, 3.0],
    !     # [0.0, 1.0, 2.0, 3.0],
    !     # [0.0, 1.0, 2.0, 3.0]]
    !     #
    !     #y
    !     #[[5.0, 5.0, 5.0, 5.0],
    !     # [6.0, 6.0, 6.0, 6.0],
    !     # [7.0, 7.0, 7.0, 7.0],
    !     # [8.0, 8.0, 8.0, 8.0]]




    ! arguments
    real(8), intent(out), allocatable :: x(:,:)
    real(8), intent(out), allocatable :: y(:,:)
    real(8), intent(in) :: xgv(:)
    real(8), intent(in),  optional  :: ygv(:)
    integer, intent(out), optional :: ierr !return the error code
    ! local variables
    integer:: i
    integer:: sv
    integer:: nx
    integer:: ny
    integer:: m
    integer:: n
    real(8):: dx
    real(8):: dy
    logical:: only_x_available

! initial setting
    m=25
    n=25
    only_x_available=.false.
    sv=0 !assume no error
    select case(size(xgv))
    case (2)
        nx=m
        dx=(xgv(2)-xgv(1))/dble(nx-1)
    case (3)
        dx=xgv(3)
        nx=int((xgv(2)-xgv(1))/dx)+1
    case default
        print*, "wrong x vector to meshgrid"
        sv=1
    end select

    if (present(ygv)) then
        select case(size(ygv))
        case (2)
            dy=ygv(2)-ygv(1)
            ny=n
        case (3)
            dy=ygv(3)
            ny=int((ygv(2)-ygv(1))/dy)+1
        case default
            print*, "wrong y vector to meshgrid"
            sv=1
        end select
    else
        only_x_available=.true.
        ny=nx
    end if

    allocate(x(ny,nx),y(ny,nx),stat=sv)
    if (sv /=0) then
        print*, "allocataion erro in meshgrid"
    end if

    x(1,:)=[(xgv(1)+dble(i-1)*dx, i=1, nx)]
    x(2:ny,:)=spread(x(1,:),dim=1,ncopies=ny-1)

    if (only_x_available) then
        y=transpose(x)
    else
        y(:,1)=[(ygv(1)+dble(i-1)*dy, i=1, ny)]
        y(:,2:nx)=spread(y(:,1),dim=2,ncopies=nx-1)
    end if

    if (present(ierr)) then
        ierr=sv
    end if

    end subroutine meshgrid
!!}}}
    end module ogpf
