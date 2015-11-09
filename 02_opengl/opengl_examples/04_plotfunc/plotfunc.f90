! This program plots a function of two variables.  The function, called
! func_to_plot, is an external procedure at the end of the file.  
! This begins with the same module used in the modview example, followed by
! another module for plotting the function, called function_plotter.
! You might want to change default initial settings in modules modview and
! function_plotter.

! William F. Mitchell
! william.mitchell@nist.gov
! Mathematical and Computational Sciences Division
! National Institute of Standards and Technology
! August, 1999

!---------------------------------------------------------------------------


!---------------------------------------------------------------------------

program plotfunc

use opengl_gl
use opengl_glut
use view_modifier
use function_plotter
implicit none

integer :: winid, menuid, submenuid

! Initializations

call glutInit
call glutInitDisplayMode(ior(GLUT_DOUBLE,ior(GLUT_RGB,GLUT_DEPTH)))
call glutInitWindowSize(500_glcint,500_glcint)

! Create a window

winid = glutCreateWindow("Function plotter")

! initialize view_modifier, receiving the id for it's submenu

submenuid = view_modifier_init()

! create the menu

call make_menu(submenuid)

! Set the display callback

call glutDisplayFunc(display)

! set the lighting conditions

call glClearColor(0.9_glclampf, 0.9_glclampf, 0.9_glclampf, 1.0_glclampf)
call glLightfv(gl_light0, gl_diffuse, (/1.,1.,1.,1./))
call glLightfv(gl_light0, gl_position, (/1.5,-.5,2.,0.0/))
call glEnable(gl_lighting)
call glEnable(gl_light0)
call glLightModelfv(gl_light_model_ambient, (/.5,.5,.5,1./))
call glDepthFunc(gl_lequal)
call glEnable(gl_depth_test)

! Create the image

call draw_func

! Let glut take over

call glutMainLoop

end program plotfunc

!---------------------------------------------------------------------------

! The function to plot

function func_to_plot(x,y)
use opengl_gl
real(GLDOUBLE) :: func_to_plot
real(GLDOUBLE), intent(in) :: x,y

func_to_plot = 0.5_GLDOUBLE * ( &
       cos(0.3_GLDOUBLE*sqrt(80.0_GLDOUBLE*x)-16.0_GLDOUBLE*y/3.0_GLDOUBLE)* &
       cos(16.0_GLDOUBLE*x/3.0_GLDOUBLE) + x-y )

end function func_to_plot
