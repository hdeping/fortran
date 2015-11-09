!--------------------------------------------------------------------------

! This is a simple program to demonstrate the use of the view_modifier module
! It consists of a module with the callback functions and a main program.

module view_demo_callbacks
use opengl_gl
use opengl_glut
use view_modifier
private
public :: display

contains

subroutine display() bind(c)

! This gets called when the display needs to be redrawn

call reset_view

call glClear(ior(GL_COLOR_BUFFER_BIT,GL_DEPTH_BUFFER_BIT))
call glCallList(1_GLint)
call glutSwapBuffers

return
end subroutine display

end module view_demo_callbacks

