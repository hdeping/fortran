!  glut fortran 90 program to draw line



program simple1
use opengl_glut
use opengl_gl
use opengl_glu
implicit none
interface
     subroutine display() bind(c)
     end subroutine display
end interface
integer :: i

call glutinit()
!call glutinitdisplaymode(glut_single+glut_rgb)
call glutinitwindowsize( 550, 550 );
call glutinitwindowposition( 300, 300 );
i =  glutcreatewindow( " simple1 " );

!-- init
call glclearcolor(0.,0.,0.,0.) !-- background color
call glmatrixmode(gl_projection)
call glloadidentity()
call glortho(0.0d0,1.0d0,0.0d0,1.0d0,-1.0d0,1.0d0)
    
call glutdisplayfunc(display)  !-- callback to draw 
call glutmainloop


contains
!subroutine display{{{
subroutine display() bind(c)

call glclear(gl_color_buffer_bit)
call glcolor3f(1.0,1.0,1.0)

call glbegin(gl_polygon);
call glvertex3f(0.25, 0.25, 0.0);
call glvertex3f(0.75, 0.25, 0.0);
call glvertex3f(0.75, 0.75, 0.0);
call glvertex3f(0.25, 0.75, 0.0);
call glend();

!--call glbegin(gl_lines);
!--    call glvertex3f(0.25, 0.25, 0.0)
!--    call glvertex3f(0.75, 0.75, 0.0)
!--call glend( );

call glflush()

end subroutine display

!}}}
end program simple1

