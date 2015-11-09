    !  glut fortran 90 program to draw line 
program simple3
    use opengl_glut
    use opengl_gl
    implicit none
    integer           :: i
    real              :: x
    
    call glutinit()
    call glutinitdisplaymode(glut_single+glut_rgb)
    call glutinitwindowsize( 250, 250 );
    call glutinitwindowposition( 200, 200 );
    i =  glutcreatewindow( " simple1 " );
    
    call glclearcolor(0.,0.,0.,0.) !-- background color
    call glmatrixmode(gl_projection)
    call glloadidentity()
    call glortho(0.0d0,1.0d0,0.0d0,1.0d0,-1.0d0,1.0d0)
    
    i = 4
    x = dble(i)*6D-2
    call glutdisplayfunc(display)  !-- callback to draw 
    call glutmainloop
    
    contains
subroutine display !bind(c)
    implicit none
    
    call glclear(gl_color_buffer_bit)
    call glcolor3f(1.0,1.0,1.0)
    call glbegin(gl_lines);
    call glvertex3f(0.25, 0.25, 0.0)
    call glvertex3f(0.75, x, 0.0)
    call glend();
    call glflush()
end subroutine display
end program simple3
