!  glut fortran 90 program to draw line
!  nasser m. abbasi  072312, gfortran 4.6.3

!------------------------------------
program simple2
use opengl_glut
use opengl_gl
use iso_c_binding
implicit none

integer :: i
real    :: g_spin= 0.0

call glutinit()
call glutinitdisplaymode(glut_double+glut_rgb)
call glutinitwindowsize( 350, 350 );
call glutinitwindowposition( 100, 100 );
i =  glutcreatewindow( "left mouse to rotate, middle to stop" );

!-- init
call glclearcolor(0.,0.,0.,0.) !-- background color
call glshademodel(gl_flat)
    
call glutdisplayfunc(display)  !-- callback to draw 
call glutreshapefunc(reshape) 
call glutmousefunc(mouse) 
call glutmainloop

contains
!subroutine display{{{
subroutine display() bind(c)
character(len=20) :: st
integer :: i
call glclear(gl_color_buffer_bit+gl_depth_buffer_bit)
call glpushmatrix()
call glrotatef(g_spin, 0.0 ,0.0 ,1.0)
call glcolor3f(1., 1.,1.)
call glrectf(-25.0,-25.0,25.0, 25.0)
call glpopmatrix()
call glrasterpos3f(0.,40.,0.0);

write( st, '( f12.8 )' ) g_spin
do i=1,len(trim(st))   
   call glutbitmapcharacter(glut_bitmap_9_by_15,ichar(st(i:i)));
end do

call glutswapbuffers()
end subroutine display
!}}}
!subroutine reshape{{{
subroutine reshape(w,h) bind(c)
integer (c_int), value :: w
integer (c_int), value :: h     

call glviewport(0, 0,  w,  h);
call glmatrixmode(gl_projection);
call glloadidentity();
call glortho(-50.0d0, 50.0d0, -50.0d0, 50.0d0, -1.0d0, 1.0d0);
call glmatrixmode(gl_modelview);
call glloadidentity();
end subroutine reshape
!}}}
!subroutine mouse{{{
subroutine mouse(button,state,x,y) !bind(c)
integer (c_int), value :: button
integer (c_int), value :: state
integer (c_int), value :: x
integer (c_int), value :: y    

select case(button) 
   case(glut_left_button) 
     if (state .eq. glut_down) then
        call glutidlefunc(spindisplay)
     end if
   case(glut_middle_button)
     if (state .eq. glut_down) then
        call glutidlefunc()
     end if
end select       
end subroutine mouse
!}}}
!subroutine spindisplay{{{
subroutine spindisplay() bind(c)
g_spin = g_spin + 0.5
if (g_spin>360.0) then
   g_spin = g_spin - 360.0
end if

call glutpostredisplay()       
end subroutine spindisplay
!}}}
end program simple2

