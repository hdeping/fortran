module function_plotter
use opengl_gl
use opengl_glut
use view_modifier, ONLY : reset_view
implicit none
private
public :: display, draw_func, make_menu

! symbolic constants

integer, parameter :: surfgrid_toggle = 1, &
                      surfsolid_toggle = 2, &
                      contour_toggle = 3, &
                      quit_selected = 4

integer, parameter :: set_nx = 1, &
                      set_ny = 2, &
                      set_ncontour = 3, &
                      set_contour_val = 4, &
                      set_xrange = 5, &
                      set_yrange = 6, &
                      reset_params = 7

integer, parameter :: black_contour = 1, &
                      rainbow_contour = 2

integer, parameter :: white_surface = 1, &
                      red_surface = 2, &
                      rainbow_surface = 3

! Default initial settings

integer, parameter :: init_ngridx = 40, &
                      init_ngridy = 40, &
                      init_num_contour = 20, &
                      init_contour_color = black_contour, &
                      init_surface_color = rainbow_surface

real(GLDOUBLE), parameter :: init_minx = 0.0_GLDOUBLE, &
                             init_maxx = 1.0_GLDOUBLE, &
                             init_miny = 0.0_GLDOUBLE, &
                             init_maxy = 1.0_GLDOUBLE

logical, parameter :: init_draw_surface_grid = .false., &
                      init_draw_surface_solid = .true., &
                      init_draw_contour = .true.

! Current settings

integer :: ngridx = init_ngridx, &
           ngridy = init_ngridy, &
           num_contour = init_num_contour, &
           contour_color = init_contour_color, &
           surface_color = init_surface_color

real(GLDOUBLE) :: minx = init_minx, &
                  maxx = init_maxx, &
                  miny = init_miny, &
                  maxy = init_maxy, &
                  minz = 0.0_GLDOUBLE, &
                  maxz = 0.0_GLDOUBLE

logical :: draw_surface_grid = init_draw_surface_grid, &
           draw_surface_solid = init_draw_surface_solid, &
           draw_contour = init_draw_contour, &
           contour_values_given = .false.

real(GLDOUBLE), allocatable :: actual_contours(:)

contains

subroutine display() bind(c)

! This gets called when the display needs to be redrawn

call reset_view

call glClear(ior(GL_COLOR_BUFFER_BIT,GL_DEPTH_BUFFER_BIT))
call glCallList(1)
call glutSwapBuffers

return
end subroutine display

subroutine draw_func
real(GLDOUBLE) :: gridx(0:ngridx),gridy(0:ngridy),zval(0:ngridy,2)
integer :: i,j,k,cont
real(GLDOUBLE) :: x1,x2,x3,xt,y1,y2,y3,yt,z1,z2,z3,zt
real(GLDOUBLE) :: frac,xcross1,xcross2,ycross1,ycross2
real(GLDOUBLE) :: contour_value(num_contour)
real(GLFLOAT) :: color(4), normal(3), &
                 red(4) = (/1.0,0.0,0.0,1.0/), &
                 black(4) = (/0.0,0.0,0.0,1.0/), &
                 white(4) = (/1.0,1.0,1.0,1.0/)
real(GLDOUBLE), external :: func_to_plot

! prepare to make a new display list

call reset_view
call glDeleteLists(1_gluint, 1_glsizei)
call glNewList(1_gluint, gl_compile_and_execute)

! set the grid points

gridx = (/ ((minx + i*(maxx-minx)/ngridx),i=0,ngridx) /)
gridy = (/ ((miny + i*(maxy-miny)/ngridy),i=0,ngridy) /)

! if this is the first call and either rainbow coloring of a solid surface
! or contours are being drawn, minz and maxz need to be set

if (minz == 0.0_GLDOUBLE .and. maxz == 0.0_GLDOUBLE .and. &
    ( (draw_surface_solid .and. surface_color == rainbow_surface) .or. &
      draw_contour ) ) then
   do i=0,ngridx
      do j=0,ngridy
         z1 = func_to_plot(gridx(i),gridy(j))
         minz = min(z1,minz)
         maxz = max(z1,maxz)
      end do
   end do
endif

! draw the solid surface

if (draw_surface_solid) then

   call glPolygonMode(gl_front_and_back, gl_fill)
   call glBegin(gl_triangles)

! set the color for a red or white surface.  For white, set the lighting
! such that there is uniform brightness.

   if (surface_color == red_surface) then
      call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,red)
   elseif (surface_color == white_surface) then
      call glDisable(gl_light0)
      call glLightModelfv(gl_light_model_ambient, (/1.,1.,1.,1./))
      call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,white)
   endif

! compute the function values for the first line of points

   do j=0,ngridy
      zval(j,2) = func_to_plot(gridx(0),gridy(j))
   end do

! for each x grid interval...

   do i=1,ngridx

! copy left side function values from the right side of the previous interval

      zval(:,1) = zval(:,2)

! compute the function values for the right side of the interval

      do j=0,ngridy
         zval(j,2) = func_to_plot(gridx(i),gridy(j))
      end do

      minz = min(minz,minval(zval))
      maxz = max(maxz,maxval(zval))

! for each y grid interval ...

      do j=1,ngridy

! add two triangles to the display list
! Before each triangle, set the normal.  Before each vertex, set the
! color if we're coloring by height

         normal = normcrossprod((/gridx(i-1),gridx(i),gridx(i)/), &
                                (/gridy(j-1),gridy(j-1),gridy(j)/), &
                                (/zval(j-1,1),zval(j-1,2),zval(j,2)/))
         call glNormal3fv(normal)
         if (surface_color == rainbow_surface) then
            call get_rainbow(zval(j-1,1),minz,maxz,color)
            call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,color)
         endif
         call glvertex3d(gridx(i-1),gridy(j-1),zval(j-1,1))
         if (surface_color == rainbow_surface) then
            call get_rainbow(zval(j-1,2),minz,maxz,color)
            call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,color)
         endif
         call glvertex3d(gridx(i  ),gridy(j-1),zval(j-1,2))
         if (surface_color == rainbow_surface) then
            call get_rainbow(zval(j,2),minz,maxz,color)
            call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,color)
         endif
         call glvertex3d(gridx(i  ),gridy(j  ),zval(j  ,2))
         normal = normcrossprod((/gridx(i-1),gridx(i-1),gridx(i)/), &
                                (/gridy(j),gridy(j-1),gridy(j)/), &
                                (/zval(j,1),zval(j-1,1),zval(j,2)/))
         call glNormal3fv(normal)
         if (surface_color == rainbow_surface) then
            call get_rainbow(zval(j,2),minz,maxz,color)
            call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,color)
         endif
         call glvertex3d(gridx(i  ),gridy(j  ),zval(j  ,2))
         if (surface_color == rainbow_surface) then
            call get_rainbow(zval(j-1,1),minz,maxz,color)
            call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,color)
         endif
         call glvertex3d(gridx(i-1),gridy(j-1),zval(j-1,1))
         if (surface_color == rainbow_surface) then
            call get_rainbow(zval(j,1),minz,maxz,color)
            call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,color)
         endif
         call glvertex3d(gridx(i-1),gridy(j  ),zval(j  ,1))

      end do
   end do

   call glEnd

! if the surface is white, reset the lighting conditions

   if (surface_color == white_surface) then
      call glEnable(gl_light0)
      call glLightModelfv(gl_light_model_ambient, (/.5,.5,.5,1./))
   endif
   

endif ! draw_surface_solid

! draw the surface grid

if (draw_surface_grid) then

   call glPolygonMode(gl_front_and_back, gl_line)
   call glBegin(gl_triangles)

! draw surface grid in black

   call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,black)

! compute the function values for the first line of points

   do j=0,ngridy
      zval(j,2) = func_to_plot(gridx(0),gridy(j))
   end do

! for each x grid interval...

   do i=1,ngridx

! copy left side function values from the right side of the previous interval

      zval(:,1) = zval(:,2)

! compute the function values for the right side of the interval

      do j=0,ngridy
         zval(j,2) = func_to_plot(gridx(i),gridy(j))
      end do

      minz = min(minz,minval(zval))
      maxz = max(maxz,maxval(zval))

! for each y grid interval ...

      do j=1,ngridy

! add two triangles to the display list

         call glvertex3d(gridx(i-1),gridy(j-1),zval(j-1,1))
         call glvertex3d(gridx(i  ),gridy(j-1),zval(j-1,2))
         call glvertex3d(gridx(i  ),gridy(j  ),zval(j  ,2))
         call glvertex3d(gridx(i  ),gridy(j  ),zval(j  ,2))
         call glvertex3d(gridx(i-1),gridy(j-1),zval(j-1,1))
         call glvertex3d(gridx(i-1),gridy(j  ),zval(j  ,1))

      end do
   end do

   call glEnd

endif ! draw_surface_grid

! draw the contour plot

if (draw_contour) then

   call glPolygonMode(gl_front_and_back, gl_line)
   call glBegin(gl_lines)
   call glNormal3fv((/0.0_glfloat, 0.0_glfloat, 1.0_glfloat/))

! draw the domain in black, also sets color to black for black_contour

   call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,black)
   call glVertex3d(minx,miny,0.0_GLDOUBLE)
   call glVertex3d(maxx,miny,0.0_GLDOUBLE)
   call glVertex3d(maxx,miny,0.0_GLDOUBLE)
   call glVertex3d(maxx,maxy,0.0_GLDOUBLE)
   call glVertex3d(maxx,maxy,0.0_GLDOUBLE)
   call glVertex3d(minx,maxy,0.0_GLDOUBLE)
   call glVertex3d(minx,maxy,0.0_GLDOUBLE)
   call glVertex3d(minx,miny,0.0_GLDOUBLE)

! set the contour values

   if (contour_values_given) then
      contour_value = actual_contours
   else
      do i=1,num_contour
         contour_value(i) = minz+(maxz-minz)*(i-1)/real(num_contour-1,GLDOUBLE)
      end do
   endif

! compute the function values for the first line of points

   do j=0,ngridy
      zval(j,2) = func_to_plot(gridx(0),gridy(j))
   end do

! for each x grid interval...

   do i=1,ngridx

! copy left side function values from the right side of the previous interval

      zval(:,1) = zval(:,2)

! compute the function values for the right side of the interval

      do j=0,ngridy
         zval(j,2) = func_to_plot(gridx(i),gridy(j))
      end do

      minz = min(minz,minval(zval))
      maxz = max(maxz,maxval(zval))

! for each y grid interval ...

      do j=1,ngridy

! for two triangles

         do k=1,2

! set the vertices

            if (k==1) then
               x1 = gridx(i-1); y1 = gridy(j-1); z1 = zval(j-1,1)
               x2 = gridx(i  ); y2 = gridy(j-1); z2 = zval(j-1,2)
               x3 = gridx(i  ); y3 = gridy(j  ); z3 = zval(j  ,2)
            else
               x1 = gridx(i-1); y1 = gridy(j-1); z1 = zval(j-1,1)
               x2 = gridx(i-1); y2 = gridy(j  ); z2 = zval(j  ,1)
               x3 = gridx(i  ); y3 = gridy(j  ); z3 = zval(j  ,2)
            endif

! order the vertices by z value

            xt = x1; yt = y1; zt = z1
            if (z2 < z1) then
               xt = x1; yt = y1; zt = z1
               if (z3 < z1) then
                  if (z2 < z3) then
                     x1 = x2; y1 = y2; z1 = z2
                     x2 = x3; y2 = y3; z2 = z3
                  else
                     x1 = x3; y1 = y3; z1 = z3
                  endif
                  x3 = xt; y3 = yt; z3 = zt
               else
                  x1 = x2; y1 = y2; z1 = z2
                  x2 = xt; y2 = yt; z2 = zt
               endif
            elseif (z3 < z1) then
               x1 = x3; y1 = y3; z1 = z3
               x3 = x2; y3 = y2; z3 = z2
               x2 = xt; y2 = yt; z2 = zt
            elseif (z3 < z2) then
               xt = x2; yt = y2; zt = z2
               x2 = x3; y2 = y3; z2 = z3
               x3 = xt; y3 = yt; z3 = zt
            endif

! if z1==z3, the triangle is horizontal and no contours pass through it

            if (z1==z3) cycle

! for each contour value

            do cont = 1,num_contour

! see if it passes through this triangle

               if (contour_value(cont) < z1) cycle
               if (contour_value(cont) > z3) exit

! set the color for contours colored by solution

               if (contour_color == rainbow_contour) then
                  call get_rainbow(contour_value(cont),minz,maxz,color)
                  call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse, &
                                    color)
               endif

! see where it crosses the 1-3 edge

               frac = (contour_value(cont)-z1)/(z3-z1)
               xcross1 = (1.0_GLDOUBLE - frac)*x1 + frac*x3
               ycross1 = (1.0_GLDOUBLE - frac)*y1 + frac*y3

! see where it crosses one of the other edges

               if (contour_value(cont) == z2) then
                  xcross2 = x2
                  ycross2 = y2
               elseif (contour_value(cont) < z2) then
                  frac = (contour_value(cont)-z1)/(z2-z1)
                  xcross2 = (1.0_GLDOUBLE - frac)*x1 + frac*x2
                  ycross2 = (1.0_GLDOUBLE - frac)*y1 + frac*y2
               else
                  frac = (contour_value(cont)-z2)/(z3-z2)
                  xcross2 = (1.0_GLDOUBLE - frac)*x2 + frac*x3
                  ycross2 = (1.0_GLDOUBLE - frac)*y2 + frac*y3
               endif

! add the line segment to the display list

               call glVertex3d(xcross1,ycross1,0.0_GLDOUBLE)
               call glVertex3d(xcross2,ycross2,0.0_GLDOUBLE)

            end do
         end do
      end do
   end do

   call glEnd

endif ! draw_contour

! finish off display list

call glEndList
call glutPostRedisplay

end subroutine draw_func

subroutine get_rainbow(val,minval,maxval,c)
real(GLDOUBLE), intent(in) :: val,maxval,minval
real(GLFLOAT), intent(out) :: c(4)

real(GLFLOAT) :: f

if (maxval > minval) then
   f = (val-minval)/(maxval-minval)
else ! probably maxval==minval
   f = 0.5_glfloat
endif

if (f < .25) then
   c(1) = 0.0_glfloat
   c(2) = 4.0_glfloat * f
   c(3) = 1.0_glfloat
   c(4) = 1.0_glfloat
elseif (f < .5) then
   c(1) = 0.0_glfloat
   c(2) = 1.0_glfloat
   c(3) = 2.0_glfloat - 4.0_glfloat*f
   c(4) = 1.0_glfloat
elseif (f < .75) then
   c(1) = 4.0_glfloat * f - 2.0_glfloat
   c(2) = 1.0_glfloat
   c(3) = 0.0_glfloat
   c(4) = 1.0_glfloat
else
   c(1) = 1.0_glfloat
   c(2) = 4.0_glfloat - 4.0_glfloat*f
   c(3) = 0.0_glfloat
   c(4) = 1.0_glfloat
endif

end subroutine get_rainbow

function normcrossprod(x,y,z)
real(glfloat), dimension(3) :: normcrossprod
real(gldouble), dimension(3), intent(in) :: x,y,z
real(glfloat) :: t1(3),t2(3),norm
t1(1) = x(2) - x(1)
t1(2) = y(2) - y(1)
t1(3) = z(2) - z(1)
t2(1) = x(3) - x(1)
t2(2) = y(3) - y(1)
t2(3) = z(3) - z(1)
normcrossprod(1) = t1(2)*t2(3) - t1(3)*t2(2)
normcrossprod(2) = t1(3)*t2(1) - t1(1)*t2(3)
normcrossprod(3) = t1(1)*t2(2) - t1(2)*t2(1)
norm = sqrt(dot_product(normcrossprod,normcrossprod))
if (norm /= 0._glfloat) normcrossprod = normcrossprod/norm
end function normcrossprod

subroutine menu_handler_fp(selection) bind(c,name="")
integer(kind=glcint), value :: selection

select case (selection)

case (surfgrid_toggle)
   draw_surface_grid = .not. draw_surface_grid
   call draw_func

case (surfsolid_toggle)
   draw_surface_solid = .not. draw_surface_solid
   call draw_func

case (contour_toggle)
   draw_contour = .not. draw_contour
   call draw_func

case (quit_selected)
   stop

end select

return
end subroutine menu_handler_fp

subroutine param_handler(selection) bind(c,name="")
integer(kind=glcint), value :: selection

select case (selection)

case (set_nx)
   print *,"Enter number of x subintervals:"
   read *, ngridx
   call draw_func

case (set_ny)
   print *,"Enter number of y subintervals:"
   read *, ngridy
   call draw_func

case (set_ncontour)
   print *,"Enter number of contour lines:"
   read *, num_contour
   contour_values_given = .false.
   call draw_func

case (set_contour_val)
   print *,"enter number of contours:"
   read *, num_contour
   if (allocated(actual_contours)) deallocate(actual_contours)
   allocate(actual_contours(num_contour))
   print *,"enter ",num_contour," contour values:"
   read *,actual_contours
   contour_values_given = .true.
   call draw_func

case (set_xrange)
   print *,"Enter minimum and maximum x values:"
   read *,minx,maxx
   call draw_func

case (set_yrange)
   print *,"Enter minimum and maximum y values:"
   read *,miny,maxy
   call draw_func

case (reset_params)
   ngridx = init_ngridx
   ngridy = init_ngridy
   num_contour = init_num_contour
   contour_color = init_contour_color
   surface_color = init_surface_color
   minx = init_minx
   maxx = init_maxx
   miny = init_miny
   maxy = init_maxy
   draw_surface_grid = init_draw_surface_grid
   draw_surface_solid = init_draw_surface_solid
   draw_contour = init_draw_contour
   call draw_func

end select

end subroutine param_handler

subroutine contour_color_handler(selection) bind(c,name="")
integer(kind=glcint), value :: selection

contour_color = selection
call draw_func

end subroutine contour_color_handler

subroutine surface_color_handler(selection) bind(c,name="")
integer(kind=glcint), value :: selection

surface_color = selection
call draw_func

end subroutine surface_color_handler

subroutine make_menu(submenuid)
integer, intent(in) :: submenuid
integer :: menuid, param_id, contour_color_menu, surface_color_menu

contour_color_menu = glutCreateMenu(contour_color_handler)
call glutAddMenuEntry(CString("black"),black_contour)
call glutAddMenuEntry(CString("contour value"),rainbow_contour)

surface_color_menu = glutCreateMenu(surface_color_handler)
call glutAddMenuEntry(CString("red"),red_surface)
call glutAddMenuEntry(CString("white"),white_surface)
call glutAddMenuEntry(CString("surface height"),rainbow_surface)

param_id = glutCreateMenu(param_handler)
call glutAddMenuEntry(CString("number of x grid intervals"),set_nx)
call glutAddMenuEntry(CString("number of y grid intervals"),set_ny)
call glutAddMenuEntry(CString("number of uniform contour lines"),set_ncontour)
call glutAddMenuEntry(CString("contour values"),set_contour_val)
call glutAddMenuEntry(CString("x range"),set_xrange)
call glutAddMenuEntry(CString("y range"),set_yrange)
call glutAddSubMenu(CString("contour color"),contour_color_menu)
call glutAddSubMenu(CString("surface color"),surface_color_menu)
call glutAddMenuEntry(CString("reset to initial parameters"),reset_params)

menuid = glutCreateMenu(menu_handler_fp)
call glutAddSubMenu(CString("View Modifier"),submenuid)
call glutAddMenuEntry(CString("toggle surface grid"),surfgrid_toggle)
call glutAddMenuEntry(CString("toggle solid surface"),surfsolid_toggle)
call glutAddMenuEntry(CString("toggle contour"),contour_toggle)
call glutAddSubMenu(CString("plotting parameters"),param_id)
call glutAddMenuEntry(CString("quit"),quit_selected)
call glutAttachMenu(GLUT_RIGHT_BUTTON)
end subroutine make_menu

end module function_plotter
