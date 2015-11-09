    module view_modifier
!     Read Me{{{
    ! this module provides facilities to modify the view in an opengl window.
    ! the mouse buttons and keyboard arrow keys can be used to zoom, pan,
    ! rotate and change the scale.  a menu or submenu can be used to select which
    ! buttons perform which function and to reset the view to the initial settings.
    ! this is limited to one window.
    
    ! william f. mitchell
    ! william.mitchell@nist.gov
    ! mathematical and computational sciences division
    ! national institute of standards and technology
    ! april, 1998
    
    ! to use this module:
    !
    ! 1) put a use view_modifier statement in any program unit that calls a
    !    procedure in this module
    !
    ! 2) set the initial operation assignments, view and scale below the
    !    "initial configuration" comment below
    !
    ! 3) call view_modifier_init after glutcreatewindow
    !    this is a function that returns integer(kind=glcint) menuid.  the menuid
    !    is the id returned by glutcreatemenu.  you can either use the view_modifier
    !    menu as your menu by calling glutattachmenu immediately after
    !    view_modifier_init, as in
    !       menuid = view_modifier_init()
    !       call glutattachmenu(glut_right_button)
    !    or by using the menuid to attach a submenu to your own menu, as in
    !       call glutaddsubmenu(cstring("view modifier"),menuid)
    !    note that string arguments have to be converted to c format
    !    (null-terminated character array). the cstring function does the
    !    conversion.
    
    !
    ! 4) in any callback functions that update the display, put
    !       call reset_view
    !    as the first executable statement
    !
    ! note that view_modifier_init sets the callback functions for glutmousefunc,
    ! glutmotionfunc and glutspecialfunc, so don't call these yourself
    !
    ! the menu allows you to select what operation is attached to the left and
    ! middle mouse buttons and arrow keys, reset to the initial view, and quit.
    ! the right mouse button should be used for the menu.
!}}}
    use opengl_gl
    use opengl_glu
    use opengl_glut

    implicit none
!parameters, initial configuration{{{
    private
    public :: view_modifier_init, reset_view
    !  private :: zoom, pan, rotate, scalex, scaley, scalez, reset, above, quit, &
    !             pi, &
    !             left_button_func, middle_button_func, arrow_key_func, &
    !             init_lookat, init_lookfrom, &
    !             init_xscale_factor, init_yscale_factor, init_zscale_factor, &
    !             angle, shift, xscale_factor, yscale_factor, zscale_factor, &
    !             moving_left, moving_middle, begin_left, begin_middle, &
    !             cart2sphere, sphere2cart, cart3d_plus_cart3d, cart3d_minus_cart3d, &
    !             reset_to_init, mouse, motion, arrows, &
    !             menu_handler, set_left_button, set_middle_button, set_arrow_keys
    
    integer(kind=glcint), parameter :: zoom = 1, pan = 2, rotate = 3, scalex = 4, &
                          scaley = 5, scalez = 6
    integer(kind=glcint), parameter :: reset = 10, above = 11, quit = 12
    real(kind=gldouble), parameter :: pi = 3.141592653589793_gldouble
    
    type, private :: cart2d ! 2d cartesian coordinates
       real(kind=gldouble) :: x, y
    end type cart2d
    
    type, private :: cart3d ! 3d cartesian coordinates
       real(kind=gldouble) :: x, y, z
    end type cart3d
    
    type, private :: sphere3d ! 3d spherical coordinates
       real(kind=gldouble) :: theta, phi, rho
    end type sphere3d
    
    type(cart2d), save :: angle
    type(cart3d), save :: shift
    real(kind=gldouble), save :: xscale_factor, yscale_factor, zscale_factor
    logical, save :: moving_left, moving_middle
    type(cart2d), save :: begin_left, begin_middle
    
    interface operator(+)
       module procedure cart3d_plus_cart3d
    end interface
    interface operator(-)
       module procedure cart3d_minus_cart3d
    end interface
    
    ! ------- initial configuration -------
    
    ! set the initial operation performed by each button and the arrow keys.
    ! the operations are zoom, pan, rotate, scalex, scaley, and scalez
    
    integer, save ::   left_button_func = rotate, &
                       middle_button_func = zoom, &
                       arrow_key_func = pan
    
    ! set the initial view as the point you are looking at, the point you are
    ! looking from, and the scale factors
    
    type(cart3d), parameter :: &
         init_lookat = cart3d(0.5_gldouble, 0.5_gldouble, 0.0_gldouble), &
       init_lookfrom = cart3d(5.0_gldouble, 10.0_gldouble, 2.5_gldouble)
    
    real(kind=gldouble), parameter :: &
       init_xscale_factor = 1.0_gldouble, &
       init_yscale_factor = 1.0_gldouble, &
       init_zscale_factor = 1.0_gldouble
    ! -------- end of initial configuration ------
!}}}
    contains
!subroutine reset_to_init{{{
subroutine reset_to_init
    ! this resets the view to the initial configuration
    type(sphere3d) :: slookfrom
    
    slookfrom     = cart2sphere(init_lookfrom-init_lookat)
    angle%x       = -180.0_gldouble*slookfrom%theta/pi - 90.0_gldouble
    angle%y       = -180.0_gldouble*slookfrom%phi/pi
    shift%x       = 0.0_gldouble
    shift%y       = 0.0_gldouble
    shift%z       = -slookfrom%rho
    xscale_factor = init_xscale_factor
    yscale_factor = init_yscale_factor
    zscale_factor = init_zscale_factor
    
    call glutpostredisplay
    
    return
end subroutine reset_to_init
!}}}
!subroutine view_from_above{{{
subroutine view_from_above()
    ! this sets the view to be from straight above
    
    type(sphere3d) :: slookfrom
    
    slookfrom = cart2sphere(cart3d(0.0,0.0,1.0))
    angle%x = -180.0_gldouble*slookfrom%theta/pi
    angle%y = -180.0_gldouble*slookfrom%phi/pi
    
    call glutpostredisplay
    
    return
end subroutine view_from_above
!}}}
!subroutine reset_view{{{
subroutine reset_view()
    
    ! this routine resets the view to the current orientation and scale
    
    call glmatrixmode(gl_modelview)
    call glpopmatrix
    call glpushmatrix
    call gltranslated(shift%x, shift%y, shift%z)
    call glrotated(angle%x, 0.0_gldouble, 0.0_gldouble, 1.0_gldouble)
    call glrotated(angle%y, cos(pi*angle%x/180.0_gldouble), &
                   -sin(pi*angle%x/180.0_gldouble), 0.0_gldouble)
    call gltranslated(-init_lookat%x, -init_lookat%y, -init_lookat%z)
    call glscaled(xscale_factor,yscale_factor,zscale_factor)
    
    return
    end subroutine reset_view
!}}}
!subroutine mouse{{{
subroutine mouse(button, state, x, y) bind(c,name="")
    integer(kind=glcint), value :: button, state, x, y
    
    ! this gets called when a mouse button changes
     
      if (button == glut_left_button .and. state == glut_down) then
        moving_left = .true.
        begin_left = cart2d(x,y)
      endif
      if (button == glut_left_button .and. state == glut_up) then
        moving_left = .false.
      endif
      if (button == glut_middle_button .and. state == glut_down) then
        moving_middle = .true.
        begin_middle = cart2d(x,y)
      endif
      if (button == glut_middle_button .and. state == glut_up) then
        moving_middle = .false.
      endif
end subroutine mouse
!}}}
!subroutine motion{{{
subroutine motion(x, y) bind(c,name="")
    integer(kind=glcint), value :: x, y
    ! this gets called when the mouse moves
    
    integer :: button_function
    type(cart2d) :: begin
    real(kind=gldouble) :: factor
    
    ! determine and apply the button function
    
    if (moving_left) then
       button_function = left_button_func
       begin = begin_left
    else if(moving_middle) then
       button_function = middle_button_func
       begin = begin_middle
    endif
    
    select case(button_function)
    case (zoom)
       if (y < begin%y) then
          factor = 1.0_gldouble/(1.0_gldouble + .002_gldouble*(begin%y-y))
       else if (y > begin%y) then
          factor = 1.0_gldouble + .002_gldouble*(y-begin%y)
       else
          factor = 1.0_gldouble
       endif
       shift%z = factor*shift%z
    case (pan)
       shift%x = shift%x + .01*(x - begin%x)
       shift%y = shift%y - .01*(y - begin%y)
    case (rotate)
       angle%x = angle%x + (x - begin%x)
       angle%y = angle%y + (y - begin%y)
    case (scalex)
       if (y < begin%y) then
          factor = 1.0_gldouble + .002_gldouble*(begin%y-y)
       else if (y > begin%y) then
          factor = 1.0_gldouble/(1.0_gldouble + .002_gldouble*(y-begin%y))
       else
          factor = 1.0_gldouble
       endif
       xscale_factor = xscale_factor * factor
    case (scaley)
       if (y < begin%y) then
          factor = 1.0_gldouble + .002_gldouble*(begin%y-y)
       else if (y > begin%y) then
          factor = 1.0_gldouble/(1.0_gldouble + .002_gldouble*(y-begin%y))
       else
          factor = 1.0_gldouble
       endif
       yscale_factor = yscale_factor * factor
    case (scalez)
       if (y < begin%y) then
          factor = 1.0_gldouble + .002_gldouble*(begin%y-y)
       else if (y > begin%y) then
          factor = 1.0_gldouble/(1.0_gldouble + .002_gldouble*(y-begin%y))
       else
          factor = 1.0_gldouble
       endif
       zscale_factor = zscale_factor * factor
    end select
    
    ! update private variables and redisplay
    
    if (moving_left) then
       begin_left = cart2d(x,y)
    else if(moving_middle) then
       begin_middle = cart2d(x,y)
    endif
    
    if (moving_left .or. moving_middle) then
       call glutpostredisplay
    endif
    
    return
end subroutine motion
!}}}
!subroutine arrows{{{
subroutine arrows(key, x, y) bind(c,name="")
    integer(glcint), value :: key, x, y
    
    ! this routine handles the arrow key operations
    
    real(kind=gldouble) :: factor
    
    select case(arrow_key_func)
    case(zoom)
       select case(key)
       case(glut_key_down)
          factor = 1.0_gldouble + .02_gldouble
       case(glut_key_up)
          factor = 1.0_gldouble/(1.0_gldouble + .02_gldouble)
       case default
          factor = 1.0_gldouble
       end select
       shift%z = factor*shift%z
    case(pan)
       select case(key)
       case(glut_key_left)
          shift%x = shift%x - .02
       case(glut_key_right)
          shift%x = shift%x + .02
       case(glut_key_down)
          shift%y = shift%y - .02
       case(glut_key_up)
          shift%y = shift%y + .02
       end select
    case(rotate)
       select case(key)
       case(glut_key_left)
          angle%x = angle%x - 1.0_gldouble
       case(glut_key_right)
          angle%x = angle%x + 1.0_gldouble
       case(glut_key_down)
          angle%y = angle%y + 1.0_gldouble
       case(glut_key_up)
          angle%y = angle%y - 1.0_gldouble
       end select
    case(scalex)
       select case(key)
       case(glut_key_down)
          factor = 1.0_gldouble/(1.0_gldouble + .02_gldouble)
       case(glut_key_up)
          factor = 1.0_gldouble + .02_gldouble
       case default
          factor = 1.0_gldouble
       end select
       xscale_factor = xscale_factor * factor
    case(scaley)
       select case(key)
       case(glut_key_down)
          factor = 1.0_gldouble/(1.0_gldouble + .02_gldouble)
       case(glut_key_up)
          factor = 1.0_gldouble + .02_gldouble
       case default
          factor = 1.0_gldouble
       end select
       yscale_factor = yscale_factor * factor
    case(scalez)
       select case(key)
       case(glut_key_down)
          factor = 1.0_gldouble/(1.0_gldouble + .02_gldouble)
       case(glut_key_up)
          factor = 1.0_gldouble + .02_gldouble
       case default
          factor = 1.0_gldouble
       end select
       zscale_factor = zscale_factor * factor
    
    end select
       
    call glutpostredisplay
    return
end subroutine arrows
!}}}
!subroutine menu_handler{{{
subroutine menu_handler(value) bind(c,name="")
    integer(kind=glcint), value :: value
    
    ! this routine handles the first level entries in the menu
    
    select case(value)
    
    case(reset)
       call reset_to_init
    case(above)
       call view_from_above
    case(quit)
       stop
    
    end select
    
    return
end subroutine menu_handler
!}}}
!subroutine set_left_button{{{
subroutine set_left_button(value) bind(c,name="")
    integer(kind=glcint), value :: value
    ! this routine sets the function of the left button as given by menu selection
    
    left_button_func = value
    
    return
    end subroutine set_left_button
    
    !          -----------------
    subroutine set_middle_button(value)  bind(c,name="")
    !          -----------------
    integer(kind=glcint), value :: value
    
    ! this routine sets the function of the middle button as given by menu selection
    
    middle_button_func = value
    
    return
end subroutine set_middle_button
!}}}
!subroutine set_arrow_keys{{{
subroutine set_arrow_keys(value)  bind(c,name="")

    integer(kind=glcint), value :: value
    
    ! this routine sets the function of the arrow keys as given by menu selection
    
    arrow_key_func = value
    
    return
end subroutine set_arrow_keys
!}}}
!function view_modifier_init{{{
function view_modifier_init() result(menuid)

    integer(kind=glcint) :: menuid
    
    ! this initializes the view modifier variables and sets initial view.
    ! it should be called immediately after glutcreatewindow
    
    integer(kind=glcint) :: button_left, button_middle, arrow_keys
    
    ! set the callback functions
    
    call glutmousefunc(mouse)
    call glutmotionfunc(motion)
    call glutspecialfunc(arrows)
    
    ! create the menu
    
    button_left = glutcreatemenu(set_left_button)
    call glutaddmenuentry(cstring("rotate"),rotate)
    call glutaddmenuentry(cstring("zoom"),zoom)
    call glutaddmenuentry(cstring("pan"),pan)
    call glutaddmenuentry(cstring("scale x"),scalex)
    call glutaddmenuentry(cstring("scale y"),scaley)
    call glutaddmenuentry(cstring("scale z"), scalez)
    button_middle = glutcreatemenu(set_middle_button)
    call glutaddmenuentry(cstring("rotate"),rotate)
    call glutaddmenuentry(cstring("zoom"),zoom)
    call glutaddmenuentry(cstring("pan"),pan)
    call glutaddmenuentry(cstring("scale x"),scalex)
    call glutaddmenuentry(cstring("scale y"),scaley)
    call glutaddmenuentry(cstring("scale z"), scalez)
    arrow_keys = glutcreatemenu(set_arrow_keys)
    call glutaddmenuentry(cstring("rotate"),rotate)
    call glutaddmenuentry(cstring("zoom"),zoom)
    call glutaddmenuentry(cstring("pan"),pan)
    call glutaddmenuentry(cstring("scale x"),scalex)
    call glutaddmenuentry(cstring("scale y"),scaley)
    call glutaddmenuentry(cstring("scale z"), scalez)
    menuid = glutcreatemenu(menu_handler)
    call glutaddsubmenu(cstring("left mouse button"),button_left)
    call glutaddsubmenu(cstring("middle mouse button"),button_middle)
    call glutaddsubmenu(cstring("arrow keys"),arrow_keys)
    call glutaddmenuentry(cstring("reset to initial view"),reset)
    call glutaddmenuentry(cstring("view from above"),above)
    call glutaddmenuentry(cstring("quit"),quit)
    
    ! set the perspective
    
    call glmatrixmode(gl_projection)
    call gluperspective(10.0_gldouble, 1.0_gldouble, 0.1_gldouble, 200.0_gldouble)
    
    ! set the initial view
    
    call glmatrixmode(gl_modelview)
    call glpushmatrix
    call reset_to_init
    
    return
end function view_modifier_init
!}}}
!function sphere2cart{{{
function sphere2cart(spoint) result(cpoint)

    type(sphere3d), intent(in) :: spoint
    type(cart3d) :: cpoint
    
    ! this converts a 3d point from spherical to cartesean coordinates
    
    real(kind=gldouble) :: t,p,r
    
    t=spoint%theta
    p=spoint%phi
    r=spoint%rho
    
    cpoint%x = r*cos(t)*sin(p)
    cpoint%y = r*sin(t)*sin(p)
    cpoint%z = r*cos(p)
    
    return
end function sphere2cart
!}}}
!function cart2sphere{{{
function cart2sphere(cpoint) result(spoint)

    type(cart3d), intent(in) :: cpoint
    type(sphere3d) :: spoint
    
    ! this converts a 3d point from cartesean to spherical coordinates
    
    real(kind=gldouble) :: x,y,z
    
    x=cpoint%x
    y=cpoint%y
    z=cpoint%z
    
    spoint%rho = sqrt(x*x+y*y+z*z)
    if (x==0.0_gldouble .and. y==0.0_gldouble) then
       spoint%theta = 0.0_gldouble
    else
       spoint%theta = atan2(y,x)
    endif
    if (spoint%rho == 0.0_gldouble) then
       spoint%phi = 0.0_gldouble
    else
       spoint%phi = acos(z/spoint%rho)
    endif
    
    return
end function cart2sphere
!}}}
!function cart3d_plus_cart3d{{{
function cart3d_plus_cart3d(cart1,cart2) result(cart3)

    type(cart3d), intent(in) :: cart1, cart2
    type(cart3d) :: cart3
    
    ! compute the sum of two 3d cartesean points
    
    cart3%x = cart1%x + cart2%x
    cart3%y = cart1%y + cart2%y
    cart3%z = cart1%z + cart2%z
    
    return
end function cart3d_plus_cart3d
!}}}
!function cart3d_minus_cart3d{{{
function cart3d_minus_cart3d(cart1,cart2) result(cart3)

    type(cart3d), intent(in) :: cart1, cart2
    type(cart3d) :: cart3
    
    ! compute the difference of two 3d cartesean points
    
    cart3%x = cart1%x - cart2%x
    cart3%y = cart1%y - cart2%y
    cart3%z = cart1%z - cart2%z
    
    return
end function cart3d_minus_cart3d
!}}}
end module view_modifier
