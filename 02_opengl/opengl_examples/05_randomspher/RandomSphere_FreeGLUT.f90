module opengl_example
   use iso_c_binding
   use opengl_gl
   use opengl_glu
   use opengl_glut
   implicit none
   private

   type, public :: spinningsphere
      type(c_ptr) :: quadric=c_null_ptr
      integer(kind=glint) :: gl_list=-1
      real(glfloat) :: radius=1.0_glfloat ! it will change randomly
   end type
   
   public :: testgl

contains

   subroutine display() bind(c) ! private so no binding label
      ! display glut callback
      
      type(c_ptr) :: handle
      type(spinningsphere), pointer :: sphere
      
      handle=glutgetwindowdata() ! a glut extension
      call c_f_pointer(cptr=handle, fptr=sphere)
   
      call glclear(ior(gl_color_buffer_bit,gl_depth_buffer_bit))
      call glpushmatrix()
      call glscalef(sphere%radius, sphere%radius, sphere%radius)                  
      call glcalllist(sphere%gl_list)      
      call glpopmatrix()
      call glutswapbuffers()
      
   end subroutine      
   
   subroutine idle() bind(c) ! private so no binding label
      ! idle glut callback
      
      type(c_ptr) :: handle
      type(spinningsphere), pointer :: sphere
      real(glfloat) :: dice
      
      handle=glutgetwindowdata() ! a glut extension
      call c_f_pointer(cptr=handle, fptr=sphere)
      
      call random_number(dice)
      sphere%radius=abs(1.0_glfloat+0.01_glfloat*(dice-0.5_glfloat))*sphere%radius
      call glutpostredisplay()   
      
   end subroutine      

   subroutine reshape(width, height) bind(c)
      ! reshape glut callback
      integer(glsizei), value :: width, height
      
      type(c_ptr) :: handle
      type(spinningsphere), pointer :: sphere
      
      handle=glutgetwindowdata() ! a glut extension
      call c_f_pointer(cptr=handle, fptr=sphere)

      call glviewport (0_glint, 0_glint, width, height)

   end subroutine    
      
   subroutine testgl(sphere)
      type(spinningsphere), intent(inout), target :: sphere
      character(kind=c_char, len=10) :: window_name="sphere"//c_null_char
      integer(glint) :: gl_window
      
      ! we do not pass command arguments for simplicity
      call glutinit()
      call glutinitdisplaymode(ior(glut_double,glut_rgb))
      gl_window=glutcreatewindow(window_name)
      call glutsetwindowdata(c_loc(sphere))

      sphere%gl_list=glgenlists(1)
      call glnewlist(sphere%gl_list, gl_compile)      
      sphere%quadric=glunewquadric()
      call gluquadricdrawstyle(sphere%quadric, glu_fill)
      call glusphere(sphere%quadric, 1.0_gldouble, 25_glint, 25_glint)
      call glendlist()

      call gllightfv(gl_light0, gl_diffuse, real((/1.0, 0.0, 0.0, 1.0/), glfloat))

      !call glenable(gl_lighting)
      !call glenable(gl_light0)
      !call glenable(gl_depth_test)

      ! set the viewing parameters (is this really needed?)
      call glmatrixmode(gl_projection)
      call gluperspective(40.0_gldouble, 1.0_gldouble, 1.0_gldouble, 10.0_gldouble)
      call glmatrixmode(gl_modelview)
      call glulookat(0.0_gldouble, 0.0_gldouble, 5.0_gldouble, &
                     0.0_gldouble, 0.0_gldouble, 0.0_gldouble, &
                     0.0_gldouble, 1.0_gldouble, 1.0_gldouble)
      call gltranslatef(0.0, 0.0, -1.0)               
      
      ! set callbacks
      call glutdisplayfunc(display) 
      call glutreshapefunc(reshape) 
      call glutidlefunc(idle)

      call glutmainloop() ! classical glut won't return
      
      call gludeletequadric(sphere%quadric) ! avoid memory leaks
      
      write(*,*) "glutmainloop returned!"

   end subroutine

end module


