program main

use opengl_gl
use opengl_glut
use view_modifier
use view_demo_callbacks
implicit none

integer :: winid, menuid
real(kind=glfloat), dimension(3) :: & ! colors for bronze from Redbook teapots
   ambient = (/ 0.2125_glfloat, 0.1275_glfloat, 0.054_glfloat /), &
   diffuse = (/ 0.714_glfloat, 0.4284_glfloat, 0.18144_glfloat /), &
  specular = (/ 0.393548_glfloat, 0.271906_glfloat, 0.166721_glfloat /)
real(kind=glfloat), dimension(4) :: &
       pos = (/ 1.0_glfloat, 1.0_glfloat, 1.0_glfloat, 0.0_glfloat /), &
     white = (/ 1.0_glfloat, 1.0_glfloat, 1.0_glfloat, 1.0_glfloat /)

! Initializations


call glutInit

CALL glutInitWindowSize( 350, 350 );
CALL glutInitWindowPosition( 200, 200 );


call glutInitDisplayMode(ior(GLUT_DOUBLE,ior(GLUT_RGB,GLUT_DEPTH)))

! Create a window

winid = glutCreateWindow("View Modifier Demo")
menuid = view_modifier_init()
call glutAttachMenu(GLUT_RIGHT_BUTTON)

! Set the display callback

call glutDisplayFunc(display)

! Create the image

call glNewList(1,GL_COMPILE)

! Draw axes so we know the orientation

call glBegin(GL_LINES)
call glVertex3f(0.0_glfloat,0.0_glfloat,0.0_glfloat)
call glVertex3f(2.0_glfloat,0.0_glfloat,0.0_glfloat)
call glVertex3f(0.0_glfloat,0.0_glfloat,0.0_glfloat)
call glVertex3f(0.0_glfloat,2.0_glfloat,0.0_glfloat)
call glVertex3f(0.0_glfloat,0.0_glfloat,0.0_glfloat)
call glVertex3f(0.0_glfloat,0.0_glfloat,2.0_glfloat)

! Draw crude x, y and z to label the axes

call glVertex3f(2.1_glfloat,-0.1_glfloat,0.1_glfloat) ! X
call glVertex3f(2.1_glfloat,0.1_glfloat,-0.1_glfloat)
call glVertex3f(2.1_glfloat,-0.1_glfloat,-0.1_glfloat)
call glVertex3f(2.1_glfloat,0.1_glfloat,0.1_glfloat)
call glVertex3f(0.1_glfloat,2.1_glfloat,0.1_glfloat) ! Y
call glVertex3f(0.0_glfloat,2.1_glfloat,0.0_glfloat)
call glVertex3f(-0.1_glfloat,2.1_glfloat,0.1_glfloat)
call glVertex3f(0.1_glfloat,2.1_glfloat,-0.1_glfloat)
call glVertex3f(-0.1_glfloat,0.1_glfloat,2.1_glfloat) ! Z
call glVertex3f(0.1_glfloat,0.1_glfloat,2.1_glfloat)
call glVertex3f(0.1_glfloat,0.1_glfloat,2.1_glfloat)
call glVertex3f(-0.1_glfloat,-0.1_glfloat,2.1_glfloat)
call glVertex3f(-0.1_glfloat,-0.1_glfloat,2.1_glfloat)
call glVertex3f(0.1_glfloat,-0.1_glfloat,2.1_glfloat)
call glEnd

! Draw a teapot

! rotate so the z-axis comes out the top, x-axis out the spout
call glRotated(90.0_gldouble,1.0_gldouble,0.0_gldouble,0.0_gldouble)
call glMaterialfv(GL_FRONT, GL_AMBIENT, ambient)
call glMaterialfv(GL_FRONT, GL_DIFFUSE, diffuse)
call glMaterialfv(GL_FRONT, GL_SPECULAR, specular)
call glMaterialf(GL_FRONT, GL_SHININESS, 25.6_glfloat)
call glutSolidTeapot(1.0_gldouble)

call glEndList

! Set the lighting

call glClearColor(0.8_glclampf, 0.8_glclampf, 0.8_glclampf, 1.0_glclampf)
call glLightfv(GL_LIGHT0, GL_DIFFUSE, white)
call glLightfv(GL_LIGHT0, GL_POSITION, pos)
call glEnable(GL_LIGHTING)
call glEnable(GL_LIGHT0)
call glEnable(GL_DEPTH_TEST)

! Let glut take over

call glutMainLoop()

end program main
