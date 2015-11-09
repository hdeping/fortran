

program blender_main
use opengl_glut
use blender
integer(glcint) i

  call glutInit()
  CALL glutInitWindowPosition( 200, 200 );
  call glutInitDisplayMode(ior(ior(GLUT_DOUBLE,GLUT_RGB),GLUT_DEPTH))
  i = glutCreateWindow("blender")
  call glutDisplayFunc(display)
  call glutVisibilityFunc(visible)

  call glNewList(1, GL_COMPILE)  ! create ico display list
  call glutSolidIcosahedron()
  call glEndList()

  call glEnable(GL_LIGHTING)
  call glEnable(GL_LIGHT0)
  call glLightfv(GL_LIGHT0, GL_AMBIENT, light0_ambient)
  call glLightfv(GL_LIGHT0, GL_DIFFUSE, light0_diffuse)
  call glLightfv(GL_LIGHT1, GL_DIFFUSE, light1_diffuse)
  call glLightfv(GL_LIGHT1, GL_POSITION, light1_position)
  call glLightfv(GL_LIGHT2, GL_DIFFUSE, light2_diffuse)
  call glLightfv(GL_LIGHT2, GL_POSITION, light2_position)
  call glEnable(GL_DEPTH_TEST)
  call glEnable(GL_CULL_FACE)
  call glEnable(GL_BLEND)
  call glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
  call glEnable(GL_LINE_SMOOTH)
  call glLineWidth(2.0)

  call glMatrixMode(GL_PROJECTION)
  call gluPerspective( 40.0_gldouble, & ! field of view in degree
                           1.0_gldouble, & ! aspect ratio
                           1.0_gldouble, & ! Z near
                          10.0_gldouble)   ! Z far
  call glMatrixMode(GL_MODELVIEW)
  call gluLookAt( &
     0.0_gldouble, 0.0_gldouble, 5.0_gldouble, & ! eye is at (0,0,5)
     0.0_gldouble, 0.0_gldouble, 0.0_gldouble, & ! center is at (0,0,0)
     0.0_gldouble, 1.0_gldouble, 0.0_gldouble)    ! up is in positive Y direction
  call glTranslatef(0.0, 0.6, -1.0)

  call glutMainLoop()
end program blender_main
