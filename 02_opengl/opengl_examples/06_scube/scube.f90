program main
  use scube_mod
  implicit none

  integer :: width = 350, height = 350
  integer i, win
  character(len=30) name
  integer fog_menu
  real(glfloat) rGL_EXP

  call glutInitWindowSize(width, height)
  call glutInit

  ! choose visual
  if (useRGB) then
    if (useDB) then
      call glutInitDisplayMode(ior(ior(GLUT_DOUBLE,GLUT_RGB),GLUT_DEPTH))
      name = windowNameRGBDB
    else
      call glutInitDisplayMode(ior(ior(GLUT_SINGLE,GLUT_RGB),GLUT_DEPTH))
      name = windowNameRGB
    endif
  else
    if (useDB) then
      call glutInitDisplayMode(ior(ior(GLUT_DOUBLE,GLUT_INDEX),GLUT_DEPTH))
      name = windowNameIndexDB
    else
      call glutInitDisplayMode(ior(ior(GLUT_SINGLE,GLUT_INDEX),GLUT_DEPTH))
      name = windowNameIndex
    endif
  endif

  win = glutCreateWindow(name)

  call buildColormap()

  call glutKeyboardFunc(keyboard)
  call glutDisplayFunc(display)
  call glutVisibilityFunc(visible)

  fog_menu = glutCreateMenu(fog_select)
  call glutAddMenuEntry(CString("Linear fog"), GL_LINEAR)
  call glutAddMenuEntry(CString("Exp fog"), GL_EXP)
  call glutAddMenuEntry(CString("Exp^2 fog"), GL_EXP2)

  i = glutCreateMenu(menu_select)
  call glutAddMenuEntry(CString("Start motion"), 1)
  call glutAddMenuEntry(CString("Stop motion"), 2)
  call glutAddMenuEntry(CString("Toggle fog"), 3)
  call glutAddMenuEntry(CString("Toggle lighting"), 4)
  call glutAddSubMenu(CString("Fog type"), fog_menu)
  call glutAddMenuEntry(CString("Quit"), 5)
  call glutAttachMenu(GLUT_RIGHT_BUTTON)

  ! setup context
  call glMatrixMode(GL_PROJECTION)
  call glLoadIdentity()
  call glFrustum(-1.0_gldouble, 1.0_gldouble, -1.0_gldouble, &
                     1.0_gldouble, 1.0_gldouble, 3.0_gldouble)

  call glMatrixMode(GL_MODELVIEW)
  call glLoadIdentity()
  call glTranslatef(0.0, 0.0, -2.0)

  call glEnable(GL_DEPTH_TEST)

  if (useLighting) then
    call glEnable(GL_LIGHTING)
  endif
  call glEnable(GL_LIGHT0)
  call glLightfv(GL_LIGHT0, GL_POSITION, lightPos)
  call glLightfv(GL_LIGHT0, GL_AMBIENT, lightAmb)
  call glLightfv(GL_LIGHT0, GL_DIFFUSE, lightDiff)
  call glLightfv(GL_LIGHT0, GL_SPECULAR, lightSpec)
  
   ! call glLightfv(GL_LIGHT0, GL_SPOT_DIRECTION, lightDir);
   ! call glLightf(GL_LIGHT0, GL_SPOT_EXPONENT, 80);
   ! call glLightf(GL_LIGHT0, GL_SPOT_CUTOFF, 25);

  call glEnable(GL_NORMALIZE)

  if (useFog) then
    call glEnable(GL_FOG)
  endif
  call glFogfv(GL_FOG_COLOR, fogColor)
  call glFogfv(GL_FOG_INDEX, fogIndex)
  rGL_EXP = GL_EXP
  call glFogf(GL_FOG_MODE, rGL_EXP)
  call glFogf(GL_FOG_DENSITY, 0.5)
  call glFogf(GL_FOG_START, 1.0)
  call glFogf(GL_FOG_END, 3.0)

  call glEnable(GL_CULL_FACE)
  call glCullFace(GL_BACK)

  call glShadeModel(GL_SMOOTH)

  call glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
  if (useLogo) then
    call glPolygonStipple(sgiPattern)
  else
    call glPolygonStipple(shadowPattern)
  endif

  call glClearColor(0.0, 0.0, 0.0, 1.0)
  call glClearIndex(0.)
  call glClearDepth(1._gldouble)

  call glutMainLoop()
end program main
