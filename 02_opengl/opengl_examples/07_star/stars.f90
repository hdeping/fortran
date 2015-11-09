program stars_prog
use opengl_gl
use opengl_glut
use stars_mod
  integer(GLenum) type
  integer(glcint) win

  call glutInitWindowSize(windW, windH)
  call glutInit()
  doubleBuffer = GL_TRUE !  Args(argc, argv);

  type = GLUT_RGB
  if (doubleBuffer == GL_TRUE) then
     type = ior(type,GLUT_DOUBLE)
  else
     type = ior(type,GLUT_SINGLE)
  end if
  call glutInitDisplayMode(type)
  win = glutCreateWindow("Stars")

  call Init()

  call glutReshapeFunc(Reshape)
  call glutKeyboardFunc(Key)
  call glutVisibilityFunc(Visible)
  call glutDisplayFunc(Display)
  call glutMainLoop()

stop
end program stars_prog
