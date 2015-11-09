program graph_main
    

    use module_sphere

    integer    ::    results
    integer    ::    windW = 400
    integer    ::    windH = 400
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
    win = glutCreateWindow("Spheres")
  
!!graphical functions{{{
!! graphical functions
!    call myinit()
!    call glutReshapeFunc_gl(myreshape)
!    call glutDisplayFunc(mydisplay)
!    call glutMainLoop()
!    !}}}
!graphical functions{{{
! graphical functions
    call Init()
    call glutReshapeFunc_gl(Reshape)
    !call glutKeyboardFunc(Key)
    call glutDisplayFunc(Display)
    call glutMainLoop()
    !}}}

end program graph_main
