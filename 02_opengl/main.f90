program main
    use module_common
    implicit none

    call glutinit_gl()
    call glutinitdisplaymode(glut_double,glut_rgb)
    call glutinitwindowsize(250,250)
    call glutinitwindowposition(100,100)
    call glutcreatewindow()
    call init()
    call glutdisplayfunc_gl(display)
    call glutreshapefunc_gl(reshape)
    call glutmousefunc(mouse)
    call glutmainloop()
end program main
