module module_common
    use opengl_gl
    use opengl_glu
    use opengl_glut
    implicit none
    integer,parameter            :: n = 100
    integer                      :: i
    integer                      :: j
    integer                      :: k
    character(10)                :: filename 
    contains
!subroutine init{{{
subroutine init()
    call glclearcolor(0.0,0.0,0.0,0.0)
    call glshademodel(gl_flat)
end subroutine init
!}}}
!subroutine display{{{
subroutine display()
    call glclear()
    call glpushmatrix()
    call glrotatef(spin,0.0,0.0,1.0)
    call glcolor3f(1.0,1.0,1.0)
    call glrectf(- 25.0, - 25.0,25.0,25.0)
    call glpopmatrix()
    call glutswapbuffers()
end subroutine display
!}}}
!subroutine spindisplay{{{
subroutine spindisplay()
    spin = spin + 2.0
    if ( spin > 360.0 )then
        spin = spin - 360.0
        call glutpostwindowredisplay()
    endif ! if ends
    
end subroutine spindisplay
!}}}
!subroutine reshape{{{
subroutine reshape(w,h)
    integer,intent(in)         :: w
    integer,intent(in)         :: h
    call glviewport(0.0,w,h)
    call glmatrixmode(gl_projection)
    call glloadidentity()
    call glortho(- 50.0,50.0,- 50.0,50.0,- 1.0,1.0)
    call glmatrixmode(gl_modelview)
    call glloadidentity()
end subroutine reshape
!}}}
!subroutine mouse{{{
subroutine mouse(button,state,x,y)
    integer,intent(in)          ::  button
    integer,intent(in)          ::  state
    integer,intent(in)          ::  x
    integer,intent(in)          ::  y

    select case(button)
        case (glut_left_button)
            if ( state == glut_down )then
                call glutidlefunc_gl(spindisplay)
            endif ! if ends
        case (glut_middle_button)   
            if ( state == glut_down )then
                call glutidlefunc_gl(null)
        endif ! if ends

    end select
end subroutine mouse
!}}}
end module module_common
