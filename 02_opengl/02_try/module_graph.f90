module graph_module
    use opengl_gl
    use opengl_glu
    use opengl_glut
    !implicit none
    integer, parameter           :: graph_width  = int(500)
    integer, parameter           :: graph_height = int(500)
    integer, parameter           :: nlocal       = 100
    integer, parameter           :: numx         = 100
    integer                      :: i
    integer                      :: j
    integer                      :: k
    character(10)                :: filename 
    contains
!subroutine myinit{{{
subroutine myinit()
    call glclearcolor(1.0, 1.0, 1.0, 1.0)
    call glshademodel(gl_smooth)
    call glenable(gl_line_smooth)
    call glenable(gl_point_smooth)
    call glenable(gl_blend)
    call glblendfunc(gl_src_alpha, gl_one_minus_src_alpha)
    call glhint(gl_line_smooth_hint, gl_dont_care)
    call glhint(gl_point_smooth_hint, gl_dont_care)
    call gllinewidth(1.0)
    call glpointsize(2.0)
end subroutine myinit
!}}}
!subroutine draw_cell{{{
subroutine draw_cell( ip )
    integer               :: ip
    real                :: xx , xy
    real                  :: cc(3) , cr(3)
    integer               :: ix , iy
    call gllinestipple( 1 , #aaaaaa )            !...»­ÐéÏß
    call glenable(gl_line_stipple)
    do ix = 1 , numx-1
        call glbegin(gl_lines)
        call glcolor3f(0.0,0.0,0.0)
        call glvertex3f(1.*ix*discell,0.0,0.0)
        call glvertex3f(1.*ix*discell,lbox,0.0)
        call glend()
    enddo
    do ix = 1 , numx-1
        call glbegin(gl_lines)
        call glcolor3f(0.0,0.0,0.0)
        call glvertex3f(0.0,1.*ix*discell,0.0)
        call glvertex3f(lbox,1.*ix*discell,0.0)
        call glend()
    enddo
    call gldisable(gl_line_stipple)
    cc(1) = 1.0 ; cc(2)=0.0d0 ; cc(3)=0.0d0
    cr(1) = 0.0 ; cr(2)=0.0d0 ; cr(3)=1.0d0
    
    !call draw_little_box( ix-1 , iy-1 , cc )
    !do il = - 1 , 1 , 1
    !    do jl = - 1 , 1 , 1
    !        if(il==jl.and.il==0)cycle
    !        call draw_little_box(ix+il-1,iy+jl-1 , cr)
    !    enddo
    !enddo
end subroutine draw_cell
!}}}
!subroutine arrow_velocity{{{
subroutine arrow_velocity()
    real               :: length, angle
    real               :: c_max, c_min
    
    c_min=-log10(1.0e-6)
    c_max=-log10(0.1)
    do i=0, nlocal, 10
        do j=0, nlocal, 10
          ! if (la(i, j)%ls==1.and. dble(i)>sp(1)%po(1)-150.0.and.dble(i)<sp(1)%po(1)+30.0) then
            length=sqrt(vloc(i,j,1)**2.+vloc(i,j,2)**2.)
            if (length==0.0) then
                angle=0.0
            elseif (vloc(i,j,2)>=0.0.and.length/=0.0) then
                angle=acos(vloc(i,j,1)/length)
            elseif (vloc(i,j,2)<0.0.and.length/=0.0) then
                angle=2.0*pi-acos(vloc(i,j,1)/length)
            endif
            call glpushmatrix()
            call gltranslatef(dble(i*vstep), dble(j*vstep), 0.0)
            call glrotatef(360.0*angle/(2.0*pi), 0.0, 0.0, 1.0)
            if (length>1.0e-6) then
                call arrow(0.1*((c_min+log10(length))/(c_min-c_max)))
            endif
            call glpopmatrix()
        enddo
    enddo
end subroutine arrow_velocity
!}}}
!subroutine arrow{{{
subroutine arrow(length)
    real                :: length
    call glbegin(gl_lines)
    call glcolor3f(0.0,0.2,0.6)
    call glvertex3f(0.0, 0.0, 0.0)
    call glvertex3f(length, 0.0, 0.0)
    call glvertex3f(length, 0.0, 0.0)
    call glvertex3f(length-length/6.0, length/6.0, 0.0)
    call glvertex3f(length, 0.0, 0.0)
    call glvertex3f(length-length/6.0, -length/6.0, 0.0)
    call glend()
end subroutine arrow
!}}}
!subroutine draw_little_box{{{
subroutine draw_little_box( ipt , jpt ,c1 )
    integer           :: ipt ,jpt
    real              :: c1(3)
    call glbegin(gl_line_strip)
    call glcolor3fv(loc(c1))
    call glvertex3f(1.*ipt*discell,1.*jpt*discell,0.0)
    call glvertex3f(1.*(ipt+1)*discell,1.*jpt*discell,0.0)
    call glvertex3f(1.*(ipt+1)*discell,1.*(jpt+1)*discell,0.0)
    call glvertex3f(1.*ipt*discell,1.*(jpt+1)*discell,0.0)
    call glvertex3f(1.*ipt*discell,1.*jpt*discell,0.0)
    call glend()
end subroutine draw_little_box
!}}}
!subroutine draw_vkticle_simple{{{
subroutine draw_vkticle_simple( )
    real               :: lenv
    real               :: radius
    integer            :: nrd
    real               :: nsd
    real               :: thetav
    real               :: drawpos(nsum,2)
    real               :: cnormal(3) , cnormal1(3)
    real               :: ctest(3)  , ctest1(3)
    real               :: vacc(3)
    cnormal(1)  = 1.0 
    cnormal(2)  =  0. 
    cnormal(3)  = 0.0
    cnormal1(1) = 1.0 
    cnormal1(2) = 0.2 
    cnormal1(3) = 0.0
    ctest(1)    = 0.0 
    ctest(2)    = 0.0 
    ctest(3)    = 1.0
    ctest1(1)   = 0.0 
    ctest1(2)   = 0.2 
    ctest1(3)   = 1.0
    do ip =  1 , nsum
        drawpos(ip,:)=vk(ip)%p(:)
        call glpushmatrix()
        call glbegin( gl_points )
        call glvertex3f(vk(ip)%p(1),vk(ip)%p(2),0.0d0)
        call glend()
        call glpopmatrix()
    enddo
end subroutine draw_vkticle_active
!}}}
!subroutine draw_vkticle_active{{{
subroutine draw_vkticle_active( radius , nrd , lenv )
    real               :: lenv
    real               :: radius
    integer              :: nrd
    real               :: nsd
    real               :: thetav
    real               :: drawpos(nsum,2)
    real                 :: cnormal(3) , cnormal1(3)
    real                 :: ctest(3)  , ctest1(3)
    real               :: vacc(3)
    cnormal(1)  = 1.0d0 
    cnormal(2)  = 0.d0 
    cnormal(3)  = 0.0d0
    cnormal1(1) = 1.0d0 
    cnormal1(2) = 0.2d0 
    cnormal1(3) = 0.0d0
    ctest(1)    = 0.0 
    ctest(2)    = 0.0 
    ctest(3)    = 1.0
    ctest1(1)   = 0.0 
    ctest1(2)   = 0.2 
    ctest1(3)   = 1.0
    do ip = 1 , nsum
        drawpos(ip,:)=vk(ip)%p(:)
        call glbegin(gl_points)
        call glcolor3f(1.0,0.0,0.0)
        call glvertex3f(drawpos(ip,1),drawpos(ip,2),0.0)
        call glend()
        call glbegin( gl_lines )
        call glcolor3f(0.0,0.0,0.0)
        call glvertex3f(drawpos(ip,1),drawpos(ip,2),0.0)
        call glvertex3f(drawpos(ip,1)+radius*cos(vk(ip)%ag),&
                         drawpos(ip,2)+radius*sin(vk(ip)%ag),0.0)
        call glvertex3f(drawpos(ip,1)+radius*cos(vk(ip)%ag),&
                         drawpos(ip,2)+radius*sin(vk(ip)%ag),0.0)
        call glvertex3f(drawpos(ip,1)+0.8*radius*cos(pi/15.+vk(ip)%ag),&
                         drawpos(ip,2)+0.8*radius*sin(pi/15.+vk(ip)%ag),0.0)
        call glvertex3f(drawpos(ip,1)+radius*cos(vk(ip)%ag),&
                         drawpos(ip,2)+radius*sin(vk(ip)%ag),0.0)
        call glvertex3f(drawpos(ip,1)+0.8*radius*cos(-pi/15.+vk(ip)%ag),&
                         drawpos(ip,2)+0.8*radius*sin(-pi/15.+vk(ip)%ag),0.0)
        call glend()
        !call v_local( vk(ip)%p(1) , vk(ip)%p(2) , vacc(1:2) )
        call glbegin( gl_lines )
        call glcolor3f(0.0,0.0,0.0)
        call glend()
    enddo
end subroutine
!}}}
!subroutine draw_box{{{
subroutine draw_box()
    call glbegin( gl_line_strip )
        call glcolor3f(0.0 ,  0.0 , 0.0)
        call glvertex3f(0.0,0.0,0.0)
        call glvertex3f(lbox,0.0,0.0)
        call glvertex3f(lbox,lbox,0.0)
        call glvertex3f(0.0,lbox,0.0)
        call glvertex3f(0.0,0.0,0.0)
    call glend()
    end subroutine draw_box
!}}}
!subroutine myreshape{{{
subroutine myreshape(w, h)
    integer        ::        w, h
    
    call glviewport(0, 0, w, h)
    call glmatrixmode(gl_projection)
    call glloadidentity()
    call gluortho2d(-0.1, dble(lbox)+0.1, -0.1, dble(lbox)+0.1)
    
    call glmatrixmode(gl_modelview)
    call glloadidentity()
end subroutine myreshape
!}}}
!subroutine mydisplay{{{
subroutine mydisplay()
    integer step
    integer it
    integer ip
    integer jp
    integer ij
    integer ix
    
    do it = 1, int(2e2)
    !    print*,it
    enddo
    do it = 1 , 1
        step=step+1
        call viscek_move()
    enddo
    if(step>1000)then
    do ip =0 , nlocal
        do jp = 0 , nlocal
    !        vloc(ip,jp,1)=0.5*sin(vstep*ip)*cos(vstep*jp)
    !        vloc(ip,jp,2)=-0.5*cos(vstep*ip)*sin(vstep*jp)
        enddo
    enddo
    endif
    !    dt=0.0001
        
    if (mod(step, 1)==0) then
        do ij = 1 , 1
            call glclearcolor(1.0, 1.0, 1.0, 1.0)
            call glclear(gl_color_buffer_bit)        
            call glcolor3f(0.0, 0.0, 0.0)
        !    call cal_mcenter( )
            call draw_vkticle_active(0.2d0 , 32, 0.0d0 )
        !    call draw_vkticle_simple( )
        !    call draw_string( 0.5d0 , 32, 0.0d0)
            call draw_box()
            call draw_cell(1)
            call arrow_velocity()
        
        !    call draw_cell(2)
            call glutswapbuffers()
        enddo
    endif
    !    write(10,*)1,pcenter(:)
end subroutine mydisplay
!}}}
!subroutine keyregister{{{
subroutine keyregister()
   call glutkeyfunc( glut_0 , loc(triger_off) )
   call glutkeyfunc( glut_1 , loc(triger_on) )
   call glutkeyfunc( glut_2 , loc(triger_mod))
   call glutkeyfunc( glut_up, loc(move_up))
   call glutkeyfunc( glut_down, loc(move_down))
   call glutkeyfunc( glut_left, loc(move_left))
   call glutkeyfunc( glut_right, loc(move_right))
   call glutkeyfunc( glut_5    , loc(save_control))
end subroutine keyregister
!}}}
end module graph_module
