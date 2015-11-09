module module_sphere
    use opengl_gl
    use opengl_glu
    use opengl_glut
    implicit none

    integer, parameter                :: pn             = 20
    integer, parameter                :: graph_width    = 800
    integer, parameter                :: graph_height   = 800
    real(4), parameter                :: pi             = 3.1415926
    real(4), parameter                :: l              = 10.0
    real(4), parameter                :: bigmax         = 1.0e8
    real(4)                           :: pp(pn, 2)   
    real(4)                           :: pv(pn, 2)   
    real(4)                           :: ra(pn)                  
    real(4)                           :: rin(pn)      
    real(4)                           :: cotime(pn)      
    integer                           :: partne(pn)      
    integer                           :: motion_step
    integer                           :: partner(pn)
    integer(GLenum)                   :: doubleBuffer
    
    contains
!subroutine initial_parameter{{{
! initial the parameters
subroutine initial_parameter()
    integer                ::                i
    integer                ::                j
    real(8)                ::                r
    
        ra       = 0.0
        rin      = 0.1
        do i=1, pn
            do j=1, 2
                call random_number(r)
                pp(i, j)=l*r
                call gaussian_noise(r)
                pv(i, j)=r
            enddo
        enddo
        
        call initial_time_list()
        motion_step=0
    
    end subroutine initial_parameter
!}}}
!subroutine evolution{{{
!  evolution of sphere collision
subroutine evolution(ii, jj)
    integer, dimension(1)      ::     it
    integer                    ::     i
    integer                    ::     j
    integer                    ::     ii
    integer                    ::     jj
    integer                    ::     k
    real(8)                    ::     rx(2)
    real(8)                    ::     irx(2)
    real(8)                    ::     vxpi(2)
    real(8)                    ::     vxpj(2)
    real(8)                    ::     vxvi(2)
    real(8)                    ::     vxvj(2)
    real(8)                    ::     vx(2)
    real(8)                    ::     aij
    real(8)                    ::     bij
    real(8)                    ::     cij
    real(8)                    ::     dij
    real(8)                    ::     tij
    real(8)                    ::     tvi
    real(8)                    ::     tvj

    !do 
        motion_step=motion_step+1
    
        it(:)=  minloc(cotime)
        tij  =  cotime(it(1))
        i    =  it(1)
        j    =  partner(i)
    
        do k = 1, pn
            cotime(k)=cotime(k)-tij
            pp(k, :)=pp(k, :)+pv(k, :)*tij
            where(pp(k, :)>l)
                pp(k, :)=pp(k, :)-l
            elsewhere(pp(k, :)<=0.0d0)
                pp(k, :)=pp(k, :)+l
            endwhere
            ra(k)=ra(k)+rin(k)*tij    
        enddo
        
        ii=i
        jj=j
        
        rx(:)=pp(ii, :)-pp(jj, :)
        where(abs(rx)>l/2.0d0)
            rx=pp(ii, :)-(pp(jj, :)-l*dsign(1.0d0, pp(jj, :)-l/2.0d0))
        endwhere
        irx(:)=rx(:)/dsqrt(rx(1)**2.+rx(2)**2.)    
        tvi=pv(ii, 1)*irx(1)+pv(ii, 2)*irx(2)
        tvj=pv(jj, 1)*irx(1)+pv(jj, 2)*irx(2)
        vxpi(:)=tvi*irx(:); vxpj(:)=tvj*irx(:)
        vxvi(:)=pv(ii, :)-vxpi(:); vxvj(:)=pv(jj, :)-vxpj(:)
        
        pv(ii, :)=(vxpj(:)+(rin(ii)+rin(jj))*irx(:))+vxvi(:)
        pv(jj, :)=(vxpi(:)-(rin(ii)+rin(jj))*irx(:))+vxvj(:)
        
        do i=1, pn
            if ((i==ii).or.(i==jj).or.(partner(i)==ii).or.(partner(i)==jj)) then
                cotime(i)=bigmax
                do j=1, pn
                    if (j/=i) then
                        rx(:)=pp(i, :)-pp(j, :)
                        where(abs(rx)>l/2.0d0)
                            rx=pp(i, :)-(pp(j, :)-l*dsign(1.0d0, pp(j, :)-l/2.0d0))
                        endwhere
                        vx(:)=pv(i, :)-pv(j, :)
                        
                        aij=(vx(1)**2.+vx(2)**2.)-(rin(i)+rin(j))**2.
                        bij=(rx(1)*vx(1)+rx(2)*vx(2))-(ra(i)+ra(j))*(rin(i)+rin(j))
                        cij=(rx(1)**2.+rx(2)**2.)-(ra(i)+ra(j))**2.
                        dij=bij**2.-aij*cij                    
                        
                        if ((bij<=0.0d0.or.aij<0.0d0).and.dij>=0.0d0) then
                            tij=(-bij-sqrt(dij))/aij
                            if (tij<cotime(i)) then
                                cotime(i)=tij
                                partner(i)=j
                            endif
                            if (tij<cotime(j)) then
                                cotime(j)=tij
                                partner(j)=i
                            endif
                        endif
    
                    endif
                enddo
            endif
        enddo
    
    !enddo
    
        !call dis_ij()   2015-09-23 10:56:21    
    
    end subroutine evolution
!}}}
!  hard sphere collision
!subroutine initial_time_list{{{
subroutine initial_time_list()
    integer                        ::                i, j
    real(8), dimension(2)          ::                rx, vx
    real(8)                        ::                aij
    real(8)                        ::                bij
    real(8)                        ::                cij
    real(8)                        ::                dij
    real(8)                        ::                tij
 
    
        cotime=bigmax
        
        do i=1, pn-1
            do j=i+1, pn
                
                rx(:)=pp(i, :)-pp(j, :)
                where(abs(rx)>l/2.0d0)
                    rx=pp(i, :)-(pp(j, :)-l*dsign(1.0d0, pp(j, :)-l/2.0d0))
                endwhere
                vx(:)=pv(i, :)-pv(j, :)
                
                aij=(vx(1)**2.+vx(2)**2.)-(rin(i)+rin(j))**2.
                bij=(rx(1)*vx(1)+rx(2)*vx(2))-(ra(i)+ra(j))*(rin(i)+rin(j))
                cij=(rx(1)**2.+rx(2)**2.)-(ra(i)+ra(j))**2.
                dij=bij**2.-aij*cij
    
                if ((bij<=0.0d0.or.aij<0.0d0).and.dij>=0.0d0) then
                    tij=(-bij-sqrt(dij))/aij
                    if (tij<cotime(i)) then
                        cotime(i)=tij
                        partner(i)=j
                    endif
                    if (tij<cotime(j)) then
                        cotime(j)=tij
                        partner(j)=i
                    endif
                endif
    
            enddo
        enddo
    
    
    end subroutine initial_time_list
    
!}}}
!subroutine gaussian_noise{{{
subroutine gaussian_noise(gauss_ran)
    real(8)                        ::                u, g
    integer                        ::                n
    integer                        ::                ij
    real(8)                        ::                bo, r
    real(8)                        ::                gauss_ran
    
    n  = 12
    bo = 0.0
    g  = 1.0
    u  = 0.0
    
        do ij=1, n
            call random_number(r)
            bo=bo+r
        enddo    
    
        gauss_ran=u+g*(bo-real(n)/2.0)                              
    
    end subroutine gaussian_noise
!}}}
!  graphical process
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
    end subroutine myinit
!}}}
!subroutine mydisplay{{{
subroutine mydisplay()
    integer                ::        i, is, js
    
            call glclear(gl_color_buffer_bit)            
            
            call evolution(is, js)
    
            !if (mod(motion_step, 400)==0) then
    
            do i=1, pn
                if (i==is.or.i==js) then
                    call glcolor3f(1.0, 1.0, 0.0)    
                else
                    call glcolor3f(1.0, 0.0, 0.0)
                endif
                call glpushmatrix()
                    call gltranslatef(pp(i, 1), pp(i, 2), 0.0)
                    call sphere_graph(ra(i), 24)
                call glpopmatrix()
                
            enddo
                print*, pp(is, :), pp(js, :), 'a'
                print*, sqrt((pp(is, 1)-pp(js, 1))**2.+(pp(is, 2)-pp(js, 2))**2.), ra(1)*2.0d0
                call glutswapbuffers()
            !    pause
        !    endif
    
    
    end subroutine mydisplay
!}}}
!subroutine Display{{{
subroutine Display() bind(c)
      !call ShowStars()
      call draw_particle_simple()
      if (doubleBuffer == GL_TRUE) then
        call glutSwapBuffers()
      else
        call glFlush()
      end if
    return
    end subroutine Display
    
!}}}
!subroutine myreshape{{{
subroutine myreshape(w, h)
    integer        ::        w, h
    
    call glviewport(0, 0, w, h)
    call glmatrixmode(gl_projection)
    call glloadidentity()
    call gluortho2d(0D0, 10D0, 0D0, 10D0)
    call glmatrixmode(gl_modelview)
    call glloadidentity()
    
    end subroutine myreshape
!}}}
!subroutine Reshape{{{
subroutine Reshape(width, height) bind(c)
    integer(glcint), intent(in)  :: width
    integer(glcint), intent(in)  :: height
    integer(glcint)              :: windW
    integer(glcint)              :: windH

      windW = width
      windH = height
    
      call glViewport(0_GLint, 0_GLint, windW, windH)
    
      call glMatrixMode(GL_PROJECTION)
      call glLoadIdentity()
      call gluOrtho2D(-0.5_gldouble, windW + 0.5_gldouble, -0.5_gldouble, windH + 0.5_gldouble)
      call glMatrixMode(GL_MODELVIEW)
    return
    end subroutine Reshape
    
!/* ARGSUSED1 */
!}}}
!subroutine sphere_graph{{{
subroutine sphere_graph(r, fr)
    real(4)        ::        r
    integer        ::        fr
    real(4)        ::        det_theta
    real(4)        ::        x
    real(4)        ::        y
    integer        ::        ss 
    integer        ::        i
    
    det_theta=2.0*pi/dble(fr)
    
    
        do i=1, fr
            call glbegin(gl_polygon)
                call glvertex3f(0.0, 0.0, 0.0)
                x  = r*cos(det_theta*dble(i-1))
                y  = r*sin(det_theta*dble(i-1))
                call glvertex3f(x,y,0.0)
                x  = r*cos(det_theta*dble(i))
                y  = r*sin(det_theta*dble(i))
                call glvertex3f(x,y,0.0)
            call glend()
        enddo
    
                call glcolor3f(0.0, 0.0, 0.0)    
    do i=1, fr
        call glbegin(gl_lines)
        x  = r*cos(det_theta*dble(i-1))
        y  = r*sin(det_theta*dble(i-1))
        call glvertex3f(x,y,0.0)
        x  = r*cos(det_theta*dble(i))
        y  = r*sin(det_theta*dble(i))
        call glvertex3f(x,y,0.0)
        call glend()
    enddo
    end subroutine sphere_graph
!}}}
!subroutine draw_particle_simple{{{
subroutine draw_particle_simple( )
real*4        ::        lenv
real*4        ::        radius
real*4        ::        nsd
real*4        ::        thetav
!real*4        ::        drawpos(nact,2)
real          ::        cnormal(3) , cnormal1(3)
real          ::        ctest(3)  , ctest1(3)
integer       ::        nrd
integer       ::        ip

cnormal(1) =   1.0d0 
cnormal(2) =   0.d0 
cnormal(3) =   0.0d0
cnormal1(1)=   1.0d0 
cnormal1(2)=   0.2d0 
cnormal1(3)=   0.0d0
ctest(1)   =   0.0 
ctest(2)   =   0.0 
ctest(3)   =   1.0
ctest1(1)  =   0.0 
ctest1(2)  =   0.2 
ctest1(3)  =   1.0

do ip =   1,pn !nstring+1 , nact
    !drawpos(ip,:)=par(ip)%p(:)
    call glpushmatrix()
        call glbegin( gl_points )
            !call glvertex3f(par(ip)%p(1),par(ip)%p(2),0.0d0)
            call glvertex3f(pp(ip,1),pp(ip,2),0.0)
        call glend()
    call glpopmatrix()
enddo

end subroutine draw_particle_simple

!}}}
end module module_sphere
