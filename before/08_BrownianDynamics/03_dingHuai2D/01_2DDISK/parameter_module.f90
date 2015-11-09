module parameter_module
    integer, parameter                        ::  pn=1000, pna=500
    real*8, parameter                         ::  pi=3.1415926
    real*8, parameter                         ::  l=11.1d0
    real*8, dimension(1:pn, 2)                ::  pp, pv
    real*8, dimension(1:pn)                   ::  ra                
    real*8, dimension(1:pn)                   ::  rin
    real*8, dimension(1:pn)                   ::  cotime
    integer, dimension(1:pn)                  ::  partner
    real*8, parameter                         ::  bigmax=1.0e8
    integer                                   ::  motion_step
    real*8                                    ::  time

    contains
!subroutine evolution{{{
subroutine evolution()
    integer, dimension(1)                ::                it
    integer                                ::                i, j, ii, jj
    real*8, dimension(2)                ::                rx, irx, vxpi, vxpj, vxvi, vxvj, vx
    real*8                                ::                aij, bij, cij, dij, tij
    real*8                                ::                as, packphi                        
    
    do 
        
        motion_step=motion_step+1
    
        if (mod(motion_step, 1000)==0) then
            as=0.0d0
            do i=1, pn
                as=as+pi*ra(i)**2.
            enddo
            packphi=as/dble(l*l)
            print*, packphi, time
        !    if (packphi>=0.80) then
        !        open(1, file='ini.txt')
        !        write(1, *) packphi
        !        do i=1, pn
        !            write(1, '(i5, 3f25.16)') i, ra(i), pp(i, :)
        !        enddo
        !        close(1)
        !        stop
        !    endif
        endif
    
    
    
        it(:)=minloc(cotime)
        tij=cotime(it(1))
        i=it(1)
        j=partner(i)
    
    
    
        do k=1, pn
            cotime(k)=cotime(k)-tij
            pp(k, :)=pp(k, :)+pv(k, :)*tij
            where(pp(k, :)>l)
                pp(k, :)=pp(k, :)-l
            elsewhere(pp(k, :)<=0.0d0)
                pp(k, :)=pp(k, :)+l
            endwhere
            ra(k)=ra(k)+rin(k)*tij    
        enddo
        
        time=time+tij
        if (time>=2.0d0) then
            as=0.0d0
            do i=1, pn
                as=as+pi*ra(i)**2.
            enddo
            packphi=as/dble(l*l)
            open(1, file='ini.txt')
            write(1, *) packphi
            do i=1, pn
                write(1, '(i5, 3f25.16)') i, ra(i), pp(i, :)
            enddo
            close(1)
            stop
        endif
    
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
    
    enddo
    
    !    call dis_ij()
    
    end subroutine evolution
!}}}
!subroutine dis_ij{{{
subroutine dis_ij()
    real*8                    ::                disf, dis
    real*8, dimension(2)    ::                lx
    integer                    ::                i, j
    
        disf=100000.0d0
        do i=1, pn-1
            do j=i+1, pn
                lx(:)=abs(pp(i, :)-pp(j, :))
                where(lx>l/2.0d0)
                    lx=l-lx
                endwhere
                dis=dsqrt(lx(1)**2.+lx(2)**2.)
                if (dis<disf) then
                    disf=dis
                endif
            enddo
        enddo
    
        print*, disf, 2*ra(1)
    
    end subroutine dis_ij
!}}}
!subroutine initial_parameter{{{
subroutine initial_parameter()
    integer                ::                i
    real*8                ::                r
    
        ra(:)=0.0d0; rin(1:pna)=0.1d0; rin(pna+1:1000)=0.070d0
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
        time=0.0d0
    
    end subroutine initial_parameter
!}}}
!subroutine initial_time_list{{{
subroutine initial_time_list()
    integer                        ::                i, j
    real*8, dimension(2)        ::                rx, vx
    real*8                        ::                aij, bij, cij, dij, tij
    
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
    real*8                        ::                u, g
    integer                        ::                n
    real*8                        ::                bo, r
    real*8                        ::                gauss_ran
    
    n=12; bo=0.0
    g=1.0; u=0.0
    
        do ij=1, n
            call random_number(r)
            bo=bo+r
        enddo    
    
        gauss_ran=u+g*(bo-real(n)/2.0)                              
    
    end subroutine gaussian_noise
!}}}
end module parameter_module
