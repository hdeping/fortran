module module_common
    integer, parameter           ::  pn=949
    !integer, parameter            ::  pna=300, pnb=pn-pna
    real(8), parameter           ::  pi=3.1415926
    real(8), parameter           ::  ll=18.0d0, vl=3.0d0
    real(8), dimension(3)        ::  le
    real(8), dimension(1:pn, 3)  ::  pp, pv
    real(8), dimension(1:pn)     ::  ra                
    real(8), dimension(1:pn)     ::  rin
    real(8), dimension(1:pn)     ::  cotime
    real(8), dimension(1:pn)     ::  botime
    integer, dimension(1:pn)     ::  partner
    integer, dimension(1:pn)     ::  bounder      !l=-1;r=1;f=2;b=-2;t=3;u=-3
    real(8), parameter            ::  bigmax=1.0e8
    integer                      ::  motion_step
    integer                      ::  iit, jjt
    real(8)                       ::  time


    contains
!  evolution
!subroutine evolution{{{
subroutine evolution()
integer, dimension(1)                ::                itc, itb
real(8)                               ::                tijc, tijb, r

    motion_step=motion_step+1

if (motion_step<10000) then
    
    itc(:)=minloc(cotime)
    tijc=cotime(itc(1))
    itb(:)=minloc(botime)
    tijb=botime(itb(1))
    
    if (tijc>=tijb) then
        time=time+tijb
        call boundary_event(itb(1), tijb)
    else
        time=time+tijc
        call hardsphere_event(itc(1), tijc)
    endif

    if (mod(motion_step, 100)==0) then
        !print*, time
        call dis_ij()
    endif
    !if (time>=5.0d0) then
    !    pause
    !endif

else

    if (motion_step==10000) then
        do i=1, pn
            do j=1, 3
                call gaussian_noise(r)
                pv(i, j)=r
            enddo    
        enddo
        ra(:)=ra(:)-0.01d0*(rin(:))
        time=time-0.01d0
        call periodic_initial_time_list()
    endif

    itc(:)=minloc(cotime)
    tijc=cotime(itc(1))
    itb(:)=minloc(botime)
    tijb=botime(itb(1))

    if (tijc>=tijb) then
        time=time+tijb
        call periodic_boundary_event(itb(1), tijb)
    else
        time=time+tijc
        call periodic_hardsphere_event(itc(1), tijc)
    endif
    if (mod(motion_step, 100)==0) then
        !print*, time
        call dis_ij()
    endif

    !if (time>=5.0d0) then
    !    pause
    !endif

endif

end subroutine evolution
!}}}
!subroutine dis_ij{{{
subroutine dis_ij()
real(8)                    ::                disf, dis
real(8), dimension(3)    ::                lx
integer                    ::                i, j, ii, jj


    disf=100000.0d0
    do i=1, pn-1
        do j=i+1, pn
            lx(:)=abs(pp(i, :)-pp(j, :))
            where(lx(1:3)>ll/2.0d0)
                lx(1:3)=ll-lx(1:3)
            endwhere
            dis=dsqrt(lx(1)**2.+lx(2)**2.+lx(3)**2.)
            if (dis<disf) then
                disf=dis
                iit=i
                jjt=j
            endif
        enddo
    enddo

    !print*, motion_step, disf, 2.0d0*ra(1), 2.0d0*ra(pn)

end subroutine dis_ij
!}}}
! initial some process
!subroutine initial_parameter{{{
subroutine initial_parameter()
integer                ::                i
real(8)                ::                r
    
    ra(:)=0.0d0; rin(1:pn)=0.1d0 ! rin(301:pn)=0.05d0

!    do i=1, pn
!       call gaussian_noise(r)
!       rin(i)=0.1d0*(1.0d0-0.1d0*r)
!   enddo

    le=(/ll, ll, vl/)
    do i=1, pn
        do j=1, 3
            call random_number(r)
            pp(i, j)=0.1d0+(le(j)-0.2d0)*r
            call gaussian_noise(r)
            pv(i, j)=r
        enddo    
    enddo
    
    call initial_time_list()
    motion_step=0

end subroutine initial_parameter
!}}}
!subroutine initial_time_list{{{
subroutine initial_time_list()
integer                        ::                i, j
real(8), dimension(3)        ::                rx, vx
real(8)                        ::                aij, bij, cij, dij, tij

    cotime=bigmax
    botime=bigmax

    do i=1, pn-1
        do j=i+1, pn

            rx(:)=pp(i, :)-pp(j, :)
            vx(:)=pv(i, :)-pv(j, :)
            
            aij=(vx(1)**2.+vx(2)**2.+vx(3)**2.)-(rin(i)+rin(j))**2.
            bij=(rx(1)*vx(1)+rx(2)*vx(2)+rx(3)*vx(3))-(ra(i)+ra(j))*(rin(i)+rin(j))
            cij=(rx(1)**2.+rx(2)**2.+rx(3)**2.)-(ra(i)+ra(j))**2.
            dij=bij**2.-aij*cij
            if (dij>=0.0d0) then
                tij=(-bij-dsqrt(dij))/aij
                if (tij>=0.0d0) then
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
    enddo

    do i=1, pn
        do j=1, 3
            if (-pv(i, j)+rin(i)/=0.0d0) then
                tij=(pp(i, j)-ra(i))/(-pv(i, j)+rin(i))
                if (tij>=0.0d0) then
                    if (tij<botime(i)) then
                        botime(i)=tij
                        bounder(i)=-j
                    endif
                endif
            endif
        enddo
        do j=1, 3
            if (pv(i, j)+rin(i)/=0.0d0) then
                tij=(le(j)-(pp(i, j)+ra(i)))/(pv(i, j)+rin(i))
                if (tij>=0.0d0) then
                    if (tij<botime(i)) then
                        botime(i)=tij
                        bounder(i)=j
                    endif
                endif
            endif
        enddo
    enddo

end subroutine initial_time_list
!}}}
!subroutine periodic_initial_time_list{{{
subroutine periodic_initial_time_list()
integer                        ::                i, j
real(8), dimension(3)        ::                rx, vx
real(8)                        ::                aij, bij, cij, dij, tij

    cotime=bigmax
    botime=bigmax

    do i=1, pn-1
        do j=i+1, pn

            rx(:)=pp(i, :)-pp(j, :)
            where(abs(rx(1:3))>ll/2.0d0)
                rx(1:3)=pp(i, 1:3)-(pp(j, 1:3)-ll*dsign(1.0d0, pp(j, 1:3)-ll/2.0d0))
            endwhere

            vx(:)=pv(i, :)-pv(j, :)
            
            aij=(vx(1)**2.+vx(2)**2.+vx(3)**2.)-(rin(i)+rin(j))**2.
            bij=(rx(1)*vx(1)+rx(2)*vx(2)+rx(3)*vx(3))-(ra(i)+ra(j))*(rin(i)+rin(j))
            cij=(rx(1)**2.+rx(2)**2.+rx(3)**2.)-(ra(i)+ra(j))**2.
            dij=bij**2.-aij*cij
            if (dij>=0.0d0) then
                tij=(-bij-dsqrt(dij))/aij
                if (tij>=0.0d0) then
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
    enddo

    do i=1, pn
        do j=3, 3
            if (-pv(i, j)+rin(i)/=0.0d0) then
                tij=(pp(i, j)-ra(i))/(-pv(i, j)+rin(i))
                if (tij>=0.0d0) then
                    if (tij<botime(i)) then
                        botime(i)=tij
                        bounder(i)=-j
                    endif
                endif
            endif
        enddo
        do j=3, 3
            if (pv(i, j)+rin(i)/=0.0d0) then
                tij=(le(j)-(pp(i, j)+ra(i)))/(pv(i, j)+rin(i))
                if (tij>=0.0d0) then
                    if (tij<botime(i)) then
                        botime(i)=tij
                        bounder(i)=j
                    endif
                endif
            endif
        enddo
    enddo

end subroutine periodic_initial_time_list
!}}}
!subroutine gaussian_noise{{{
subroutine gaussian_noise(gauss_ran)
real(8)                        ::                u, g
integer                        ::                n
real(8)                        ::                bo, r
real(8)                        ::                gauss_ran

n=12; bo=0.0
g=1.0; u=0.0

    do ij=1, n
        call random_number(r)
        bo=bo+r
    enddo    

    gauss_ran=u+g*(bo-real(n)/2.0)                              

end subroutine gaussian_noise
!}}}
! pattern module
!subroutine boundary_event{{{
subroutine boundary_event(ii, tij)
integer                        ::                    ii, jj
real(8)                        ::                    tij, aij, bij, cij, dij
real(8), dimension(3)        ::                    rx, vx



    do k=1, pn
        botime(k)=botime(k)-tij
        cotime(k)=cotime(k)-tij
        pp(k, :)=pp(k, :)+pv(k, :)*tij
        ra(k)=ra(k)+rin(k)*tij    
    enddo    

    jj=bounder(ii)

    if (jj<0) then
        if (pv(ii, abs(jj))<=0.0d0) then
            pv(ii, abs(jj))=-pv(ii, abs(jj))+rin(ii)
        elseif (pv(ii, abs(jj))>0.0d0) then
            pv(ii, abs(jj))=pv(ii, abs(jj))+rin(ii)
        endif
    elseif (jj>0) then
        if (pv(ii, jj)>=0.0d0) then
            pv(ii, jj)=-pv(ii, jj)-rin(ii)
        elseif (pv(ii, jj)<0.0d0) then
            pv(ii, jj)=pv(ii, jj)-rin(ii)
        endif
    endif
    
    do i=1, pn
        if ((i==ii).or.(partner(i)==ii)) then
            cotime(i)=bigmax
            do j=1, pn
                if (j/=i) then

                    rx(:)=pp(i, :)-pp(j, :)
                    vx(:)=pv(i, :)-pv(j, :)
                    
                    aij=(vx(1)**2.+vx(2)**2.+vx(3)**2.)-(rin(i)+rin(j))**2.
                    bij=(rx(1)*vx(1)+rx(2)*vx(2)+rx(3)*vx(3))-(ra(i)+ra(j))*(rin(i)+rin(j))
                    cij=(rx(1)**2.+rx(2)**2.+rx(3)**2.)-(ra(i)+ra(j))**2.
                    dij=bij**2.-aij*cij
                    if (dij>=0.0d0) then                    
                        tij=(-bij-dsqrt(dij))/aij
                        if (tij>=0.0d0) then
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
                endif
            enddo
        endif
    enddo
    
    do i=1, pn
        if (i==ii) then
            botime(i)=bigmax
            do j=-3, -1
                if (j/=jj) then
                    if (-pv(i, -j)+rin(i)/=0.0d0) then
                        tij=(pp(i, -j)-ra(i))/(-pv(i, -j)+rin(i))
                        if (tij>=0.0d0) then
                            if (tij<botime(i)) then
                                botime(i)=tij
                                bounder(i)=j
                            endif
                        endif
                    endif
                endif
            enddo
            do j=1, 3
                if (j/=jj) then
                    if (pv(i, j)+rin(i)/=0.0d0) then
                        tij=(le(j)-(pp(i, j)+ra(i)))/(pv(i, j)+rin(i))
                        if (tij>=0.0d0) then
                            if (tij<botime(i)) then
                                botime(i)=tij
                                bounder(i)=j
                            endif
                        endif
                    endif
                endif
            enddo
        endif
    enddo

end subroutine boundary_event
!}}}
!subroutine hardsphere_event{{{
subroutine hardsphere_event(ii, tij)
integer                        ::                ii, jj, i, j
real(8)                        ::                tij, aij, bij, cij, dij, tvi, tvj
real(8), dimension(3)        ::                rx, vx, irx, vxpi, vxpj, vxvi, vxvj

    do k=1, pn
        cotime(k)=cotime(k)-tij
        botime(k)=botime(k)-tij
        pp(k, :)=pp(k, :)+pv(k, :)*tij
        ra(k)=ra(k)+rin(k)*tij    
    enddo
    
    jj=partner(ii)
    
    rx(:)=pp(ii, :)-pp(jj, :)

    irx(:)=rx(:)/dsqrt(rx(1)**2.+rx(2)**2.+rx(3)**2.)    
    tvi=pv(ii, 1)*irx(1)+pv(ii, 2)*irx(2)+pv(ii, 3)*irx(3)
    tvj=pv(jj, 1)*irx(1)+pv(jj, 2)*irx(2)+pv(jj, 3)*irx(3)
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
                    vx(:)=pv(i, :)-pv(j, :)
                    
                    aij=(vx(1)**2.+vx(2)**2.+vx(3)**2.)-(rin(i)+rin(j))**2.
                    bij=(rx(1)*vx(1)+rx(2)*vx(2)+rx(3)*vx(3))-(ra(i)+ra(j))*(rin(i)+rin(j))
                    cij=(rx(1)**2.+rx(2)**2.+rx(3)**2.)-(ra(i)+ra(j))**2.
                    dij=bij**2.-aij*cij
                    if (dij>=0.0d0) then                    
                        tij=(-bij-dsqrt(dij))/aij
                        if (tij>=0.0d0) then
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
                endif
            enddo
        endif
    enddo

    do i=1, pn
        if (i==ii.or.i==jj) then
            botime(i)=bigmax
            do j=-3, -1
                if (-pv(i, -j)+rin(i)/=0.0d0) then
                    tij=(pp(i, -j)-ra(i))/(-pv(i, -j)+rin(i))
                    if (tij>=0.0d0) then
                        if (tij<botime(i)) then
                            botime(i)=tij
                            bounder(i)=j
                        endif
                    endif
                endif
            enddo
            do j=1, 3
                if (pv(i, j)+rin(i)/=0.0d0) then
                    tij=(le(j)-(pp(i, j)+ra(i)))/(pv(i, j)+rin(i))
                    if (tij>=0.0d0) then
                        if (tij<botime(i)) then
                            botime(i)=tij
                            bounder(i)=j
                        endif
                    endif
                endif
            enddo            
        endif
    enddo


end subroutine hardsphere_event
!}}}
!subroutine periodic_boundary_event{{{
subroutine periodic_boundary_event(ii, tij)
integer                        ::                    ii, jj
real(8)                        ::                    tij, aij, bij, cij, dij
real(8), dimension(3)        ::                    rx, vx



    do k=1, pn
        botime(k)=botime(k)-tij
        cotime(k)=cotime(k)-tij
        pp(k, :)=pp(k, :)+pv(k, :)*tij
        where(pp(k, :)>ll)
            pp(k, :)=pp(k, :)-ll
        elsewhere(pp(k, :)<=0.0d0)
            pp(k, :)=pp(k, :)+ll
        endwhere
        ra(k)=ra(k)+rin(k)*tij    
    enddo    

    jj=bounder(ii)

    if (jj<0) then
        if (pv(ii, abs(jj))<=0.0d0) then
            pv(ii, abs(jj))=-pv(ii, abs(jj))+rin(ii)
        elseif (pv(ii, abs(jj))>0.0d0) then
            pv(ii, abs(jj))=pv(ii, abs(jj))+rin(ii)
        endif
    elseif (jj>0) then
        if (pv(ii, jj)>=0.0d0) then
            pv(ii, jj)=-pv(ii, jj)-rin(ii)
        elseif (pv(ii, jj)<0.0d0) then
            pv(ii, jj)=pv(ii, jj)-rin(ii)
        endif
    endif
    
    do i=1, pn
        if ((i==ii).or.(partner(i)==ii)) then
            cotime(i)=bigmax
            do j=1, pn
                if (j/=i) then

                    rx(:)=pp(i, :)-pp(j, :)
                    where(abs(rx(1:2))>ll/2.0d0)
                        rx(1:2)=pp(i, 1:2)-(pp(j, 1:2)-ll*dsign(1.0d0, pp(j, 1:2)-ll/2.0d0))
                    endwhere

                    vx(:)=pv(i, :)-pv(j, :)
                    
                    aij=(vx(1)**2.+vx(2)**2.+vx(3)**2.)-(rin(i)+rin(j))**2.
                    bij=(rx(1)*vx(1)+rx(2)*vx(2)+rx(3)*vx(3))-(ra(i)+ra(j))*(rin(i)+rin(j))
                    cij=(rx(1)**2.+rx(2)**2.+rx(3)**2.)-(ra(i)+ra(j))**2.
                    dij=bij**2.-aij*cij
                    if (dij>=0.0d0) then                    
                        tij=(-bij-dsqrt(dij))/aij
                        if (tij>=0.0d0) then
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
                endif
            enddo
        endif
    enddo
    
    do i=1, pn
        if (i==ii) then
            botime(i)=bigmax
            do j=-3, -3
                if (j/=jj) then
                    if (-pv(i, -j)+rin(i)/=0.0d0) then
                        tij=(pp(i, -j)-ra(i))/(-pv(i, -j)+rin(i))
                        if (tij>=0.0d0) then
                            if (tij<botime(i)) then
                                botime(i)=tij
                                bounder(i)=j
                            endif
                        endif
                    endif
                endif
            enddo
            do j=3, 3
                if (j/=jj) then
                    if (pv(i, j)+rin(i)/=0.0d0) then
                        tij=(le(j)-(pp(i, j)+ra(i)))/(pv(i, j)+rin(i))
                        if (tij>=0.0d0) then
                            if (tij<botime(i)) then
                                botime(i)=tij
                                bounder(i)=j
                            endif
                        endif
                    endif
                endif
            enddo
        endif
    enddo

end subroutine periodic_boundary_event
!}}}
!subroutine periodic_hardsphere_event{{{
subroutine periodic_hardsphere_event(ii, tij)
integer                        ::                ii, jj, i, j
real(8)                        ::                tij, aij, bij, cij, dij, tvi, tvj
real(8), dimension(3)        ::                rx, vx, irx, vxpi, vxpj, vxvi, vxvj

    do k=1, pn
        cotime(k)=cotime(k)-tij
        botime(k)=botime(k)-tij
        pp(k, :)=pp(k, :)+pv(k, :)*tij
        where(pp(k, :)>ll)
            pp(k, :)=pp(k, :)-ll
        elsewhere(pp(k, :)<=0.0d0)
            pp(k, :)=pp(k, :)+ll
        endwhere
        ra(k)=ra(k)+rin(k)*tij    
    enddo
    
    jj=partner(ii)
    
    rx(:)=pp(ii, :)-pp(jj, :)

    where(abs(rx(1:3))>ll/2.0d0)
        rx(1:3)=pp(ii, 1:3)-(pp(jj, 1:3)-ll*dsign(1.0d0, pp(jj, 1:3)-ll/2.0d0))
    endwhere

    irx(:)=rx(:)/dsqrt(rx(1)**2.+rx(2)**2.+rx(3)**2.)    
    tvi=pv(ii, 1)*irx(1)+pv(ii, 2)*irx(2)+pv(ii, 3)*irx(3)
    tvj=pv(jj, 1)*irx(1)+pv(jj, 2)*irx(2)+pv(jj, 3)*irx(3)
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

                    where(abs(rx(1:3))>ll/2.0d0)
                        rx(1:3)=pp(i, 1:3)-(pp(j, 1:3)-ll*dsign(1.0d0, pp(j, 1:3)-ll/2.0d0))
                    endwhere

                    vx(:)=pv(i, :)-pv(j, :)
                    
                    aij=(vx(1)**2.+vx(2)**2.+vx(3)**2.)-(rin(i)+rin(j))**2.
                    bij=(rx(1)*vx(1)+rx(2)*vx(2)+rx(3)*vx(3))-(ra(i)+ra(j))*(rin(i)+rin(j))
                    cij=(rx(1)**2.+rx(2)**2.+rx(3)**2.)-(ra(i)+ra(j))**2.
                    dij=bij**2.-aij*cij
                    if (dij>=0.0d0) then                    
                        tij=(-bij-dsqrt(dij))/aij
                        if (tij>=0.0d0) then
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
                endif
            enddo
        endif
    enddo

!return

    do i=1, pn
        if (i==ii.or.i==jj) then
            botime(i)=bigmax
            do j=-3, -3
                if (-pv(i, -j)+rin(i)/=0.0d0) then
                    tij=(pp(i, -j)-ra(i))/(-pv(i, -j)+rin(i))
                    if (tij>=0.0d0) then
                        if (tij<botime(i)) then
                            botime(i)=tij
                            bounder(i)=j
                        endif
                    endif
                endif
            enddo
            do j=3, 3
                if (pv(i, j)+rin(i)/=0.0d0) then
                    tij=(le(j)-(pp(i, j)+ra(i)))/(pv(i, j)+rin(i))
                    if (tij>=0.0d0) then
                        if (tij<botime(i)) then
                            botime(i)=tij
                            bounder(i)=j
                        endif
                    endif
                endif
            enddo            
        endif
    enddo
end subroutine periodic_hardsphere_event
!}}}

end module module_common
