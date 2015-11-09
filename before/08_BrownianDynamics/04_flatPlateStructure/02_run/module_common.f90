module module_common
! some parameters{{{
    integer, parameter          :: pn     = 949
    integer, parameter          :: base   = 20000
    integer, parameter          :: ln2    = 10
    real*8, parameter           :: ll     = 18.0d0
    real*8, parameter           :: vl     = 3.0d0
    real*8, parameter           :: dl     = ll/dble(ln2)
    real*8, parameter           :: bigmax = 1.0e8
    real*8, dimension(3)        :: le
    real*8, dimension(pn, 3)    :: pp, pv, pp0, ppa
    integer, dimension(pn, 2)   :: pl
    real*8, dimension(pn, 3)    :: ki
    integer, dimension(pn)      :: partner, bounder
    real*8, dimension(pn)       :: cotime, botime
    real*8, dimension(pn)       :: ra
    integer                     :: motion_step
    
    
    type :: node
        integer                    ::            se
        type(node), pointer        ::            next
    end type node
    
    type :: delay_structure    
        type(node), pointer        ::            pt
    end type delay_structure
    
    type(delay_structure), dimension(0:ln2-1, 0:ln2-1)      ::        lhead
    type(delay_structure), dimension(pn)                    ::        phead
    
    real*8, parameter              ::                ddr=0.2d0, ddr2=0.1d0
    real*8, dimension(2)           ::                dpc
    real*8                         ::                det_time
    real*8                         ::                fri, v0, tau, dr, clength, d
!}}}
    contains
    ! initial some processors
!subroutine initial_parameter{{{
    !--*****************************************--!
subroutine initial_parameter()
    real*8                ::                packphi, ra1, ra2
    integer               ::                i, k1, k2
    character(len = 20)   ::                filename
    real*8, parameter     ::                pi=3.1415926
    
    !--*******************--!
    !d diffusion
    !v0  activity
    !tau time 
    d        = 0.2d0
    v0       = 10.0d0/dsqrt(3.0d0)
    tau      = 0.05
    dr       = 1.0/(2.0*tau)
    fri      = 2.0d0*dr
    det_time = 1.0e-4
    
    !--*******************--!
    
    filename = 'init_55_18_18_3'
    open(1, file = filename )
    !read(1, *) packphi, ra1, ra2
    do i=1, pn
        read(1, '(i6, 4f25.16)') k1, ra(i), pp(i, :)
    enddo
    close(1)
    
    
    do i=1, (pn-1)/2
        do j=1, 3
            call gaussian_noise(ki(i, j))
        enddo
    enddo
    
    do i=(pn-1)/2+2, pn
        do j=1, 3
            ki(i, j)=-ki(i-(pn-1)/2-1, j)
        enddo
    enddo
    
    ki((pn-1)/2+1, :)=0.0d0
    
    call create_list()
    
    call create_neighbor_list()
    
    motion_step=0
    ppa=pp
    pp0=pp
    
    end subroutine initial_parameter
!}}}
!   data writing
!subroutine write_data{{{
subroutine write_data()
    integer                  :: i
    character(len=15)        :: filename
    
        if (motion_step>=base) then
        
            if (mod(motion_step, 1)==0.and.motion_step-base<=18000) then
                write(filename, '(i6.6, a4)') motion_step-base, '.txt'
                open(1, file='1/'//filename)
                do i=1, pn
                    write(1, '(i5, 3f25.16)') i, ppa(i, 1:2), pp(i, 3)
                enddo
                close(1)
            endif
            if (mod(motion_step, 100)==0.and.motion_step-base<=1800000) then
                write(filename, '(i6.6, a4)') (motion_step-base)/100, '.txt'
                open(1, file='100/'//filename)
                do i=1, pn
                    write(1, '(i5, 3f25.16)') i, ppa(i, 1:2), pp(i, 3)
                enddo
                close(1)
                if (motion_step-base==1800000) then
                    stop
                endif
            endif
    
            if (mod(motion_step, 10000)==0) then
                call exam_dis()
            endif
        
        endif
    
    end subroutine write_data
!}}}
!subroutine exam_dis{{{
subroutine exam_dis()
    integer                    ::            i, j
    real*8, dimension(3)    ::            lx
    real*8                    ::            dis
    
    open(10, file='exam_dis.txt')
    
    do i=1, pn-1
        do j=i+1, pn
            lx=abs(pp(i, :)-pp(j, :))
            where(lx>=ll/2.0d0)
                lx=ll-lx
            endwhere
            dis=dsqrt(lx(1)**2.+lx(2)**2.+lx(3)**2.)
            if (dis<ra(i)+ra(j)) then
                write(10, *) dis, ra(i)+ra(j)
            endif
        enddo
    enddo
    
    write(10, '(a7)') 'perfect'
    close(10)
    
    end subroutine exam_dis
!}}}
!   evolution 
!subroutine evolution{{{
subroutine evolution()
    real*8, dimension(3)        ::            rx
    integer, dimension(1)        ::            it
    integer                        ::            i
    real*8                        ::            netdis
            
    motion_step=0
            
    do 
            
        motion_step=motion_step+1
    
        !if (mod(motion_step, 100)==0) then
        !    print*, motion_step
        !call dis_ij()
        !print*, 'top', vl-maxval(pp(:, 3)), minval(pp(:, 3))
        !endif
    
        call fctive_veolocity_extract()
    
        call new_collision_list()
    
        call once_collision_evolution()
    
        call write_data()
            
        do i=1, pn    
            
            rx(:)=abs(pp(i, :)-pp0(i, :))
    
            where(rx>=ll/2.0d0)
                rx=ll-rx
            endwhere
    
            netdis=sqrt(rx(1)**2.+rx(2)**2.+rx(3)**2.)
            it=minloc(dpc)
            
            if (netdis>dpc(it(1))) then
                dpc(it(1))=netdis
            endif
            
        enddo
    enddo    
    end subroutine evolution
!}}}
! process of the lists
!subroutine create_list{{{
subroutine create_list()
integer                        ::            i, j, k, ikind
integer, dimension(2)        ::            pi
type(node), pointer            ::            ver


do i=0, ln2-1
    do j=0, ln2-1
        nullify(lhead(i, j)%pt)
    enddo
enddo

do i=1, pn

    pi(:)=int(pp(i, 1:2)/dl)
    pl(i, :)=pi(:)

    if (.not.associated(lhead(pi(1), pi(2))%pt)) then
        allocate(lhead(pi(1), pi(2))%pt)
        lhead(pi(1), pi(2))%pt%se=i
        nullify(lhead(pi(1), pi(2))%pt%next)
    else
        ver=>lhead(pi(1), pi(2))%pt
        allocate(lhead(pi(1), pi(2))%pt)
        lhead(pi(1), pi(2))%pt%se=i
        lhead(pi(1), pi(2))%pt%next=>ver
    endif
enddo


end subroutine create_list
!}}}
!subroutine create_neighbor_list{{{
subroutine create_neighbor_list()
integer                            ::                i, ikind, jkind
integer, dimension(2)            ::                pi, npi
integer                            ::                ii, jj
real*8, dimension(2)            ::                per
type(node), pointer                ::                ver, tver
real*8                            ::                dis    

    do i=1, pn
        
        pi(:)=pl(i, :)    

        do ii=-1, 1
            do jj=-1, 1

                npi(1)=pi(1)+ii; npi(2)=pi(2)+jj
                per(:)=0.0d0

                where(npi==ln2)
                    npi=0
                    per=ll
                elsewhere(npi(:)==-1)
                    npi=ln2-1
                    per=-ll
                endwhere
                    
                if (.not.associated(lhead(npi(1), npi(2))%pt)) then
                    cycle
                else
                    ver=>lhead(npi(1), npi(2))%pt
                    do
                        j=ver%se
                        if (j>i) then
                            dis=sqrt((pp(i, 1)-(pp(j, 1)+per(1)))**2.+(pp(i, 2)-(pp(j, 2)+per(2)))**2.+(pp(i, 3)-pp(j, 3))**2.)
                            if (dis<=ra(i)+ra(j)+ddr) then
                                if (.not.associated(phead(i)%pt)) then
                                    allocate(phead(i)%pt)
                                    phead(i)%pt%se=j
                                    nullify(phead(i)%pt%next)
                                else
                                    tver=>phead(i)%pt
                                    allocate(phead(i)%pt)
                                    phead(i)%pt%se=j
                                    phead(i)%pt%next=>tver
                                endif
                                if (.not.associated(phead(j)%pt)) then
                                    allocate(phead(j)%pt)
                                    phead(j)%pt%se=i
                                    nullify(phead(j)%pt%next)
                                else
                                    tver=>phead(j)%pt
                                    allocate(phead(j)%pt)
                                    phead(j)%pt%se=i
                                    phead(j)%pt%next=>tver
                                endif
                            endif                                                                
                        endif
                        if (.not.associated(ver%next)) then
                            exit
                        else    
                            ver=>ver%next
                        endif
                    enddo
                endif
            enddo
        enddo

    enddo

    nullify(ver); nullify(tver)

end subroutine create_neighbor_list
!}}}
!subroutine updata_list{{{
subroutine updata_list()
integer                        ::            i
integer, dimension(2)        ::            ti, nti

    do i=1, pn
        
        nti(:)=int(pp(i, 1:2)/dl)
        ti(:)=pl(i, :)

        if (nti(1)==ti(1).and.nti(2)==ti(2)) then
            cycle
        else
            call list_change(i, ti, nti)
            pl(i, :)=nti(:)
        endif
        
    enddo

end subroutine updata_list
!}}}
!subroutine list_change{{{
subroutine list_change(i, ti, nti)
integer                            ::            i
integer, dimension(2)            ::            ti, nti
type(node), pointer                ::            ver, tver

ver=>lhead(ti(1), ti(2))%pt
tver=>lhead(ti(1), ti(2))%pt

if (ver%se==i) then
    if (.not.associated(ver%next)) then
        deallocate(ver)
        nullify(lhead(ti(1), ti(2))%pt)
    else
        lhead(ti(1), ti(2))%pt=>ver%next
        deallocate(ver)
    endif
else
    do
        ver=>ver%next
        if (i==ver%se) then
            if (.not.associated(ver%next)) then
                deallocate(ver)
                nullify(tver%next)
            else
                tver%next=>ver%next
                deallocate(ver)
            endif
            exit
        else
            tver=>ver
        endif
    enddo
endif

if (.not.associated(lhead(nti(1), nti(2))%pt)) then
    allocate(lhead(nti(1), nti(2))%pt)
    lhead(nti(1), nti(2))%pt%se=i
    nullify(lhead(nti(1), nti(2))%pt%next)
else
    ver=>lhead(nti(1), nti(2))%pt
    allocate(lhead(nti(1), nti(2))%pt)
    lhead(nti(1), nti(2))%pt%se=i
    lhead(nti(1), nti(2))%pt%next=>ver
endif

nullify(tver)
nullify(ver)

end subroutine list_change
!}}}
! pattern module
!subroutine fctive_veolocity_extract{{{
subroutine fctive_veolocity_extract()
integer                    ::            i, j
real*8                    ::            r
real*8, dimension(3)    ::            white_r
real*8                    ::            dk


do i=1, pn
    do j=1, 3
        call gaussian_noise(white_r(j))
        dk=-fri*ki(i, j)*det_time+sqrt(2.0*fri*det_time)*white_r(j)
        ki(i, j)=ki(i, j)+dk
    enddo
enddo

do i=1, pn    
        do j=1, 3
        call gaussian_noise(r)
        pv(i, j)=v0*ki(i, j)+sqrt(2.0*d/det_time)*r   
    enddo
enddo

end subroutine fctive_veolocity_extract
!}}}
!subroutine new_collision_list{{{
subroutine new_collision_list()
real*8                    ::            ac

    ac=sum(dpc(:))    
    
    if (ac<ddr2) then
        call list_collision_time()
    else
        pp0=pp
        dpc=0.0d0
        call delete_list()
        call updata_list()
        call create_neighbor_list_and_collision_time()
    endif

end subroutine new_collision_list
!}}}
!subroutine delete_list{{{
subroutine delete_list()
integer                        ::                i
type(node), pointer            ::                ver

    do i=1, pn
        if (.not.associated(phead(i)%pt)) then
            cycle
        else
            ver=>phead(i)%pt
            do
                if (.not.associated(ver%next)) then
                    deallocate(ver)
                    nullify(phead(i)%pt)
                    exit
                else
                    ver=>ver%next
                    deallocate(phead(i)%pt)
                    phead(i)%pt=>ver
                endif
            enddo
        endif
    enddo

    nullify(ver)

end subroutine delete_list
!}}}
!subroutine create_neighbor_list_and_collision_time{{{
subroutine create_neighbor_list_and_collision_time()
integer                  :: i, j
real*8, dimension(3)     :: rx, vx
real*8                   :: aij, bij, cij, dij, tij, dis
integer, dimension(2)    :: pi, npi
real*8, dimension(2)     :: per
integer                  :: ii, jj
type(node), pointer      :: ver, tver

cotime=bigmax
botime=bigmax

do i=1, pn

    pi(:)=pl(i, :)
    
    do ii=-1, 1
        do jj=-1, 1

                npi(1)=pi(1)+ii; npi(2)=pi(2)+jj
                per(:)=0.0d0
                where(npi==ln2)
                    npi=0
                    per=ll
                elsewhere(npi==-1)
                    npi=ln2-1
                    per=-ll
                endwhere
                
                if (.not.associated(lhead(npi(1), npi(2))%pt)) then
                    cycle
                else

                    ver=>lhead(npi(1), npi(2))%pt

                    do

                        j=ver%se
                        if (j/=i) then
                        rx(1:2)=pp(i, 1:2)-(pp(j, 1:2)+per(:))
                        rx(3)=pp(i, 3)-pp(j, 3)
                        dis=dsqrt(rx(1)**2.+rx(2)**2.+rx(3)**2.)
                        
                        if (dis<=ra(i)+ra(j)+ddr) then
                            
                            vx=pv(i, :)-pv(j, :)
                            bij=rx(1)*vx(1)+rx(2)*vx(2)+rx(3)*vx(3)

                            if (bij<0.0d0) then
                                aij=vx(1)**2.+vx(2)**2.+vx(3)**2.
                                cij=dis**2.-(ra(i)+ra(j))**2.
                                dij=bij**2.-aij*cij
                                
                                if (dij>=0.0d0) then
                                    tij=(-bij-dsqrt(dij))/aij
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
                        
                            if (j>i) then
                                if (.not.associated(phead(i)%pt)) then
                                    allocate(phead(i)%pt)
                                    phead(i)%pt%se=j
                                    nullify(phead(i)%pt%next)
                                else
                                    tver=>phead(i)%pt
                                    allocate(phead(i)%pt)
                                    phead(i)%pt%se=j
                                    phead(i)%pt%next=>tver
                                endif
                                if (.not.associated(phead(j)%pt)) then
                                    allocate(phead(j)%pt)
                                    phead(j)%pt%se=i
                                    nullify(phead(j)%pt%next)
                                else
                                    tver=>phead(j)%pt
                                    allocate(phead(j)%pt)
                                    phead(j)%pt%se=i
                                    phead(j)%pt%next=>tver
                                endif
                            endif                        
                        
                        endif
                        endif

                        if (.not.associated(ver%next)) then
                            exit
                        else
                            ver=>ver%next
                        endif




                    enddo
                endif

        enddo
    enddo

    if (pp(i, 3)-ra(i)<ddr) then
        j=-3
        if (pv(i, 3)<0.0d0) then
            tij=(pp(i, -j)-ra(i))/(-pv(i, -j))
            if (tij<botime(i)) then
                botime(i)=tij
                bounder(i)=j
            endif
        endif
    elseif (vl-(pp(i, 3)+ra(i))<ddr) then
        j=3
        if (pv(i, 3)>0.0d0) then
            tij=(vl-(pp(i, j)+ra(i)))/(pv(i, j))
            if (tij<botime(i)) then
                botime(i)=tij
                bounder(i)=j
            endif
        endif
    endif

enddo

nullify(ver)

end subroutine create_neighbor_list_and_collision_time
!}}}
!subroutine list_collision_time{{{
subroutine list_collision_time()
integer                        ::                    i, j
real*8                        ::                    bij, aij, cij, dij, tij
real*8, dimension(3)        ::                    rx, vx
type(node), pointer            ::                    ver

cotime=bigmax
botime=bigmax
    
    do i=1, pn

        if (.not.associated(phead(i)%pt)) then
            cycle
        else    
            ver=>phead(i)%pt
            do
                j=ver%se
                if (j/=i) then
                rx(:)=pp(i, :)-pp(j, :)
                where(abs(rx)>ll/2.0d0)
                    rx=pp(i, :)-(pp(j, :)-ll*dsign(1.0d0, pp(j, :)-ll/2.0d0))
                endwhere
                
                vx=pv(i, :)-pv(j, :)
                bij=rx(1)*vx(1)+rx(2)*vx(2)+rx(3)*vx(3)
                    if (bij<=0.0d0) then
                        aij=vx(1)**2.+vx(2)**2.+vx(3)**2.
                        cij=(rx(1)**2.+rx(2)**2.+rx(3)**2.)-(ra(i)+ra(j))**2.
                        dij=bij**2.-aij*cij
                        if (dij>=0.0d0) then
                            tij=(-bij-dsqrt(dij))/aij
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
                if (.not.associated(ver%next)) then
                    exit
                else
                    ver=>ver%next
                endif
            enddo
        endif

        if (pp(i, 3)-ra(i)<ddr) then
            j=-3
            if (pv(i, 3)<0.0d0) then
                tij=(pp(i, -j)-ra(i))/(-pv(i, -j))
                if (tij<botime(i)) then
                    botime(i)=tij
                    bounder(i)=j
                endif
            endif
        elseif (vl-(pp(i, 3)+ra(i))<ddr) then
            j=3
            if (pv(i, 3)>0.0d0) then
                tij=(vl-(pp(i, j)+ra(i)))/(pv(i, j))
                if (tij<botime(i)) then
                    botime(i)=tij
                    bounder(i)=j
                endif
            endif
        endif

    enddo

nullify(ver)

end subroutine list_collision_time
!}}}
!subroutine once_collision_evolution{{{
subroutine once_collision_evolution()
integer, dimension(1)                ::                itc, itb
integer                                ::                k, istep, kind_of_collision
real*8                                ::                actime, tij, tijc, tijb
real*8, dimension(3)                ::                ddx

actime=0.0d0
istep=0

do

    istep=istep+1
    itc(:)=minloc(cotime)
    tijc=cotime(itc(1))
    itb(:)=minloc(botime)
    tijb=botime(itb(1))
    
    if (tijc>=tijb) then
        tij=tijb
        kind_of_collision=0
    else
        tij=tijc
        kind_of_collision=1
    endif

    if (tij>=det_time.and.istep==1) then
        
        do k=1, pn
            ddx(:)=pv(k, :)*det_time
            pp(k, :)=pp(k, :)+ddx(:)
            ppa(k, 1:2)=ppa(k, 1:2)+ddx(1:2)
            where(pp(k, :)>ll)
                pp(k, :)=pp(k, :)-ll
            elsewhere(pp(k, :)<=0.0d0)
                pp(k, :)=pp(k, :)+ll
            endwhere
        enddo
            
        exit

    else
        
        actime=actime+tij

        if (actime>=det_time) then
            do k=1, pn
                ddx(:)=pv(k, :)*(det_time-(actime-tij))
                pp(k, :)=pp(k, :)+ddx(:)
                ppa(k, 1:2)=ppa(k, 1:2)+ddx(1:2)
                where(pp(k, :)>ll)
                    pp(k, :)=pp(k, :)-ll
                elsewhere(pp(k, :)<=0.0d0)
                    pp(k, :)=pp(k, :)+ll
                endwhere
            enddo

            exit

        else

            do k=1, pn
                cotime(k)=cotime(k)-tij
                botime(k)=botime(k)-tij
                ddx(:)=pv(k, :)*tij
                pp(k, :)=pp(k, :)+ddx(:)
                ppa(k, 1:2)=ppa(k, 1:2)+ddx(1:2)
                where(pp(k, :)>ll)
                    pp(k, :)=pp(k, :)-ll
                elsewhere(pp(k, :)<=0.0d0)
                    pp(k, :)=pp(k, :)+ll
                endwhere
            enddo            
                    
            if (kind_of_collision==0) then
                call boundary_event(itb(1))
            elseif (kind_of_collision==1) then
                call hardsphere_event(itc(1))
            endif

        endif

    endif

enddo


end subroutine once_collision_evolution
!}}}
!subroutine boundary_event{{{
subroutine boundary_event(ii)
integer                        ::                ii, jj, i, j
type(node), pointer            ::                ver
real*8, dimension(3)        ::                rx, vx
real*8                        ::                aij, bij, cij, dij, tij

    jj=bounder(ii)

    pv(ii, abs(jj))=-pv(ii, abs(jj))
    
    do i=1, pn
        if ((i==ii).or.(partner(i)==ii)) then
            cotime(i)=bigmax        
            if (.not.associated(phead(i)%pt)) then
                cycle
            else

                ver=>phead(i)%pt

                do 
                    j=ver%se
                    if (j/=i) then

                        rx(:)=pp(i, :)-pp(j, :)
                        where(abs(rx)>ll/2.0d0)
                            rx=pp(i, :)-(pp(j, :)-ll*dsign(1.0d0, pp(j, :)-ll/2.0d0))
                        endwhere                    
                        
                        vx=pv(i, :)-pv(j, :)
                        bij=rx(1)*vx(1)+rx(2)*vx(2)+rx(3)*vx(3)
                        if (bij<=0.0d0) then
                            aij=vx(1)**2.+vx(2)**2.+vx(3)**2.
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
                    endif
                    if (.not.associated(ver%next)) then
                        exit    
                    else
                        ver=>ver%next
                    endif                                
                enddo
            endif
        endif
    enddo    

    botime(ii)=bigmax

nullify(ver)

end subroutine boundary_event
!}}}
!subroutine hardsphere_event{{{
subroutine hardsphere_event(ii)
integer                        ::                ii, jj, i, j
real*8                        ::                tij, aij, bij, cij, dij, tvi, tvj
real*8, dimension(3)        ::                rx, vx, irx, vxpi, vxpj, vxvi, vxvj
type(node), pointer            ::                ver

    jj=partner(ii)
    
    rx(:)=pp(ii, :)-pp(jj, :)
    
    where(abs(rx)>ll/2.0d0)
        rx=pp(ii, :)-(pp(jj, :)-ll*dsign(1.0d0, pp(jj, :)-ll/2.0d0))
    endwhere

    irx(:)=rx(:)/dsqrt(rx(1)**2.+rx(2)**2.+rx(3)**2.)    
    tvi=pv(ii, 1)*irx(1)+pv(ii, 2)*irx(2)+pv(ii, 3)*irx(3)
    tvj=pv(jj, 1)*irx(1)+pv(jj, 2)*irx(2)+pv(jj, 3)*irx(3)
    vxpi(:)=tvi*irx(:); vxpj(:)=tvj*irx(:)
    vxvi(:)=pv(ii, :)-vxpi(:); vxvj(:)=pv(jj, :)-vxpj(:)
    
    pv(ii, :)=vxpj(:)+vxvi(:)
    pv(jj, :)=vxpi(:)+vxvj(:)

    do i=1, pn
        if ((i==ii).or.(i==jj).or.(partner(i)==ii).or.(partner(i)==jj)) then
            cotime(i)=bigmax
            if (.not.associated(phead(i)%pt)) then
                cycle
            else
                ver=>phead(i)%pt
                do 
                    j=ver%se
                    if (j/=i) then
                        rx(:)=pp(i, :)-pp(j, :)
    
                        where(abs(rx)>ll/2.0d0)
                            rx=pp(i, :)-(pp(j, :)-ll*dsign(1.0d0, pp(j, :)-ll/2.0d0))
                        endwhere                    
                        
                        vx=pv(i, :)-pv(j, :)

                        bij=rx(1)*vx(1)+rx(2)*vx(2)+rx(3)*vx(3)
                        if (bij<=0.0d0) then
                            aij=vx(1)**2.+vx(2)**2.+vx(3)**2.
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
                    endif
                    if (.not.associated(ver%next)) then
                        exit    
                    else
                        ver=>ver%next
                    endif                                
                enddo
            endif            
        endif
    enddo


do k=1, 2
    if (k==1) then
        i=ii
    elseif (k==2) then
        i=jj
    endif
    botime(i)=bigmax
    if (pp(i, 3)-ra(i)<ddr) then
        j=-3
        if (pv(i, 3)<0.0d0) then
            tij=(pp(i, -j)-ra(i))/(-pv(i, -j))
            if (tij<botime(i)) then
                botime(i)=tij
                bounder(i)=j
            endif
        endif
    elseif (vl-(pp(i, 3)+ra(i))<ddr) then
        j=3
        if (pv(i, 3)>0.0d0) then
            tij=(vl-(pp(i, j)+ra(i)))/(pv(i, j))
            if (tij<botime(i)) then
                botime(i)=tij
                bounder(i)=j
            endif
        endif
    endif    
enddo

end subroutine hardsphere_event
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
end module module_common
