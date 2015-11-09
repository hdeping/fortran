program main
    integer, parameter    ::                pn=1000
    real*8, parameter     ::                pi=3.1415926
    integer               ::                i, si, it, ac, k1
    real*8                ::                as, packphi, ll, vl, k2
    real*8                ::                crp, lcrp
    
    type :: node
        integer               ::                se
        type(node),pointer    ::                next
    end type node
    
    type(node), pointer                    ::    head, ver, tver
    
    real*8, dimension(1:pn, 2)             :: pp
    real*8, dimension(1:pn)                :: ra
    real*8, dimension(:, :), allocatable   :: tpp
    real*8, dimension(:), allocatable      :: tra
    
    
    
    call random_seed()
    
    crp   = 0.75d0
    lcrp  = 0.749d0
    ll    = 11.1d0
    
    open(1, file='ini_0.70_11.1.txt')
    read(1, '(f25.16)') k2
    do i=1, pn
        read(1, '(i5, 3f25.16)') k1, ra(i), pp(i, :)
    enddo
    close(1)
    
    ra=ra-0.000001d0
    
    as=0.0d0
    do i=1, pn
        as=as+pi*(ra(i)**2.)
    enddo
    packphi=as/(ll*ll)
    print*, packphi
    
    nullify(head)
    do i=1, pn
        if (i==1) then
            allocate(head)
            head%se=i
            nullify(head%next)
        else
            ver=>head
            allocate(head)
            head%se=i
            head%next=>ver
        endif
    enddo
    
    si=pn
    
    do
    
    100    call random_number(r)
        it=dble(si)*r+1
        if (it==1) then
            ii=head%se
            pa_packphi=pi*(ra(ii)**2.)/(ll*ll)
            tpackphi=packphi-pa_packphi
            if (tpackphi>=crp) then
                ver=>head%next
                deallocate(head)
                head=>ver
                packphi=tpackphi
                si=si-1
            elseif (tpackphi>lcrp.and.tpackphi<crp) then
                ver=>head%next
                deallocate(head)
                head=>ver
                packphi=tpackphi
                exit            
            elseif (tpackphi<=lcrp) then
                go to 100
            endif
        else
            ver=>head
            do i=2, it    
                tver=>ver
                ver=>ver%next
            enddo
            ii=ver%se
            pa_packphi=pi*(ra(ii)**2.)/(ll*ll)
            tpackphi=packphi-pa_packphi        
            if (tpackphi>=crp) then
                if (.not.associated(ver%next)) then
                    nullify(tver%next)
                    deallocate(ver)
                else
                    tver%next=>ver%next
                    deallocate(ver)
                endif
                packphi=tpackphi
                si=si-1
            elseif (tpackphi>lcrp.and.tpackphi<crp) then
                if (.not.associated(ver%next)) then
                    nullify(tver%next)
                    deallocate(ver)
                else
                    tver%next=>ver%next
                    deallocate(ver)
                endif
                packphi=tpackphi
                exit                        
            elseif (tpackphi<=lcrp) then
                go to 100
            endif
        endif
    
    enddo
    
    ac=0
    ver=>head
    do
        if (.not.associated(ver%next)) then
            ac=ac+1
            exit
        else
            ac=ac+1
            ver=>ver%next
        endif
    enddo
    
    allocate(tpp(ac, 2))
    allocate(tra(ac))
    
    ac=0
    ver=>head
    do
    
        if (.not.associated(ver%next)) then
            ac=ac+1
            tpp(ac, :)=pp(ver%se, :)
            tra(ac)=ra(ver%se)
            exit
        else
            ac=ac+1
            tpp(ac, :)=pp(ver%se, :)
            tra(ac)=ra(ver%se)
            ver=>ver%next
        endif
    enddo
    
    as=0.0d0
    do i=1, ac
        as=as+pi*(tra(i)**2.)/(ll*ll)
    enddo
    
    print*, as
    
    open(1, file='ini_.txt')
    do i=ac, 1, -1
        write(1, '(i5, 4f25.16)') ac-i+1, tra(i), tpp(i, :)
    enddo
    close(1)

end program main
