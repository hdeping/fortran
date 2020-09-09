module evolution_module
use parameter_module
use function_module


contains


!--*********************************************--!

subroutine initial_phie()
integer                    ::                l, i

    istep=0

    phie(:, 0)=1.0d0
    phie(0, 0)=0.0d0
    
    do l=1, n-1
        ke_mee(l, 0)=kernele(l, 0)
    enddo
    ke_mee(0, 0)=0.0d0

    do i=1, n-1
        dphie(i, 0)=-(au(i)*phie(i, 0)*det_time)/(1.0d0+bu(i)*ke_mee(i, 0)*det_time)
    enddo

    dphie(0, 0)=0.0d0
    phia(:)=0.0d0    

end subroutine initial_phie

!--*********************************************--!

subroutine eq_evolution()
integer                                        ::                        i, l, k
integer                                        ::                        istep, istart
real*8, dimension(:, :), allocatable        ::                        tphie, tke_mee

allocate(tphie(0:n-1, 1:sn/2))
allocate(tke_mee(0:n-1, 1:sn/2-1))

do istep=0, ini

    if (istep==0) then
        istart=1
    else
        istart=sn/stn
        det_time=det_time*dble(stn)
    endif

    do i=istart, sn

        phie(:, i)=phie(:, i-1)+dphie(:, i-1)

        phie(0, i)=0.0

        do l=1, n-1
            ke_mee(l, i)=kernele(l, i)
        enddo
        
        do l=1, n-1
            dphie(l, i)=(-aphief(l, i)-bu(l)*memorye(l, i))*det_time/(1.0d0+bu(l)*ke_mee(l, 0)*det_time)
        enddo
        
        if (istep==0) then
            do l=1, n-1
                phia(l)=phia(l)+(phie(l, i)**2.)*det_time
            enddo
        else
            if (i/=istart) then
                do l=1, n-1
                    phia(l)=phia(l)+(phie(l, i)**2.)*det_time
                enddo                
            endif
        endif

    enddo

    do i=1, n-1
        do j=istart, sn
            write(i, '(3e25.16e3)') det_time*dble(j), phie(i, j), phia(i)
        enddo
    enddo

    do i=1, sn/2
        tphie(:, i)=phie(:, 2*i)
        if (i/=sn/2) then
            tke_mee(:, i)=ke_mee(:, 2*i)
        endif
    enddo

    phie(:, 1:sn/2)=tphie(:, :)
    ke_mee(:, 1:sn/2-1)=tke_mee(:, :)

    do i=0, sn/2-1
        dphie(:, i)=phie(:, i+1)-phie(:, i)
    enddo

enddo

deallocate(tphie)
deallocate(tke_mee)
do i=1, n-1
    close(i)
enddo

do k=1, n-1
    sk_neq(k)=sk(k)+0.5*df*(dble(k)**2.)*(1.0-sk(k))*phia(k)
enddo

open(1, file='skneq.txt')
do i=1, n-1
    write(1, '(2f25.16)') dble(i)*dk, sk_neq(i)
enddo
close(1)

end subroutine eq_evolution

!--**********************************************--!


end module evolution_module
