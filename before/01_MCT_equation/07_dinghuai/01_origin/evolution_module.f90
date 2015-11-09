module evolution_module
use parameter_module
use initial_module
use communicate_module
use save_data_module

contains

!--******************************************--!

subroutine initial_phi()
integer                            ::                i, j, k
real*8, dimension(2, 2)            ::                ta, tc
character(len=10)                ::                filename

    phi(0:n-1)=1.0
    
    do l=inode(senode, 1), inode(senode, 2)
        ke_me(l, 0)=kernel(l)
    enddo 

    do i=inode(senode, 1), inode(senode, 2)
        dphi(i, 0)=-ak(i)*gammak(i)*phi(i)*det_time             !-ak(i)*gammak(i)*phi(i)*det_time/(1.0+(gammak(i)/ak(i))*ke_me(i, 0)*det_time)
    enddo    

    dphi(inode(senode, 1):inode(senode, 2), 0)=dphi(inode(senode, 1):inode(senode, 2), 0)

    write(filename, '(i2.2, i2.2, a3)') istep, senode, 'phi'
    open(unit=senode+1, file='data/'//trim(filename), form='binary')
    write(senode+1) phi(inode(senode, 1):inode(senode, 2))

end subroutine initial_phi

!--******************************************--!

subroutine evolution()
integer                                        ::                    i, j, k, l
real*8, dimension(:), allocatable            ::                    tem_phi
real*8                                        ::                    time
character(len=8)                            ::                    time1
character(len=10)                            ::                    time2
character(len=5)                            ::                    time3
integer(4), dimension(8)                    ::                    time4

    allocate(tem_phi(inode(senode, 1):inode(senode, 2)))

    istep=0

    do i=1, sn
    
    if (mod(i, wstep)==0.and.senode==0) then
        call date_and_time(time1, time2, time3, time4)
        write(200, '(i10, i10, i10, a2, i2, a2, i2, a2, i2, a2, i2, a2, i2, a2, i3, a2)') istep, i/wstep, &
                    time4(1), '年', time4(2), '月', time4(3), '日', time4(5), '时', time4(6), '分', time4(7), '秒', time4(8), '毫'
    endif    
        
        phi(inode(senode, 1):inode(senode, 2))=phi(inode(senode, 1):inode(senode, 2))+dphi(inode(senode, 1):inode(senode, 2), i-1) 
        
        write(senode+1) phi(inode(senode, 1):inode(senode, 2))    

        call comm_phi()
        
        do l=inode(senode, 1), inode(senode, 2)
            ke_me(l, i)=kernel(l)
        enddo 

        do l=inode(senode, 1), inode(senode, 2)
            dphi(l, i)=(-ak(l)*gammak(l)*phi(l)-memory(i, l))*det_time     !/(1.0+(gammak(l)/ak(l))*ke_me(l, 0)*det_time)
        enddo            

    enddo

    tem_phi(inode(senode, 1):inode(senode, 2))=phi(inode(senode, 1):inode(senode, 2))

    close(senode+1)

    call binary_doc(istep, 1)

    call out_data(istep)

    phi(:)=0.0; dphi(inode(senode, 1):inode(senode, 2), :)=0.0; ke_me(inode(senode, 1):inode(senode, 2), 1:sn)=0.0

    do istep=1, ini
        
        det_time=det_time*dble(stn)

        call binary_doc(istep, -1)

        call get_dphi(istep)

        phi(inode(senode, 1):inode(senode, 2))=tem_phi(inode(senode, 1):inode(senode, 2))
        
        do i=sn/stn, sn-1            
        print*, istep, i
            if (mod(i, wstep)==0.and.senode==0) then
                call date_and_time(time1, time2, time3, time4)
                write(200, '(i10, i10, i10, a2, i2, a2, i2, a2, i2, a2, i2, a2, i2, a2, i3, a2)') istep, i/wstep, &
                            time4(1), '年', time4(2), '月', time4(3), '日', time4(5), '时', time4(6), '分', time4(7), '秒', time4(8), '毫'
            endif

            call comm_phi()
            
            if (istep==ini.and.i==sn-1.and.senode==0) then
                open(154, file='k.txt')
                do ii=0, n-1
                    write(154, '(2f18.8)') dble(ii)*dk, phi(ii)
                enddo
                close(154)
            endif

            do l=inode(senode, 1), inode(senode, 2)
                ke_me(l, i)=kernel(l)
            enddo             

            do l=inode(senode, 1), inode(senode, 2)            
                dphi(l, i)=(-ak(l)*gammak(l)*phi(l)-memory(i, l))*det_time     !/(1.0+(gammak(l)/ak(l))*ke_me(l, 0)*det_time)                
            enddo                

            phi(inode(senode, 1):inode(senode, 2))=phi(inode(senode, 1):inode(senode, 2))+dphi(inode(senode, 1):inode(senode, 2), i)
            write(node+senode+1) phi(inode(senode, 1):inode(senode, 2))    

        enddo

        tem_phi(inode(senode, 1):inode(senode, 2))=phi(inode(senode, 1):inode(senode, 2))        

        close(node+senode+1)

        call binary_doc(istep, 1)

        call out_data(istep)

        phi(:)=0.0; dphi(inode(senode, 1):inode(senode, 2), :)=0.0; ke_me(inode(senode, 1):inode(senode, 2), 1:sn)=0.0

    enddo

    deallocate(tem_phi)

end subroutine evolution

!--******************************************--!

subroutine get_dphi(istep)
integer                                        ::                i, istep
real*8, dimension(:), allocatable        ::                temphi, f_temphi, b_temphi    
character(len=10)                            ::                filename

allocate(temphi(inode(senode, 1):inode(senode, 2)))
allocate(f_temphi(inode(senode, 1):inode(senode, 2)))
allocate(b_temphi(inode(senode, 1):inode(senode, 2)))

write(filename, '(i2.2, i2.2, a3)') istep-1, senode, 'phi'
open(unit=senode+1, file='data/'//trim(filename), form='binary')
write(filename, '(i2.2, i2.2, a3)') istep, senode, 'phi'
open(unit=node+senode+1, file='data/'//trim(filename), form='binary') 

do i=0, sn
    read(senode+1) temphi(inode(senode, 1):inode(senode, 2))
    if (mod(i, stn)==0) then
        f_temphi=temphi
        write(node+senode+1) f_temphi(:)
        if (i/=0) then
            dphi(:, i/stn-1)=f_temphi(:)-b_temphi(:)
        endif
        b_temphi=f_temphi
    endif
enddo

    close(senode+1)

deallocate(temphi)
deallocate(f_temphi)
deallocate(b_temphi)

end subroutine get_dphi

!--*****************************************--!

end module evolution_module
