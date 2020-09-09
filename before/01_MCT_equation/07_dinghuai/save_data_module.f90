module save_data_module
use parameter_module

contains

!--***************************************--!

subroutine binary_doc(istep, state)
integer                    ::                istep
integer                    ::                state
character(len=10)        ::                filename

if (state==1) then
    write(filename, '(i2.2, i2.2, a5)') istep, senode, 'ke_me'
    open(unit=senode+1, file='data/'//trim(filename), form='binary', status='new')
        do j=1, sn
            if (mod(j, stn)==0) then
                write(senode+1) ke_me(inode(senode, 1):inode(senode, 2), j)
            endif
        enddo
    close(senode+1)
else
    write(filename, '(i2.2, i2.2, a5)') istep-1, senode, 'ke_me'
    open(unit=senode+1, file='data/'//trim(filename), form='binary', status='old')
        do j=1, sn
            if (mod(j, stn)==0) then
                read(senode+1) ke_me(inode(senode, 1):inode(senode, 2), j/stn)
            endif
        enddo
    close(senode+1)
endif

end subroutine binary_doc

!--***************************************--!

subroutine out_data(istep)
integer                                        ::                i, j, istep
character(len=10)                            ::                filename_1, filename_2
real*8, dimension(:), allocatable        ::                temphi                

    allocate(temphi(inode(senode, 1):inode(senode, 2)))
    write(filename_1, '(i2.2, i2.2, a3)') istep, senode, 'phi' 
    open(unit=senode+1, file='data/'//trim(filename_1), form='binary')
    
if (istep==0) then
    do i=inode(senode, 1), inode(senode, 2)
        write(filename_2, '(i3.3, a4)')  i, '.txt'
        open(2*node+1+i, file='result/'//trim(filename_2), status='new')
    enddo
    do j=0, sn
        read(senode+1) temphi(inode(senode, 1):inode(senode, 2)) 
        if (mod(j, 1)==0) then
            do i=inode(senode, 1), inode(senode, 2)
                write(2*node+1+i, '(e22.12e3, f25.16)') det_time*dble(j), temphi(i)
            enddo
        endif
    enddo
    do i=inode(senode, 1), inode(senode, 2)
        close(2*node+1+i)
    enddo
elseif (istep/=0) then
    do i=inode(senode, 1), inode(senode, 2)
        write(filename_2, '(i3.3, a4)') i, '.txt'
        open(2*node+1+i, file='result/'//trim(filename_2), status='old', position='append')        
    enddo
    do j=0, sn
        read(senode+1) temphi(inode(senode, 1):inode(senode, 2))
        if (j>sn/stn.and.mod(j, 1)==0) then
            do i=inode(senode, 1), inode(senode, 2)
                write(2*node+1+i, '(e22.12e3, f25.16)') det_time*dble(j), temphi(i)
            enddo
        endif
    enddo
    do i=inode(senode, 1), inode(senode, 2)
        close(2*node+1+i)
    enddo
endif

    deallocate(temphi)
    close(senode+1)

end subroutine out_data

!--***************************************--!

end module save_data_module
