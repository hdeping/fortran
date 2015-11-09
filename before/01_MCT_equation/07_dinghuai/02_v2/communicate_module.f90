module communicate_module
use parameter_module

contains

!--******************************************--!
!subroutine comm_initial_parameter{{{
subroutine comm_initial_parameter()
integer                        ::                i
integer, parameter             ::                nn=(4*n)*8
integer(1), dimension(nn)      ::                buff
integer                        ::                position, size    

    if (senode==0) then
        position=0
        call mpi_pack(sk(0), n, mpi_double_precision, buff, nn, position, mpi_comm_world, ierror)
        call mpi_pack(ck(0), n, mpi_double_precision, buff, nn, position, mpi_comm_world, ierror)
        call mpi_pack(ak(0), n, mpi_double_precision, buff, nn, position, mpi_comm_world, ierror)
        call mpi_pack(gammak(0), n, mpi_double_precision, buff, nn, position, mpi_comm_world, ierror)
        do i=1, node-1
            call mpi_send(buff, position, mpi_packed, i*thread, 0, mpi_comm_world, ierror)
        enddo
    else
        call mpi_recv(buff, nn/8, mpi_double_precision, 0, 0, mpi_comm_world, status, ierror)
        position=0
        call mpi_unpack(buff, nn, position, sk(0), n, mpi_double_precision, mpi_comm_world, ierror)
        call mpi_unpack(buff, nn, position, ck(0), n, mpi_double_precision, mpi_comm_world, ierror)
        call mpi_unpack(buff, nn, position, ak(0), n, mpi_double_precision, mpi_comm_world, ierror)
        call mpi_unpack(buff, nn, position, gammak(0), n, mpi_double_precision, mpi_comm_world, ierror)
    endif

end subroutine comm_initial_parameter
!}}}
!subroutine comm_phi{{{
subroutine comm_phi()
integer                        ::                i, j

do i=0, node-1
    if (senode==i) then
        do j=0, node-1
            if (j/=i) then
                call mpi_send(phi(inode(i, 1)), nnode, mpi_double_precision, j*thread, i, mpi_comm_world, ierror)
            endif
        enddo
    else
        call mpi_recv(phi(inode(i, 1)), nnode, mpi_double_precision, i*thread, i, mpi_comm_world, status, ierror)
    endif
enddo

end subroutine comm_phi
!}}}
end module communicate_module
