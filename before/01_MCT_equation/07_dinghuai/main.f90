program main
use parameter_module
use initial_module
use evolution_module

character(len=30)    ::        cmdfolder
logical                ::        istatus
integer                ::        sys
    
    call mpi_init(ierror)
    call mpi_comm_rank(mpi_comm_world, rank, ierror)
    call mpi_comm_size(mpi_comm_world, size, ierror)

    sys=1            !1 is windows
!    sys=0            !0 is linux    
    senode=rank
    if (senode==1) then
        open(500, file='dphi.txt')
    endif
    if (senode==0) then
        open(200, file='count.txt')
        open(400, file='err.txt')

        inquire(file='data', exist=istatus)
        if (.not.istatus) then
            cmdfolder='mkdir data'
            call system(trim(cmdfolder))
        else
            select case (sys)
            case(0)
                cmdfolder='rm -rf data'; call system(trim(cmdfolder))                !linux
                cmdfolder='mkdir data'; call system(trim(cmdfolder))                !linux
            case(1)
                cmdfolder='del data /q'; call system(trim(cmdfolder))                !windows
            case default
                open(300, file='remind.txt')
                    write(300, '(a15)') 'system is wrong'
                    call mpi_abort(mpi_comm_world, 99, ierror)
                close(300)
            end select
        endif
        inquire(file='result', exist=istatus)
        if (.not.istatus) then
            cmdfolder='mkdir result'
            call system(trim(cmdfolder))
        else
            select case (sys)
            case(0)
                cmdfolder='rm -rf result'; call system(trim(cmdfolder))                !linux
                cmdfolder='mkdir result'; call system(trim(cmdfolder))                !linux    
            case(1)
                cmdfolder='del result /q'; call system(trim(cmdfolder))                !windows
            end select
        endif
    endif    
    
    call mpi_barrier(mpi_comm_world, ierror)

    call initial_parameter()
    call initial_phi()
    call evolution()

    call mpi_barrier(mpi_comm_world, ierror)

    call mpi_finalize(ierror)

end program main
