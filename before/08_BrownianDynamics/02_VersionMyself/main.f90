program main
    use module_common
    implicit none

    call random_seed()

    ! get the initial data

    call getInitialParameter()

    ! evolution the process

    call evolution()
end program main
