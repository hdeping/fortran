program main
    implicit none
    integer,parameter        :: n     = 100
    integer,parameter        :: m     = int(1E8)
    real(8),parameter        :: delta = 1E-2
    real(8),parameter        :: x1    = 1E-2
    integer                  :: i
    integer                  :: j
    integer                  :: k
    
    interface cal
        module procedure ical,dcal
    end interface
    function ical(n)
        integer    n

    end function ical


end program main
