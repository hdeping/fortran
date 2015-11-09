module module_common
    implicit none
    integer,parameter        :: n     = 4
    integer,parameter        :: m     = int(1E8)
    integer,parameter        :: fre   = 100
    real(8),parameter        :: delta = 1E-2
    real(8),parameter        :: error = 1E-6
    integer                  :: i
    !integer                  :: j
    !integer                  :: k
    integer                  :: times
    real(8)                  :: t1
    real(8)                  :: t2
    real(8)                  :: t3
    real(8)                  :: t4
    real(8)                  :: lambda
    real(8)                  :: x(n)
    real(8)                  :: deltax(n)
    real(8)                  :: a(n,n)
    real(8)                  :: b(n)
    character(20)            :: filename

end module module_common
