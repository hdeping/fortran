module module_common
    implicit none
    integer,parameter            :: n     = 100
    real(8),parameter            :: epsi0 = 8.85D-12
    real(8),parameter            :: ev    = 1.6D-19
    real(8),parameter            :: emass = 9.1D-31
    real(8),parameter            :: h     = 6.626D-34
    real(8),parameter            :: pi    = 3.1415927
    real(8),parameter            :: c     = 2.98D8
    real(8),parameter            :: g     = 6.61D-11
    integer                      :: i
    integer                      :: j
    integer                      :: k
    real(8)                      :: r
    real(8)                      :: v
    real(8)                      :: x
    character(10)                :: filename 
end module module_common
