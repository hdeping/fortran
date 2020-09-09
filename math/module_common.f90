module module_common
    implicit none
    integer,parameter            :: n = 100
    integer,parameter            :: d = 2
    integer,parameter            :: m = d*(d + 3)/2 
    integer                      :: i
    integer                      :: j
    integer                      :: k
    integer                      :: ierror
    integer                      :: num
    integer                      :: total
    real(8)                      :: x(m,d)
    real(8)                      :: a(m,m)
    real(8)                      :: b(m)
    real(8)                      :: y(m)
    character(10)                :: filename 
end module module_common
