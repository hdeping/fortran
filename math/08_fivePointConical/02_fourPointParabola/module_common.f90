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
    real(8)                      :: a(m - 1,m)
    real(8)                      :: b(m)
    real(8)                      :: y(m)
    real(8)                      :: x_tmp(m - 1,2)
    real(8)                      :: coef(3) ! coefficient
    real(8)                      :: bvalue(2)
    character(10)                :: filename 
end module module_common
