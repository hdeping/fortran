module module_common
    implicit none
!variables{{{
    integer,parameter         :: n = 100
    integer,parameter         :: m = 4
    integer                   :: ierror
    integer                   :: i
    integer                   :: j
    real(8)                   :: x(n,m)
    real(8)                   :: y(n)
    real(8)                   :: a(m)
    real(8)                   :: a_new(m)
    real(8)                   :: b
    real(8)                   :: b_new
    real(8)                   :: r
    real(8)                   :: ax
    real(8)                   :: bx
    real(8)                   :: rx
    real(8)                   :: x1
    real(8)                   :: sx(n)
    real(8)                   :: sy(n)
    real(8)                   :: slope 
    real(8)                   :: intercept
    real(8)                   :: mat_trya(3,3)
    real(8)                   :: object(3)
    real(8)                   :: answer(3)
    character(20)             :: filename
!}}}
end module module_common
