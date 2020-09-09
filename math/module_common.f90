module module_common
    implicit none
    integer,parameter            :: n = 9
    integer,parameter            :: total = 3**(n - 1) - 1
    integer                      :: i
    integer                      :: j
    integer                      :: k
    integer                      :: a(n)
    integer                      :: b(n)
    character(len = 10)          :: filename
    
    
    contains
! get ternary notation
!function getTernary{{{
function getTernary(ii)
    integer ,intent(in)    :: ii
    integer                :: getTernary(n - 1)
    integer                :: tmp
    integer                :: jj

    tmp = ii

    do jj = n - 1,1,- 1
        getTernary(jj) = mod(tmp,3)
        tmp = tmp / 3
    enddo !cycle ends
end function getTernary
!}}}
end module module_common
