module module_common
    implicit none
    integer,parameter            :: n = 1000
    integer,parameter            :: m = n/10
    integer,parameter            :: total = int(1E6)
    integer                      :: i
    integer                      :: j
    integer                      :: k
    character(10)                :: filename 
    contains
!function getSigma{{{
function getSigma(a,n)
    integer ,intent(in)     :: n
    integer ,intent(in)     :: a(n)
    real                    :: getSigma
    real                    :: x
    real                    :: tmp
    integer                 :: jj

    x = sum(a)/ dble(n)
    tmp = 0.0
    do jj = 1,n
        tmp = tmp + (a(jj) - x)**2.0 
    enddo !cycle ends
    getSigma  = sqrt(tmp)
end function getSigma
!}}}
end module module_common
