module module_common
    implicit none
    !integer,parameter            :: n = 42
    integer,parameter            :: m = 2
    integer                      :: i
    integer                      :: j
    integer                      :: k
    integer,allocatable          :: a(:)  ! denotes status 
    integer,allocatable          :: b(:)  
    integer                      :: c(m-1)
    integer                      :: times ! denotes the cycle number
    integer                      :: countNum
    character(10)                :: filename 
    contains
!function getThrees{{{
function getThrees(a,n)
    integer,intent(in)       :: n
    integer,intent(in)       :: a(n)
    integer                 :: getThrees
    integer                 :: ii

    getThrees = 0
    do ii = 1,n
        if ( a(ii) == 3 )then
            getThrees = getThrees + 1 
        endif ! if ends
    enddo !cycle ends
end function getThrees
!}}}
!function getFinal{{{
function getFinal(n,m)
    integer,intent(in)     :: n
    integer,intent(in)     :: m
    integer                :: getFinal(m - 1)

    a = 0
    times = 0
    countNum = 0
    do 
        times = times + 1
        if ( times > n )then
            times = mod(times,n)
        endif ! if ends
        
        if ( a(times) == m )then
            cycle
        endif ! if ends

        k = k + 1
        if ( k > m )then
            k = mod(k,m) 
        endif ! if ends
        a(times) = k
        if ( k == m )then
            countNum = countNum + 1
            b(countNum) = times
        endif ! if ends
        if (countNum == n)then
            exit
        endif ! if ends
    enddo !cycle ends

    do k = n+2-m,n
        getFinal(k-n+m-1) = b(k)
    enddo !cycle ends

end function getFinal

!}}}

end module module_common
