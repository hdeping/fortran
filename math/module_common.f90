module module_common
    implicit none
    integer,parameter            :: n = 49
    integer                      :: i
    integer                      :: j
    integer                      :: k
    integer                      :: kk
    integer                      :: jj
    character(10)                :: filename 

    contains
!function tail{{{
function tail(a,b,c,d,e)
    integer,intent(in)        ::  a
    integer,intent(in)        ::  b
    integer,intent(in)        ::  c
    integer,intent(in)        ::  d
    integer,intent(in)        ::  e
    integer                   ::  tail
    integer                   ::  judge(5)
    integer                   ::  ii
    integer                   ::  num1
    integer                   ::  num2

    judge = (/a,b,c,d,e/)

    num1 = 0
    num2 = 0
    do ii = 1,5
        num1 = num1 + judge(ii)
        num2 = num2 + mod(judge(ii),10)
    enddo !cycle ends
    if ( mod(num1,100) == num2 )then
        tail = num2
    else
        tail = 0
    endif ! if ends


    

end function tail
!}}}
end module module_common
