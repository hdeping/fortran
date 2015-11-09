module module_common
    implicit none
    integer,parameter            :: n = 100
    integer                      :: i
    integer                      :: j
    integer                      :: k
    real(8)                      :: a(n,n)
    real(8)                      :: b(n,n)
    real(8)                      :: c(n,n)
    integer                      :: ii
    integer                      :: jj
    character(10)                :: filename 

    contains
!subroutine getLeastTime{{{
subroutine getLeastTime()
    integer                  :: i
    integer                  :: j

    do i = 1,10
        print *,i
    enddo !cycle ends
     
end subroutine getLeastTime
!}}}
!subroutine getInitialParameter{{{
subroutine getInitialParameter()
    real(8)              :: x1
    real(8)              :: x2
    integer              :: ii
    integer              :: jj

    do ii = 1,n
        call random_number(x1) 
    enddo !cycle ends
    <++> 
end subroutine getInitialParameter
!}}}
end module module_common
