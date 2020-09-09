module module_common
    integer,parameter        ::   n = 20
    integer,parameter        ::   m = 30
    integer                  ::   a(n)
    integer                  ::   b(m)
    integer                  ::   c(m + n)
    integer                  ::   i
    integer                  ::   j
    integer                  ::   stmp
    ! 9*9 multiplication table
    integer                  ::   table(9,9) 
    ! check the big number algorithm
    integer(8)               ::   anum
    integer(8)               ::   bnum
    integer(8)               ::   cnum
    real                     ::   t1
    real                     ::   t2
    real                     ::   x1
    character(10)            ::   filename

      
    contains
!function multiply{{{
!****** calculate the multiplication    **********
!****** of two big numbers  with digits **********
!******             n and m             **********
function multiply(a,b,n,m)
    integer,intent(in)         :: n
    integer,intent(in)         :: m
    integer,intent(in)         :: a(n)
    integer,intent(in)         :: b(m)
    integer                    :: multiply(n+m)
    integer                    :: ii
    integer                    :: jj

    multiply = 0
    do ii = n,1,-1 
        do jj = m,1,-1
           multiply(ii + jj) = multiply(ii + jj) + &
                               table(a(ii),b(jj)) 
        end do
    end do
    call dot_carry_one(multiply,n+m)

end function multiply
!}}}
!function getMultiplyTable{{{
! get 9*9 multiplication table
function getMultiplyTable()
    integer              ::   getMultiplyTable(9,9)
    integer              ::   ii
    integer              ::   jj
    integer              ::   tmp
    integer              ::   tmp1
    integer              ::   tmp2
    do ii = 1,9
        do jj = 1,9
            tmp  = ii*jj
            getMultiplyTable(ii,jj) = tmp
        end do
    end do
    
end function getMultiplyTable
!}}}
!subroutine dot_carry_one{{{
! get dot and carry one 
subroutine dot_carry_one(a,n)
    integer,intent(in)        :: n
    integer,intent(inout)     :: a(n)
    integer                   :: tmp

    do i = n,2,-1
       tmp      = a(i)
       a(i)     = mod(tmp,10)
       a(i - 1) = a(i - 1) + tmp/10
    end do
end subroutine dot_carry_one
!}}}

end module module_common
