module module_laplace
    implicit none
    integer,parameter            :: n   = 10
    integer,parameter            :: num = 1000
    real(8),parameter            :: dt  = 1E-2
    real(8),parameter            :: ln2 = log(2.0)
    real(8)                      :: facn
    integer                      :: i
    integer                      :: j
    integer                      :: k
    integer                      :: a(2*n)
    real(8)                      :: inver_f(num)
    real(8)                      :: t1
    real(8)                      :: t2
    real(8)                      :: tvalue
    character(10)                :: filename 

    contains
!subroutine get_inver{{{
subroutine get_inver()
    integer                      :: ii
    integer                      :: jj
    integer                      :: kk
    real(8)                      :: tmp
    real(8)                      :: x
    
    filename = "data.txt"
    open(10,file = filename)

    facn  = factorial(n)
    call getmatA()
    x   = dt
    tmp = lap_inver(x)
    write(10,"(2f18.6)")x,1.0
    do ii = 2,num
        x           = dt*ii
        inver_f(ii) = lap_inver(x)
        inver_f(ii) = inver_f(ii)/tmp
        write(10,"(2f18.6)")x,inver_f(ii)
    end do
    close(10)

end subroutine get_inver
!}}}
!function lap_inver{{{
! get n1!
function lap_inver(x)
    real(8),intent(in)               :: x
    real(8)                          :: lap_inver
    real(8)                          :: tmp
    real(8)                          :: tmp1
    integer                          :: ii

    lap_inver = 0.0
    tmp       = ln2/x

    !print *,"laplace = ",lap_inver
    do ii = 1,n*2
        tmp1      = fun(dble(ii)*tmp)
        lap_inver = lap_inver + a(ii)*tmp1
    end do
    lap_inver = tmp*lap_inver
    !print *,"x       = ",x
    !print *,"laplace = ",lap_inver
    !pause
end function lap_inver
!}}}
!subroutine getmatA{{{
subroutine getmatA()
    integer                      :: ii
    integer                      :: jj
    real(8)                      :: tmp

    do ii = 1,2*n
        tmp = 0.0
        do jj = (ii + 1)/2,min(ii,n)
        tmp   = tmp + jj**(n+1)*combi(n-jj,jj)*combi(jj,jj)&
                *combi(ii-jj,2*j-ii)
        end do
        a(ii) = tmp*symbol(n+ii)/facn
        !print *,"ii = ",ii
        !print *,"tmp   = ",tmp
        !print *,"a(ii) = ",a(ii)
        !print *,"facn  = ",facn
        !pause
    end do

end subroutine getmatA
!}}}
!function fun{{{
function fun(x)
    real(8),intent(in)               :: x
    real(8)                          :: fun
    fun = 1.0/(x**3.0)
end function fun
!}}}
!function symbol{{{
function symbol(n)
    integer,intent(in)               :: n
    real(8)                          :: symbol
    if(mod(n,2) == 0)then
        symbol = 1.0
    else
        symbol = - 1.0
    endif
end function symbol
!}}}
!function combi{{{
! get the combination 
!number of n1 and n2
! (n1+n2)!/n1!/n2!
function combi(n1,n2)
    integer,intent(in)               :: n1
    integer,intent(in)               :: n2
    integer                          :: ii
    real(8)                          :: combi

    combi = 1.0
    do ii = 1,n1
        combi = combi*dble((n2 + ii))/dble(ii)
    end do
end function combi
!}}}
!function factorial{{{
! get n1!
function factorial(n1)
    integer,intent(in)               :: n1
    integer                          :: ii
    real(8)                          :: factorial

    factorial = 1.0
    do ii = 1,n1
        factorial = factorial*dble(ii)
    end do
end function factorial
!}}}
end module module_laplace
