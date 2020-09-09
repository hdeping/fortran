module module_regre
    interface getbar
     module procedure getbar8, getbar4
    end interface getbar
    contains
!subroutine regre{{{
!  common linear regression
subroutine regre(a,b,r,x,y,n)
    integer,intent(in)        :: n    ! rank of the array
    real(8),intent(in)        :: x(n) ! data input x  
    real(8),intent(in)        :: y(n) ! data input y
    real(8),intent(out)       :: a  ! for slope
    real(8),intent(out)       :: b  ! for intercept
    real(8),intent(out)       :: r  ! for linear factor
    real(8)                   :: x_bar  ! mean value of x 
    real(8)                   :: y_bar  ! mean value of y
    real(8)                   :: xtmp   ! value for I_xx   
    real(8)                   :: ytmp   ! value for I_xy
    real(8)                   :: ztmp   ! value for I_yy
    integer                   :: ii     ! for the cycle
    ! get the mean value of x and y
    x_bar = getbar(x,n) 
    y_bar = getbar(y,n) 
    ! get the value of summations
    ytmp = getmat(x,y,n) 
    xtmp = getmat(x,x,n) 
    ztmp = getmat(y,y,n) 
    ! get the slope intercept and linear factor
    a = ytmp/xtmp
    b = y_bar - a*x_bar
    r = ytmp/sqrt(xtmp*ztmp)
    r = r**2
end subroutine regre
!}}}
!function getmat{{{
function getmat(sa,sb,nn)
    integer,intent(in)        :: nn
    real(8),intent(in)        :: sa(nn)
    real(8),intent(in)        :: sb(nn)
    real(8)                   :: getmat
    real(8)                   :: tmp
    real(8)                   :: valsa  ! mean value  
    real(8)                   :: valsb  ! mean value
    integer                   :: itmp

    ! get mean value
    valsa = getbar(sa,nn) 
    valsb = getbar(sb,nn) 
    !  calculate the function

    tmp = 0.0
    do itmp = 1,nn
        tmp = tmp + (sa(itmp) - valsa)*(sb(itmp) - valsb)
    end do
    getmat = tmp
end function getmat
!}}}
!function getbar8{{{
function getbar8(s,n_len)
    integer,intent(in)           :: n_len
    real(8),intent(in)           :: s(n_len)
    real(8)                      :: getbar8
    real(8)                      :: tmp
    integer                      :: ii

    tmp = 0.0
    do ii = 1,n_len
        tmp = tmp + s(ii)
    end do
    getbar8 = tmp/dble(n_len)
end function getbar8
!}}}
!function getbar4{{{
function getbar4(s,n_len)
    integer,intent(in)           :: n_len
    real   ,intent(in)           :: s(n_len)
    real                         :: getbar4
    real                         :: tmp
    integer                      :: ii

    tmp = 0.0
    do ii = 1,n_len
        tmp = tmp + s(ii)
    end do
    getbar4 = tmp/dble(n_len)
end function getbar4
!}}}
!function getgamma{{{
function getgamma(x)
    real(8),intent(in)       :: x
    real(8)                  ::getgamma
    getgamma = gamma(1 - x)**2.0/gamma(1 - 2.0*x)
end function getgamma
!}}}
!function getgamma2{{{
function getgamma2(x)
    real(8),intent(in)       :: x
    real(8)                  ::getgamma2
    getgamma2 = gamma(1 + x)**2.0/gamma(1 + 2.0*x)
end function getgamma2
!}}}
!function fun{{{
function fun(x,y)
    real(8),intent(in)       :: x
    real(8),intent(in)       :: y
    real(8)                  ::fun
    fun = getgamma2(x) - y
end function fun
!}}}
!function sol_b{{{
function sol_b(x)
    real(8),intent(in)   :: x
    real(8)              :: sol_b
    real(8)              :: y
    real(8)              :: ra
    real(8)              :: rb
    real(8)              :: rc
    integer              :: times 

    ra    = bg
    rb    = ed
    !if(x > 0.499)then
    !    ra = ed
    !endif
    times = 0
    y     = getgamma(x)

    do 
        rc = (ra + rb)/2.0 
        if(fun(rc,y) > 0.0)then
            ra = rc
        else
            rb = rc
        endif
        !  convergence
        if(abs(ra - rb) < delta)exit
    end do
    sol_b = rc


end function sol_b
!}}}
end module module_regre
