module module_regre
    use module_algebra
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
!subroutine multi_regre{{{
! multi-dimensional linear regression
! see the formulas in the tex document
subroutine multi_regre(a,b,r,x,y,n,m)
    integer,intent(in)        :: n
    integer,intent(in)        :: m
    real(8),intent(in)        :: x(n,m)
    real(8),intent(in)        :: y(n)
    real(8),intent(out)       :: a(m)
    real(8),intent(out)       :: b
    real(8),intent(out)       :: r
    real(8)                   :: mat(m,m)
    real(8)                   :: maty(m)
    real(8)                   :: y_bar
    real(8)                   :: x_bar(m)
    real(8)                   :: yytmp
    real(8)                   :: tmp
    integer                   :: ii
    integer                   :: jj

    !  get the value of matrix
    do ii = 1,m
       do jj = ii,m
         mat(ii,jj) = getmat(x(:,ii),x(:,jj),n) 
       end do
         maty(ii) = getmat(x(:,ii),y,n) 
    end do
    do ii = 2,m
        do jj = 1,ii - 1
          mat(ii,jj) = mat(jj,ii)
        end do
    end do
    !  get the value of array a
    a = sol_equ(mat,maty,m)
    ! get the mean value of x,y
    do ii = 1,m
        x_bar(ii) = getbar(x(:,ii),n)
    end do
    y_bar = getbar(y,n)
    !  get the value of b
    b = y_bar
    do ii = 1,m
        b = b - a(ii)*x_bar(ii)
    end do
    !  get the value of r
    yytmp = getmat(y,y,n)
    tmp = 0.0
    do ii = 1,m
        tmp = tmp + a(ii)*maty(ii)
    end do
    r = sqrt(tmp/yytmp)
!*********** FINISHED *******************************
end subroutine multi_regre
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
end module module_regre
