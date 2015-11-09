module module_common
!variables{{{
    implicit none
    integer,parameter            :: n     = 400
    integer,parameter            :: m     = 40 ! rank of polynomial
    real(8),parameter            :: pi    = 3.141592653
    real(8),parameter            :: dt    = 0.01
    real(8),parameter            :: error = 1E-6
    integer                      :: i
    integer                      :: j
    integer                      :: k
    character(10)                :: filename 
    real(8)                      :: x
    real(8)                      :: y
    real(8)                      :: r
    real(8)                      :: theta
    real(8)                      :: tmp
    real(8)                      :: a(m+1) ! coefficient
    real(8)                      :: t1
    real(8)                      :: t2
    real(8)                      :: solx(2)
    real(8)                      :: roots(m,2)
    integer                      :: icount
    integer                      :: cyTimes
    integer                      :: times
    real(8)                      :: judge
    integer                      :: ierror
!}}}
    contains
!function getall{{{
function getall()
    real(8)                     :: getall(m,2)
    real(8)                     :: c(m,2)
    integer                     :: ii

    call cpu_time(t1)
    times   = 0
    cyTimes = 0
    do   ii = 1,1000
        call random_number(r)
        call random_number(theta)
        r     = 100*r
        theta = 2.0*pi*theta
        solx  = getroot(r,theta)
        !print *,solx(:)
        !pause
        !  judge different
        if(times == 0)then
            times      = times + 1
            c(times,:) = solx(:)
            if(abs(solx(2)) > 100.0*error)then 
                times = times + 1 
                c(times,1) =   solx(1) 
                c(times,2) = - solx(2) 
            endif
        else
            icount = 0
            do j = 1,times - 1
                judge = geterror(c(j,:),solx,2)
                if(judge < error)icount = icount + 1
            end do
            if(icount == 0)then
                 times      = times + 1
                 c(times,:) = solx(:)
                 if(abs(solx(2)) > 100.0*error)then 
                     times = times + 1 
                     c(times,1) =   solx(1) 
                     c(times,2) = - solx(2) 
                 endif
            endif
        endif
        if(times == m)exit
        cyTimes = cyTimes + 1
    end do
    call cpu_time(t2)

    getall = c
end function getall
!}}}
!function getroot{{{
function getroot(r,theta)
    real(8),intent(in)          :: r
    real(8),intent(in)          :: theta
    real(8)                     :: getroot(2)
    !integer                     :: getroot
    real(8)                     :: sol(2)
    real(8)                     :: mat_a(2,2)
    real(8)                     :: mat_b(2)
    real(8)                     :: delsol(2)
    real(8)                     :: test(2)
    real(8)                     :: judge
    integer                     :: ii
    integer                     :: jj
    integer                     :: itmp
    integer                     :: times

    times = 0
    sol   = (/r,theta/)
    !print *,sol
    do 
        times = times + 1
        mat_a(1,1) = getg_r(sol(1),sol(2),a)
        mat_a(1,2) = getg_theta(sol(1),sol(2),a)
        mat_a(2,1) = geth_r(sol(1),sol(2),a)
        mat_a(2,2) = geth_theta(sol(1),sol(2),a)
        mat_b(1)   = getg(sol(1),sol(2),a)
        mat_b(2)   = geth(sol(1),sol(2),a)
        ! get delsol
        delsol     = sol_equ(mat_a,mat_b,2)
        judge = 0.0
        do ii = 1,2
            judge  = judge + abs(delsol(ii))
        end do
        !print *,"judge = ",judge
        !pause
        !do ii = 1,2
        !    print *,mat_a(ii,:),mat_b(ii)
        !end do
        sol    = sol - delsol
        !print *,sol(:)
        if(judge < error)exit
    end do
    getroot(1) = sol(1)*cos(sol(2))
    getroot(2) = sol(1)*sin(sol(2))
    !pause
    !getroot = times


end function getroot
!}}}
!subroutine getr_theta{{{
subroutine getr_theta(r,theta,x,y)
    real(8),intent(in)          :: x
    real(8),intent(in)          :: y
    real(8),intent(out)         :: r
    real(8),intent(out)         :: theta
    integer                     :: ii

    r = sqrt(x**2.0 + y**2.0)
    !  first judge about theta
    if(x == 0.0)then
        if(y > 0.0)then
            theta = 0.0
        else
            theta = pi
        endif
    else
        theta  = atan(y/x) 
    endif
    
    !  second judge about theta
    if(x < 0.0)then
        if(y > 0.0)then
            theta = theta + pi/2.0
        else
            theta = theta - pi/2.0
        endif
    endif
end subroutine getr_theta
!}}}
!function getNorm{{{
function getNorm(x,y,a)
    real(8),intent(in)          :: a(m+1)
    real(8),intent(in)          :: x
    real(8),intent(in)          :: y
    real(8)                     :: getNorm
    real(8)                     :: tmpx
    real(8)                     :: tmpy
    integer                     :: ii
    call getr_theta(r,theta,x,y)
    !tmpx    = getg(x,y,a)
    !tmpy    = geth(x,y,a)
    tmpx    = getg(r,theta,a)
    tmpy    = geth(r,theta,a)
    getNorm = sqrt(tmpx**2.0 + tmpy**2.0)

end function getNorm
!}}}
!function getPoly{{{
function getPoly(x,a,dima)
    real(8),intent(in)          :: x
    integer,intent(in)          :: dima
    real(8),intent(in)          :: a(dima)
    real(8)                     :: getPoly
    integer                     :: ii

    getPoly = 0.0
    do ii = dima,1,- 1
        getPoly = getPoly*x + a(ii)
    end do
end function getPoly
!}}}
!function getg{{{
function getg(r,theta,a)
    real(8),intent(in)          :: r
    real(8),intent(in)          :: theta
    real(8),intent(in)          :: a(m+1)
    real(8)                     :: b(m+1)
    real(8)                     :: getg
    integer                     :: ii

    ! get the coefficients of g(r,theta)
    do ii = 1,m+1
         b(ii) = a(ii)*cos((ii - 1)*theta) 
         !print *,theta,b(ii)
         !pause
    end do
    getg = getPoly(r,b,m+1)

end function getg
!}}}
!function geth{{{
function geth(r,theta,a)
    real(8),intent(in)          :: r
    real(8),intent(in)          :: theta
    real(8),intent(in)          :: a(m+1)
    real(8)                     :: b(m+1)
    real(8)                     :: geth
    integer                     :: ii

    ! get the coefficients of g(r,theta)
    do ii = 1,m+1
         b(ii) = a(ii)*sin((ii - 1)*theta) 
    end do
    geth = getPoly(r,b,m+1)

end function geth
!}}}
!function getg_r{{{
function getg_r(r,theta,a)
    real(8),intent(in)          :: r
    real(8),intent(in)          :: theta
    real(8),intent(in)          :: a(m+1)
    real(8)                     :: b(m)
    real(8)                     :: getg_r
    integer                     :: ii

    ! get the coefficients of g(r,theta)
    do ii = 1,m
         b(ii) = ii*a(ii+1)*cos(ii*theta) 
    end do
    getg_r = getPoly(r,b,m)

end function getg_r
!}}}
!function getg_theta{{{
function getg_theta(r,theta,a)
    real(8),intent(in)          :: r
    real(8),intent(in)          :: theta
    real(8),intent(in)          :: a(m+1)
    real(8)                     :: b(m+1)
    real(8)                     :: getg_theta
    integer                     :: ii

    ! get the coefficients of g(r,theta)
    do ii = 1,m+1
         b(ii) = (ii - 1)*a(ii)*sin((ii - 1)*theta) 
    end do
    getg_theta = - getPoly(r,b,m+1)

end function getg_theta
!}}}
!function geth_r{{{
function geth_r(r,theta,a)
    real(8),intent(in)          :: r
    real(8),intent(in)          :: theta
    real(8),intent(in)          :: a(m+1)
    real(8)                     :: b(m)
    real(8)                     :: geth_r
    integer                     :: ii

    ! get the coefficients of g(r,theta)
    do ii = 1,m
         b(ii) = ii*a(ii+1)*sin(ii*theta) 
    end do
    geth_r = getPoly(r,b,m)

end function geth_r
!}}}
!function geth_theta{{{
function geth_theta(r,theta,a)
    real(8),intent(in)          :: r
    real(8),intent(in)          :: theta
    real(8),intent(in)          :: a(m+1)
    real(8)                     :: b(m+1)
    real(8)                     :: geth_theta
    integer                     :: ii

    ! get the coefficients of g(r,theta)
    do ii = 1,m+1
         b(ii) = (ii - 1)*a(ii)*cos((ii - 1)*theta) 
    end do
    geth_theta = getPoly(r,b,m+1)

end function geth_theta
!}}}
!function sol_equ{{{
function sol_equ(a,b,n)
    ! in and out dummies
    integer,intent(in)    :: n
    real(8),intent(in)    :: a(n,n)
    real(8),intent(in)    :: b(n)
    real(8)               :: sol_equ(n)

    real(8)               :: c(n,n+1)
    real(8)               :: tmp
    real(8)               :: rate
    integer               :: s
    integer               :: i
    integer               :: j
    integer               :: k
    s = 0
    ! initial the new array c(n,n+1)
    do i = 1,n
        do j = 1,n
            c(i,j) = a(i,j)
        end do
        c(i,n+1) = b(i)
    end do
    !do i = 1,n
    !    print "(<n+1>f8.3)",c(i,:)
    !end do
    
     do  i = 1,n
         if(c(i,i) == 0)then
            do j = i+1,n
                if(c(j,i) /= 0)exit
            end do     !  i
            if(j == n+1)then
                s = 1
                cycle
            else
                do k = i,n+1          !  the ith column to  (n+1)th  column
                      tmp = c(i,k)
                   c(i,k) = c(j,k)
                   c(j,k) = tmp
               end do  !k
            endif
         endif
         if(s == 0)then
            do j = 1,n
                if(j == i)cycle
                if(c(j,i) == 0)cycle
                rate = c(j,i)/c(i,i)
                do k = i,n+1         !  the ith column to  (n+1)th  column
                    c(j,k) = c(j,k)-rate*c(i,k)
                end do   !k
            end do   !j
         end if     
     end do   !i
     do i = 1,n
         sol_equ(i) = c(i,n+1)/c(i,i)
     end do  !i
end function sol_equ
!}}}
!function multi_vec{{{
function multi_vec(mata,matb,s)
    integer,intent(in)          :: s
    real(8),intent(in)          :: mata(s,s)
    real(8),intent(in)          :: matb(s)
    real(8)                     :: multi_vec(s)
    integer                     :: ii
    integer                     :: jj

    do ii = 1,s
       multi_vec(ii) = 0.0
       do jj = 1,s
           multi_vec(ii) = multi_vec(ii) + mata(ii,jj)*matb(jj)
       end do
    end do

end function multi_vec
!}}}
!function gettheta{{{
function gettheta(theta)
    real(8),intent(in)          :: theta
    real(8)                     :: gettheta
    integer                     :: ii
    integer                     :: jj

    gettheta = theta
    if (gettheta < 0.0)then
        do while (gettheta < 0.0)
            gettheta = gettheta + 2.0*pi
        enddo  
    elseif(theta > 2.0*pi)then
        do while (gettheta > 2.0*pi)
            gettheta = gettheta - 2.0*pi
        enddo  
    endif


end function gettheta
!}}}
!function geterror{{{
function geterror(array_a,array_b,s)
    integer,intent(in)          :: s
    real(8),intent(in)          :: array_a(s)
    real(8),intent(in)          :: array_b(s)
    real(8)                     :: geterror
    integer                     :: ii

    geterror = 0.0
    do ii = 1,s
        geterror = geterror + abs(array_a(ii) - array_b(ii))
    end do


end function geterror
!}}}
!subroutine check{{{
subroutine check()
    integer                     :: ii
    real(8)                     :: valuetmp

    !filename = "data.txt"
    !open(10,file = filename,status = "old",iostat = ierror)
    !ii = 0
    !do 
    !    ii = ii + 1
    !    if(ii <= m + 1)then
    !        read(10,*,iostat = ierror)a(ii)
    !    else
    !        read(10,*,iostat = ierror)x,y
    !        call getr_theta(r,theta,x,y)
    !        print *,"norm is => ",getNorm(r,theta,a)
    !    endif 
    !    if(ierror /= 0)exit
    !end do

    !close(10)

    do ii = 1,m
        valuetmp = getNorm(roots(ii,1),roots(ii,2),a)
        if(abs(valuetmp) < error)print *,roots(ii,:)
    end do



end subroutine check
!}}}

end module module_common
