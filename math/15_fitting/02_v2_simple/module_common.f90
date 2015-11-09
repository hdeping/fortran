module module_common
    implicit none
    integer,parameter            :: n     = 20
    integer,parameter            :: m     = 6
    integer,parameter            :: fre   = 100
    real(8),parameter            :: dt    = 1D-2
    real(8),parameter            :: error = 1D-8
    integer                      :: i
    integer                      :: j
    integer                      :: k
    real(8)                      :: x(n)
    real(8)                      :: y(n)
    real(8)                      :: a_para(m)
    real(8)                      :: b_para(m)
    real(8)                      :: anew(m)
    real(8)                      :: bnew(m)
    real(8)                      :: lambda
    character(10)                :: filename 

    contains
!subroutine getpara{{{
!  get the parameter with 
!  Newton-Raphason method
subroutine getpara(a,b,x,y)
    ! in and out dummies
    real(8),intent(out)   :: a(m)
    real(8),intent(out)   :: b(m)
    real(8),intent(in)    :: x(n)
    real(8),intent(in)    :: y(n)
    real(8)               :: para(2*m)  !  parameter
    real(8)               :: med(2*m)   !  medium variable
    real(8)               :: jacob(2*m,2*m)  !  jacobi matrix
    real(8)               :: equ(2*m)
    integer               :: ii
    integer               :: jj
    integer               :: times


    ! initial parameter
    call random_number(para(:))
    times = 0
    do 
        equ   = equations(para(1:m),para((m+1):2*m),x,y)
        jacob = getpartial(para(1:m),para((m+1):2*m),x,y)
        med   = sol_equ(jacob,equ,2*m)
        !do ii = 1,2*m
        !    print *,jacob(ii,:)
        !end do
        !print *,equ(:)
        !print *,med(:)
        !pause
        !  judge convergence
        lambda = judge(med,2*m)
        if(lambda < error)exit
        times = times + 1
        if(mod(times,fre) == 0)then
            print *,"lambda = ",lambda
            !print *,med(:)
            !pause
        endif
        !  get new para
        para = para - med
    end do
    print *,"lambda = ",lambda

    a = para(1:m)
    b = para((1+m):2*m)
    
end subroutine getpara
!}}}
!function getpartial{{{
function getpartial(a,b,x,y)
    ! in and out dummies
    real(8),intent(in)    :: a(m)
    real(8),intent(in)    :: b(m)
    real(8),intent(in)    :: x(n)
    real(8),intent(in)    :: y(n)
    real(8)               :: ynew(n)
    real(8)               :: equ(2*m)
    real(8)               :: var_x
    real(8)               :: getpartial(2*m,2*m)
    real(8)               :: tmp = 0.0
    integer               :: ii
    integer               :: jj
    integer               :: i1
    integer               :: j1

    
    !  get temporary varibles
    ynew  = funx(x,a,b)
    do ii = 1,n
        ynew(ii) =  y(ii) - ynew(ii)
    end do
    !  get equations
    equ = equations(a,b,x,y)
    !  part 1
    do i1 = 1,m
        do j1 = 1,m
           tmp = 0.0
           do ii = 1,n
                tmp = tmp - exp(- (b(i1) + b(j1))*x(ii))
           end do
           getpartial(i1,j1) = tmp 
        end do
    end do
    !  part 2
    do i1 = 1,m
        do j1 = 1,m
           tmp = 0.0
           do ii = 1,n
               tmp = tmp + a(j1)*x(ii)*exp(- (b(i1) + b(j1))*x(ii))
           end do
           if(i1 == j1)then
               tmp = tmp - equ(m + i1)
           endif
           getpartial(i1,j1 + m) = tmp 
        end do
    end do
    !  part 3
    do i1 = 1,m
        do j1 = 1,m
           tmp = 0.0
           do ii = 1,n
                tmp = tmp - x(ii)*exp(- (b(i1) + b(j1))*x(ii))
           end do
           getpartial(i1 + m,j1) = tmp 
        end do
    end do
    !  part 4
    do i1 = 1,m
        do j1 = 1,m
           tmp = 0.0
           do ii = 1,n
               tmp = tmp + a(j1)*x(ii)**2.0*exp(- (b(i1) + b(j1))*x(ii))
           end do
           if(i1 == j1)then
               do jj = 1,n
                   tmp = tmp - x(jj)**2.0*exp(- b(ii)*x(jj))*ynew(jj)
               end do
           endif
           getpartial(i1 + m,j1 + m) = tmp 
        end do
    end do

end function getpartial
!}}}
!function equations{{{
function equations(a,b,x,y)
    ! in and out dummies
    real(8),intent(in)    :: a(m)
    real(8),intent(in)    :: b(m)
    real(8),intent(in)    :: x(n)
    real(8),intent(in)    :: y(n)
    real(8)               :: ynew(n)
    real(8)               :: equations(2*m)
    real(8)               :: tmp 
    integer               :: ii
    integer               :: jj
    integer               :: i1
    integer               :: j1

    ynew  = funx(x,a,b)
    do ii = 1,n
        ynew(ii) =  y(ii) - ynew(ii)
        print *,ynew(ii)
        pause
    end do
    !  get 1~m
    do ii = 1,m
        tmp = 0.0
        do jj = 1,n
            tmp = tmp + exp(- b(ii)*x(jj))*ynew(jj)
        end do
        equations(ii) = tmp
    end do
    !  get m+1~2*m
    do ii = 1,m
        tmp = 0.0
        do jj = 1,n
            tmp = tmp + x(jj)*exp(- b(ii)*x(jj))*ynew(jj)
        end do
        equations(ii + m) = tmp
    end do

end function equations
!}}}
!function funx{{{
function funx(x,a,b)
    ! in and out dummies
    real(8),intent(in)    :: x(n)
    real(8),intent(in)    :: a(m)
    real(8),intent(in)    :: b(m)
    real(8)               :: var_x
    real(8)               :: funx(n)
    integer               :: ii
    integer               :: jj

    do jj = 1,n
        funx(jj)  = 0.0
        do ii = 1,m
            funx(jj)  = funx(jj) + a(ii)*exp(- b(ii)*x(jj))
        end do
    end do
end function funx
!}}}
!function judge{{{
function judge(a,n)
    integer,intent(in)        :: n
    real(8),intent(in)        :: a(n)
    real(8)                   :: judge
    integer                   :: ii

    judge = 0.0
    do ii = 1,n
       judge = judge + abs(a(ii))
    end do

end function judge
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
end module module_common
