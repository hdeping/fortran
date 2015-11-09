module module_algebra
    use module_common
    implicit none
    contains
! get Ax=b (A for matrix, b, x for vec)
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
!function iter_sol_equ{{{
function iter_sol_equ(a,b,n)
    ! in and out dummies
    integer,intent(in)    :: n
    real(8),intent(in)    :: a(n,n)
    real(8),intent(in)    :: b(n)
    real(8)               :: iter_sol_equ(n)
    real(8)               :: x(n)
    real(8)               :: x_new(n)

    real(8)               :: c(n,n)
    real(8)               :: tmp
    real(8)               :: rate
    integer               :: s
    integer               :: ii
    integer               :: jj
    integer               :: kk

    ! calculate the iteration matrix C
    do i = 1,n
        do j = 1,n
            if(i == j)then
                c(i,j) = 0.0
            else
                c(i,j) = - a(i,j)/a(i,i)
            endif
        end do
    end do
    !do i = 1,n
    !    print "(<n>f12.6)",a(i,:)
    !end do
    !print *,"                 "
    !do i = 1,n
    !    print "(<n>f12.6)",c(i,:)
    !end do
    !pause
    x = b
    do 
        x_new  = vec_multi_add(c,x,b,n)
        lambda = 0.0
        do i = 1,n
            lambda = lambda + abs(x(i) - x_new(i))
        end do
        if(lambda < error)exit
        print *,"lambda = ",lambda
        x = x_new
    end do
    iter_sol_equ = x

end function iter_sol_equ
!}}}
! get AX=B (A, B, X for matrix)
!function sol_mat{{{
function sol_mat(a,b,n)
    ! in and out dummies
    integer,intent(in)    :: n
    real(8),intent(in)    :: a(n,n)
    real(8),intent(in)    :: b(n,n)
    real(8)               :: sol_mat(n,n)

    real(8)               :: c(n,n*2)
    real(8)               :: tmp
    real(8)               :: rate
    integer               :: s
    integer               :: i
    integer               :: j
    integer               :: k
    s = 0
    ! initial the new array c(n,n*2)
    do i = 1,n
        do j = 1,n
            c(i,j)   = a(i,j)
            c(i,n+j) = b(i,j)
        end do
    end do
    !do i = 1,n
    !    print "(<n*2>f8.3)",c(i,:)
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
                do k = i,n*2          !  the ith column to  (n*2)th  column
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
                do k = i,n*2         !  the ith column to  (n*2)th  column
                    c(j,k) = c(j,k)-rate*c(i,k)
                end do   !k
            end do   !j
         end if     
     end do   !i
     do i = 1,n
        do j = 1,n
          sol_mat(i,j) = c(i,n+j)/c(i,i)
        end do
     end do  !i
end function sol_mat
!}}}
!function multi_mat{{{
function multi_mat(a,b,nlen,mlen,llen)
    integer,intent(in)        :: nlen
    integer,intent(in)        :: mlen
    integer,intent(in)        :: llen
    real(8),intent(in)        :: a(nlen,mlen)
    real(8),intent(in)        :: b(mlen,llen)
    real(8)                   :: multi_mat(nlen,llen)
    integer                   :: ii
    integer                   :: jj
    integer                   :: kk

    do ii = 1,nlen
        do jj = 1,llen
            ! get the value of every element
            multi_mat(ii,jj) = 0.0
            do kk = 1,mlen
                multi_mat(ii,jj) =  multi_mat(ii,jj) &
                                    + a(ii,kk)*b(kk,jj)
            end do
        end do
    end do
end function multi_mat
!}}}
!function multi_vec{{{
function multi_vec(a,b,nlen,mlen)
    integer,intent(in)        :: nlen
    integer,intent(in)        :: mlen
    real(8),intent(in)        :: a(nlen,mlen)
    real(8),intent(in)        :: b(mlen)
    real(8)                   :: multi_vec(nlen)
    integer                   :: ii
    integer                   :: kk

    do ii = 1,nlen
        ! get the value of every element
        multi_vec(ii) = 0.0
        do kk = 1,mlen
            multi_vec(ii) =  multi_vec(ii) &
                             + a(ii,kk)*b(kk)
        end do
    end do
end function multi_vec
!}}}
!function mat_add{{{
function mat_add(a,b,nlen,mlen)
    integer,intent(in)        :: nlen
    integer,intent(in)        :: mlen
    real(8),intent(in)        :: a(nlen,mlen)
    real(8),intent(in)        :: b(nlen,mlen)
    real(8)                   :: mat_add(nlen,mlen)
    integer                   :: ii
    integer                   :: kk

    do ii = 1,nlen
        do kk = 1,mlen
            mat_add(ii,kk) =  a(ii,kk) + b(ii,kk) 
        end do
    end do
end function mat_add
!}}}
!function vec_multi_add{{{
! x(n) = a(n,n)*x(n) + b(n) 
function vec_multi_add(a,x,b,nlen)
    integer,intent(in)        :: nlen
    real(8),intent(in)        :: a(nlen,nlen)
    real(8),intent(in)        :: b(nlen)
    real(8),intent(in)        :: x(nlen)
    real(8)                   :: vec_multi_add(nlen)
    real(8)                   :: multi(nlen)
    integer                   :: ii
    integer                   :: kk

    multi = multi_vec(a,x,nlen,nlen)
    do ii = 1,nlen
        vec_multi_add(ii) = multi(ii) + b(ii)
    end do
end function vec_multi_add 
!}}}
end module module_algebra
