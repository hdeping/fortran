module module_common
    implicit none
    integer,parameter            :: n     = 3
    real(8),parameter            :: error = 1D-8
    integer                      :: i
    integer                      :: j
    integer                      :: k
    character(10)                :: filename 
    contains
!!function mat_sqare_root{{{
!!  if X^2=A, how to get A (matrices)
!function mat_sqare_root(a,n)
!    integer,intent(in)         :: n
!    real(8),intent(in)         :: a(n,n)
!    real(8)                    :: mat_sqare_root(n,n)
!    real(8)                    :: tmp(n,n)
!    real(8)                    :: tmpa(n,n)
!    real(8)                    :: tmpb(n,n)
!    real(8)                    :: lambda
!    integer                    :: ii
!    integer                    :: jj
!
!    !  get the random value of the matrix
!    do ii = 1,n
!        do jj = 1,n
!            mat_sqare_root(ii,jj) = sqrt(abs(a(ii,jj))/dble(n))
!        enddo !cycle ends
!    enddo !cycle ends
!    do ii = 1,n
!        print *,a(ii,:)
!    enddo !cycle ends
!     
!
!    !  iteration begins
!    do 
!       tmpa   = sol_mat(mat_sqare_root,a,n)
!       tmpb   = tran(sol_mat(tran(mat_sqare_root,n),tran(a,n),n),n)
!       tmp    = (tmpa + tmpb)/2.0
!       lambda = diff(mat_sqare_root,n)
!       if ( lambda < error )then
!           exit
!       endif ! if ends
!       ! get a new one
!       do ii = 1,n
!            do jj = 1,n
!                mat_sqare_root(ii,jj) = (mat_sqare_root(ii,jj) + &
!                     tmp(ii,jj))/2.0
!            enddo !cycle ends
!       enddo !cycle ends
!       !mat_sqare_root = tmp
!       print *,"lambda = ",lambda
!       pause
!    enddo !cycle ends
!     
!
!end function mat_sqare_root
!!}}}
!function mat_sqare_root{{{
!  if X^2=A, how to get A (matrices)
function mat_sqare_root(a,n)
    integer,intent(in)         :: n
    real(8),intent(in)         :: a(n,n)
    real(8)                    :: mat_sqare_root(n,n)
    real(8)                    :: tmp(n,n)
    real(8)                    :: lambda
    integer                    :: ii
    integer                    :: jj

    !  get the random value of the matrix
    do ii = 1,n
        do jj = 1,n
            mat_sqare_root(ii,jj) = sqrt(abs(a(ii,jj))/dble(n))
        enddo !cycle ends
    enddo !cycle ends
    do ii = 1,n
        print *,a(ii,:)
    enddo !cycle ends
     

    !  iteration begins
    do 
       tmpa   = sol_mat(mat_sqare_root,)
       lambda = diff(mat_sqare_root,n)
       if ( lambda < error )then
           exit
       endif ! if ends
       ! get a new one
       do ii = 1,n
            do jj = 1,n
                mat_sqare_root(ii,jj) = mat_sqare_root(ii,jj) - &
                                   tmp(ii,jj)
            enddo !cycle ends
       enddo !cycle ends
       print *,"lambda = ",lambda
       pause
    enddo !cycle ends
     

end function mat_sqare_root
!}}}
!!function mat_sqare_root{{{
!!  if X^2=A, how to get A (matrices)
!function mat_sqare_root(a,n)
!    integer,intent(in)         :: n
!    real(8),intent(in)         :: a(n,n)
!    real(8)                    :: mat_sqare_root(n,n)
!    real(8)                    :: tmp(n,n)
!    real(8)                    :: tmpa(n,n)
!    real(8)                    :: tmpb(n,n)
!    real(8)                    :: lambda
!    integer                    :: ii
!    integer                    :: jj
!
!    !  get the random value of the matrix
!    do ii = 1,n
!        do jj = 1,n
!                mat_sqare_root(ii,jj) = sqrt(abs(a(ii,jj))/dble(n))
!        enddo !cycle ends
!    enddo !cycle ends
!    do ii = 1,n
!        print *,a(ii,:)
!    enddo !cycle ends
!     
!
!    !  iteration begins
!    do 
!       tmpa   = minus(multi(mat_sqare_root,&
!                     mat_sqare_root,n),a,n)
!       tmp    = 0.5*sol_mat(mat_sqare_root,tmpa,n)
!       lambda = diff(mat_sqare_root,n)
!       if ( lambda < error )then
!           exit
!       endif ! if ends
!       ! get a new one
!       do ii = 1,n
!            do jj = 1,n
!                mat_sqare_root(ii,jj) = mat_sqare_root(ii,jj) - &
!                                   tmp(ii,jj)
!            enddo !cycle ends
!       enddo !cycle ends
!       print *,"lambda = ",lambda
!       pause
!    enddo !cycle ends
!     
!
!end function mat_sqare_root
!!}}}
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
!function diff{{{
function diff(a,n)
    ! in and out dummies
    integer,intent(in)    :: n
    real(8),intent(in)    :: a(n,n)
    real(8)               :: diff
    integer               :: ii
    integer               :: jj
    
    diff = 0.0
    do ii = 1,n
        do jj = 1,n
            !diff = diff + abs(a(ii,jj) - b(ii,jj))
            diff = diff + abs(a(ii,jj))
        enddo !cycle ends
    enddo !cycle ends
     
end function diff
!}}}
!function multi{{{
function multi(a,b,n)
    ! in and out dummies
    integer,intent(in)    :: n
    real(8),intent(in)    :: a(n,n)
    real(8),intent(in)    :: b(n,n)
    real(8)               :: multi(n,n)
    integer               :: ii
    integer               :: jj
    integer               :: kk
    
    multi = 0.0
    do ii = 1,n
        do jj = 1,n
            do kk = 1,n
                multi(ii,jj) = multi(ii,jj) + a(ii,kk)*b(kk,jj)
            enddo !cycle ends
        enddo !cycle ends
    enddo !cycle ends
     
end function multi
!}}}
!function minus{{{
function minus(a,b,n)
    ! in and out dummies
    integer,intent(in)    :: n
    real(8),intent(in)    :: a(n,n)
    real(8),intent(in)    :: b(n,n)
    real(8)               :: minus(n,n)
    integer               :: ii
    integer               :: jj
    integer               :: kk
    
    do ii = 1,n
        do jj = 1,n
            minus(ii,jj) = a(ii,jj) - b(ii,jj)
        enddo !cycle ends
    enddo !cycle ends
     
end function minus
!}}}
!function tran{{{
function tran(a,n)
    ! in and out dummies
    integer,intent(in)    :: n
    real(8),intent(in)    :: a(n,n)
    real(8)               :: tran(n,n)
    integer               :: ii
    integer               :: jj
    integer               :: kk
    
    do ii = 1,n
        do jj = 1,n
            tran(ii,jj) = a(jj,ii)
        enddo !cycle ends
    enddo !cycle ends
     
end function tran
!}}}
end module module_common
