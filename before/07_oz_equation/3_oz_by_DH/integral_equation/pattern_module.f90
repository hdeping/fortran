module pattern_module
use parameter_module
use imslf90
use function_module

contains
!subroutine basis_function{{{
subroutine basis_function()
integer, dimension(0:niu)    ::            i_alpha 

i_alpha=(/1, 4, 8, 12, 16, 20, 24/)

p(:, 0)=0.0
p(1, 0)=1.0

do i=1, n-1
    do alpha =1, niu
        if (alpha==1) then
            if (i>=1.and.i<=i_alpha(1)) then
                p(alpha, i)=(dble(i_alpha(1))-dble(i))/dble(i_alpha(1))
            elseif (i>=i_alpha(1).and.i<=n) then
                p(alpha, i)=0.0
            endif
        else
            if (i>=1.and.i<=i_alpha(alpha-2)) then
                p(alpha, i)=0.0
            elseif (i>=i_alpha(alpha-2).and.i<=i_alpha(alpha-1)) then
                p(alpha, i)=(dble(i)-dble(i_alpha(alpha-2)))/(dble(i_alpha(alpha-1))-dble(i_alpha(alpha-2)))
            elseif (i>=i_alpha(alpha-1).and.i<=i_alpha(alpha)) then
                p(alpha, i)=(dble(i_alpha(alpha))-dble(i))/(dble(i_alpha(alpha))-dble(i_alpha(alpha-1)))
            elseif (i>=i_alpha(alpha).and.i<=n-1) then
                p(alpha, i)=0.0
            endif
        endif
    enddo
enddo

end subroutine basis_function
!}}}
!subroutine conjugates_function{{{
subroutine conjugates_function()
integer                            ::            i
integer                            ::            j
integer                            ::            k
real*8, dimension(niu, niu)        ::            rab
real*8, dimension(niu, niu)        ::            bab    

do i=1, niu 
    do j=1, niu
        rab(i, j)=0.0
        do k=ib, n-1
            rab(i, j)=rab(i, j)+p(i, k)*p(j, k)
        enddo
    enddo
enddo

bab=.i.rab

do i=1, niu
    do j=ib, n-1
        q(i, j)=0.0
        do k=1, niu
            q(i, j)=q(i, j)+bab(i, k)*p(k, j)
        enddo
    enddo
enddo

end subroutine conjugates_function
!}}}
!subroutine jacobian{{{
subroutine jacobian()
integer                ::        i, j, ii, jj
real*8                 ::        tq

do i=1, niu
    do j=1, niu
        tq=0.0
        do ii=ib, n-1
            do jj=ib, n-1
                if (q(i, ii)==0.0.or.p(j, jj)==0.0) then
                    cycle
                else
                    tq=tq+q(i, ii)*gij(ii, jj)*p(j, jj)
                endif
            enddo
        enddo
        jac(i, j)=delta(i, j)-tq
    enddo
enddo

end subroutine jacobian
!}}}

end module pattern_module
