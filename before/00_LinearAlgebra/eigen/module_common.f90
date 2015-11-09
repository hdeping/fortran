module module_algebra
    implicit none
    integer                 ::  m
    integer                 ::  ierror                    !  for the control of the file
    integer                 ::  i,j,k                     ! for the cycle  in the subroutines
    integer                 ::  itmp,jtmp                 !  for the cycle  in the main  
    !real                    ::  a(n,n)                    !  array 
    !real                    ::  b(n,n)                    !  array 
    !real                    ::  c(n,n)                    !  array 
    !real                    ::  re(n)                     !  array 
    !real                    ::  x(n),test(n)
    integer                 ::   ii,jj
    real                    ::  t1,t2
    real*8                  ::  x1,x2
    
    !integer   order(n)
    real  alpha
    real  beta
    real  gamma
    real  s
    character*(20)  filename
    
    contains
!subroutine multi_mat{{{
subroutine multi_mat(x,a,b,n)
    integer  n
    real  a(n,n),b(n,n)
    real  x(n,n)
    x=0
    do i=1,n
        do j=1,n
            do k=1,n
                x(i,j)=x(i,j)+a(i,k)*b(k,j)
                
            end do  !k
        end do  !j
    end do  !i
    end subroutine multi_mat
!}}}
!subroutine multi_mat_vec{{{
subroutine multi_mat_vec(b,a,x,n)
    integer  n
    real x(n),a(n,n)
    real b(n)
    do i=1,n
        b(i)=0
        do j=1,n
            b(i)=b(i)+a(i,j)*x(j)
        end do  !j
    end do  !i
    end subroutine multi_mat_vec
!}}}
!subroutine inverse{{{
subroutine inverse(b,a,n)
    integer   n
    real  a(n,n),b(n,n)
    real  c(n,2*n),tmp,rate
     integer  s,order(n),value
     c=0
     do i=1,n
         order(i)=i
     end do  !i
     s=0
    do i=1,n
        do j=1,n
            c(i,j)=a(i,j)
        end do  !j
        c(i,i+n)=1
    end do   !i
    !do i=1,n
    !print '(<2*n>f10.3)',c(i,:)
    !end do  !i
    !read*
     do  i=1,n
         if(c(i,i)==0)then
         do j=i+1,n
             if(c(j,i)/=0)exit
         end do     !  i
         if(j==n+1)then
             s=1
             cycle
         else
             do k=i,n*2          !  the ith column to  (2n)th  column
                 tmp=c(i,k)
                 c(i,k)=c(j,k)
                 c(j,k)=tmp
                 value=order(i)
                 order(i)=order(j)
                 order(j)=value
            end do  !k
         endif
         end if
         if(s==0)then
         do j=1,n
             if(j==i)cycle
             if(c(j,i)==0)cycle
             rate=c(j,i)/c(i,i)
             do k=i,n*2         !  the ith column to  (2n)th  column
                 c(j,k)=c(j,k)-rate*c(i,k)
             end do   !k
         end do   !j
         end if     
     end do   !i
     do i=1,n
         do j=n+1,2*n
             c(i,j)=c(i,j)/c(i,i)
             b(i,j-n)=c(i,j)
         end do  !j
         c(i,i)=1
     end do   !i
    ! write(*,*)"the final form is "
    ! do i=1,n
    !print '(<2*n>f10.3)',c(i,:)
    !end do  !i
    !read*
    end subroutine inverse
!}}}
!subroutine tr{{{
    !  calculate the trace of the matrix
subroutine tr(s,a,n)
    integer n
    real  s,a(n,n),tmp
    
    tmp=0
    do ii=1,n
        tmp=tmp+a(ii,ii)
    end do    ! ii
    s=tmp
    end subroutine tr
!}}}
!subroutine  eigen_vec{{{
    !  calculate the eigen vector
subroutine  eigen_vec(b,a,n)
    integer  n
    real a(n,n),b(n,n),c(n),tmp(n,n)
    call eigen(c,a,n)
    do ii=1,n
        tmp=a
        do jj=1,n
            tmp(jj,jj)=tmp(jj,jj)-c(ii)
        end do   !   jj
        call diagonal2(tmp,tmp,n)
        c(n)=1
        do jj=1,n-1
            c(jj)=-tmp(jj,n)/tmp(jj,jj)
        end do   !   jj
        call normal(c,c,n)
        b(:,ii)=c
    end do  !  jj
    end subroutine eigen_vec
!}}}
!subroutine ortho{{{
subroutine ortho(b,a,n,m)
    integer  n,m
    real  b(n,m),a(n,m)
    real  c(n),stmp
    call normal(c,a(:,1),n)
    a(:,1)=c
    do ii=2,m
        do jj=1,ii-1
            stmp=dot_product(a(:,ii),a(:,jj))
            a(:,ii)=a(:,ii)-stmp*a(:,jj)
        end do   ! jj
        call normal(c,a(:,ii),n)
        a(:,ii)=c
    end do   ! ii
    end subroutine ortho
!}}}
!subroutine  normal{{{
subroutine  normal(b,a,n)
    integer  n
    real  b(n),a(n)
    real  tmp
    tmp=0
    do ii=1,n
        tmp=tmp+a(ii)**2
    end do   ! ii
    tmp=tmp**0.5
    do ii=1,n
        b(ii)=a(ii)/tmp
    end do  !  ii
    end subroutine 
!}}}
!subroutine eigen{{{
    ! calculate the eigen value of the symmetric matrix
subroutine eigen(b,a,n)
    integer  n
    real  a(n,n),b(n),c(n,n)
    call diagonal2(c,a,n)
    do i=1,n
        b(i)=c(i,i)
    end do  !i
    end subroutine eigen
!}}}
!subroutine det{{{
    ! calculate the determination of the matrix
subroutine det(s,a,n)
    integer  n
    real  a(n,n),b(n,n)
    real  s,tmp
    call diagonal(b,a,n)
    tmp=1
    do i=1,n
        tmp=tmp*b(i,i)
    end do  !i
    s=tmp
    end subroutine det
!}}}
!subroutine diagonal2{{{
    ! diagonalize the matrix  (full)
subroutine diagonal2(b,a,n)
     integer  n
     real  b(n,n),a(n,n)
     real  tmp
     call diagonal(b,a,n)
     do i=1,n
         do j=1,i-1
             tmp=b(i,j)
             b(i,j)=b(j,i)
             b(j,i)=tmp
         end do
     end do  !i
     call diagonal(b,b,n)
    end subroutine diagonal2
!}}}
!subroutine diagonal{{{
    ! diagonalize the matrix  (half)
subroutine diagonal(b,a,n)
     integer  n
     real  b(n,n),a(n,n)
     real  tmp,rate
     integer  s
     s=0
     b=a
     do  i=1,n-1
         if(b(i,i)==0)then
         do j=i+1,n
             if(b(j,i)/=0)exit
         end do     !  i
         if(j==n+1)then
             s=1
             cycle
         else
             do k=i,n
                 tmp=b(i,k)
                 b(i,k)=b(j,k)
                 b(j,k)=tmp
            end do  !k
         endif
         end if
         if(s==0)then
         do j=i+1,n
             if(b(j,i)==0)cycle
             rate=b(j,i)/b(i,i)
             do k=i,n
                 b(j,k)=b(j,k)-rate*b(i,k)
             end do   !k
         end do   !j
         end if     
     end do   !i
     
     
     end subroutine diagonal
!}}}
    end module module_algebra
