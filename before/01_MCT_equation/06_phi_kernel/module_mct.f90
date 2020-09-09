module module_mct
    use module_common

contains
!subroutine getmatrixes{{{
subroutine getmatrixes()
    !real(8)             :: tmpvalue(m,m)
    integer              :: ii
    integer              :: jj
    real(8)              :: tovl ! total value
    real(8)              :: tmpa

    !  get r and k
    do i = 1,n
        dk(i) = h*dble(i)
    end do
    filename = "sk.txt"
    open(20,file = filename, status = "old", iostat = ierror)
    ! read Sk and Ck
    do i = 1,ncut
        read(20,*,iostat = ierror)tmpa,sk(i)
        !print *,sk(i)
        !pause
        if(ierror /= 0)exit
    end do
    close(20)
    filename = "ck.txt"
    open(20,file = filename, status = "old",iostat = ierror)
    do i = 1,ncut
        read(20,*,iostat = ierror)tmpa,ck(i)
        if(ierror /= 0)exit
    end do
    close(20)


    diffu = 1.0
    !  get K and R
    do q = 1,ncut
        mat_K(q)  = dk(q)**2.0*diffu*v
        mat_R(q)  = mat_K(q)/sk(q)
        !print *,mat_K(q),mat_R(q)
    end do
    !  get matrix A
    do q = 1,ncut
        tovl = h**7.0/(16.0*pi**2.0*dble(q)*rho)
       ! print *,"tovl = ",tovl
       ! pause
        do k = 1,ncut
            do p = 1,ncut
                l1 = dble(q)**2.0 + dble(k)**2.0 - dble(p)**2.0
                l2 = dble(q)**2.0 + dble(p)**2.0 - dble(k)**2.0
                mat_A(p,q,k) = tovl*p*k*l1*ck(k)*(l1*ck(k) + &
                               l2*ck(p))
                !print *,mat_A(p,q,k)
                !pause
            end do
        end do
    end do

end subroutine getmatrixes
!}}}
!subroutine getmct{{{
!get MCT equation
subroutine getmct()
   integer             :: ii 
   integer             :: jj
   real(8)             :: tmpa
   real(8)             :: rec1  ! for the time cost
   real(8)             :: rec2  ! for the time cost


   ! get initial f(q,1)
   do ii = 1,ncut
       f(ii,1) = sk(ii)
   end do
   ! get matrix U
   ! get memory(q,1)
   memory(:,1) = getmemory(1)
   !print *,memory(10:20,1)
   !pause
   do q = 1,ncut
       mat_U(q) = 1.0 + dt*mat_K(q)*memory(q,1) 
       !print *,mat_U(q)
       !pause
   end do
   !print *,mat_U(10:20)
   !print *,mat_K(10:20)
   !pause
   do t = 1,tmnum - 1
       !call cpu_time(rec1)
       if(t /= 1)then
           memory(:,t) = getmemory(t)
       endif
       !call cpu_time(rec2)
       !print *,"time cost is ",rec2 - rec1
       do q  = 1,ncut
           tmpa  = - dt*mat_R(q)*f(q,t)
           !print *,tmpa,dt
           !pause
           do ii = 1,t - 1
               tmpa = tmpa - dt*mat_K(q)*memory(q,&
                      t + 1 - ii)*par_f(q,ii)
           end do
           par_f(q,t) = tmpa/mat_U(q)
           f(q,t + 1) = f(q,t) + par_f(q,t)
       end do
       !print *,par_f(90:100,t)
       !pause
   end do
end subroutine getmct
!}}}
!subroutine getmct2{{{
!get MCT equation
subroutine getmct2()
   integer             :: ii 
   integer             :: jj
   real(8)             :: tmpa
   real(8)             :: rec1  ! for the time cost
   real(8)             :: rec2  ! for the time cost


   do jj = 1,ncut
    do ii = 1,tmnum/2
        f(jj,ii)    = f(jj,ii*2)
     memory(jj,ii)  = memory(jj,ii*2)
    end do
   end do
   
   do t = tmnum/2,tmnum - 1
       !call cpu_time(rec1)
           memory(:,t) = getmemory(t)
       !call cpu_time(rec2)
       !print *,"time cost is ",rec2 - rec1
       do q  = 1,ncut
           tmpa  = - dt*mat_R(q)*f(q,t)
           !print *,tmpa,dt
           !pause
           do ii = 1,t - 1
               tmpa = tmpa - dt*mat_K(q)*memory(q,&
                      t + 1 - ii)*par_f(q,ii)
           end do
           par_f(q,t) = tmpa/mat_U(q)
           f(q,t + 1) = f(q,t) + par_f(q,t)
       end do
       !print *,par_f(90:100,t)
       !pause
   end do
end subroutine getmct2
!}}}
!get memory(q,t){{{
function getmemory(t)
   integer,intent(in)  :: t
   real(8)             :: getmemory(ncut)

   do q = 1,ncut
       tmp = 0.0
       do k = 1,ncut !  cut-off
           !do p = abs(q - k),q + k
           do p = abs(q - k),ncut !   cut-off
               if(p == 0)cycle
               tmp  = tmp + mat_A(p,q,k)*f(k,t)*f(p,t)
           end do
       end do
       getmemory(q) = tmp
   end do
end function getmemory
!}}}
!get newmemory(q,t){{{
function getnew_memory(input)
   real(8)             :: input(ncut)
   real(8)             :: getnew_memory(ncut)

   do q = 1,ncut
       tmp = 0.0
       do k = 1,ncut !  cut-off
           !do p = abs(q - k),q + k
           do p = abs(q - k),ncut !   cut-off
               if(p == 0)cycle
               tmp  = tmp + mat_A(p,q,k)*input(k)*input(p)
           end do
       end do
       getnew_memory(q) = tmp
   end do
end function  getnew_memory
!}}}
!function getfinal{{{
function getfinal()
    real(8)              :: getfinal(ncut)
    real(8)              :: newmemory(ncut)
    real(8)              :: tmpa(ncut)
    real(8)              :: tmpb(ncut)
    integer              :: ii
    integer              :: times 

    tmpa = 1.0 
    do 
        newmemory = getnew_memory(tmpa)
        do ii = 1,ncut
            !tmpb(ii) = newmemory(ii)*(sk(ii) - tmpa(ii))/sk(ii)
            tmpb(ii) = newmemory(ii)*sk(ii)/(sk(ii) + newmemory(ii))
            !print *,"memory = ",newmemory(ii)
            !print *,"tmpb = ",tmpb(ii),"tmpa = ",tmpa(ii)
            !pause
        end do
        lambda = judge(tmpa,tmpb,ncut)
        if(lambda < error)exit
        if(mod(times,fre) == 0)then
            print *,"lambda = ",lambda
        endif
        tmpa = tmpb
    end do


    getfinal = tmpa
     

end function getfinal
!}}}
!function judge{{{
! compare the difference between two arrays
! judge the convergence
function judge(a,b,n)
    integer,intent(in)      :: n
    real(8),intent(in)      :: a(n)
    real(8),intent(in)      :: b(n)
    real(8)                 :: judge
    integer                 :: ii
    
    judge = 0
    do ii = 1,n 
        judge = judge + abs(a(ii) - b(ii))
    end do
    !judge = judge/dble(n)
    
end function judge
!}}}

end module module_mct
