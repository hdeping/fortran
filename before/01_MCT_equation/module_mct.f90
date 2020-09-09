module module_mct
    use module_common

contains
!function getfinal{{{
function getfinal()
    real(8)              :: getfinal(ncut)
    real(8)              :: newmemory(ncut)
    real(8)              :: tmpa(ncut)
    real(8)              :: tmpb(ncut)
    real(8)              :: tmpc(ncut)
    integer              :: ii
    integer              :: times 

    ! get deltafun
    deltafun = 0.0
    do ii = 1,ncut
        deltafun(ii,ii) = 1.0
    end do
    tmpa     = 1.0 
    do 
        newmemory = getnew_memory(tmpa)
        do ii = 1,ncut
            tmpc(ii) = tmpa(ii) + newmemory(ii)*&
                       (tmpa(ii) - sk(ii))
            !pause
        end do
        tmpb = sol_equ(jacobi(tmpa,newmemory),tmpc,ncut)
        lambda = newjudge(tmpb,ncut)
        if(lambda < error)exit
        if(mod(times,fre) == 0)then
            print "('lambda = ',D12.5)",lambda
            !print *,"fre = " ,fre
            !pause
        endif
        tmpa = tmpa - tmpb
    end do
    print "('lambda = ',D12.5)",lambda


    getfinal = tmpa
     

end function getfinal
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
               tmp  = tmp + input(k)*input(p)
           end do
       end do
       getnew_memory(q) = tmp/dble(ncut*smul)
   end do
end function  getnew_memory
!}}}
!get partial memory(q,t){{{
function jacobi(input,newmemory)
   real(8),intent(in)  :: input(ncut)
   real(8),intent(in)  :: newmemory(ncut)
   real(8)             :: jacobi(ncut,ncut)
   real(8)             :: tcost1
   real(8)             :: tcost2

   do l = 1,ncut
      call cpu_time(tcost1)
      do q = 1,ncut
          tmp = 0.0
          do k = 1,ncut !  cut-off
              !do p = abs(q - k),q + k
              do p = abs(q - k),ncut !   cut-off
                  if(p == 0)cycle
                  tmp  = tmp + input(k)*deltafun(p,l)&
                         + input(p)*deltafun(k,l)
              end do
          end do
          tmp  = tmp*(input(q) - sk(q))/dble(ncut*smul)
          if(q == l)then
              tmp = tmp + memory(q) + 1.0
          endif
          jacobi(q,l) = tmp
      end do
      call cpu_time(tcost2)
      print *,"time cost is ",tcost2 - tcost1
      pause
   end do
end function jacobi 

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
!function newjudge{{{
! compare the difference between two arrays
! newjudge the convergence
function newjudge(a,n)
    integer,intent(in)      :: n
    real(8),intent(in)      :: a(n)
    real(8)                 :: newjudge
    integer                 :: ii
    
    newjudge = 0.0
    do ii = 1,n 
        newjudge = newjudge + abs(a(ii))
    end do
end function newjudge
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

end module module_mct
