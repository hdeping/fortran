module module_mct
    use module_common

contains
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
            tmpb(ii) = newmemory(ii)*sk(ii)/(1 + newmemory(ii))  !/sk(ii)
            !print *,tmpb(ii),tmpa(ii)
            !pause
        end do
        lambda = judge(tmpa,tmpb,ncut)
        if(lambda < error)exit
        if(mod(times,fre) == 0)then
            print "('lambda = ',D12.5)",lambda
            !print *,"fre = " ,fre
            !pause
        endif
        tmpa = tmpb
    end do
    print "('lambda = ',D12.5)",lambda


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
