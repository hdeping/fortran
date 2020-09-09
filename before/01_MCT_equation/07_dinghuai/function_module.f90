module function_module
use parameter_module
    
contains
    
!--********************************************--!
    
real*8 function aphief(l, i)
integer                ::                i, l
    
    aphief=au(l)*phie(l, i)
    
return    
end function aphief

!--********************************************--!

real*8 function memorye(l, i)
integer                ::                i, l, m        
 
!$omp parallel private(m) if (i>5000)
!$omp do reduction(+:memorye)
memorye=0.0d0
do m=1, i
    memorye=memorye+ke_mee(l, m)*dphie(l, i-m)
enddo
!$omp end do
!$omp end parallel

return
end function memorye

!--********************************************--!

real*8 function kernele(l, i)
integer                ::                i, l, k, p

if (l==0) then
    kernele=0.0d0
else
    kernele=0.0d0
    do k=1, n-1
        do p=max(1, abs(l-k)), l+k
            if (p<=n-1) then
                kernele=kernele+(dble(k*p))*(((dble(k)**2.+dble(l)**2.-dble(p)**2.)*ck(k)+ &
                       (dble(p)**2.+dble(l)**2.-dble(k)**2.)*ck(p))**2.)*sk(k)*sk(p)*phie(k, i)*phie(p, i)
            else
                kernele=kernele+(dble(k*p))*(((dble(k)**2.+dble(l)**2.-dble(p)**2.)*ck(k))**2.)&
                    *sk(k)*1.0*phie(k, i)*phie(n-1, i)
            endif
        enddo
    enddo
    kernele=aus(l)*sk(l)*kernele
endif

end function kernele

!--********************************************--!

end module function_module
