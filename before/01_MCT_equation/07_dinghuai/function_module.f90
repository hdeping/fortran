module function_module
use parameter_module


contains

!--****************************************************************--!

real*8 function memory(i, l)
integer            ::            i, l, m


memory=0.0
!$omp parallel private(m) if (i>500)
!$omp do reduction(+:memory) 
    do m=1, i
        memory=memory+ke_me(l, m)*dphi(l, i-m)
    enddo
!$omp end do
!$omp end parallel
memory=(gammak(l)/ak(l))*memory
return

end function memory

!--******************************************************************--!

real*8 function kernel(l)
integer                ::            l, k, p

kernel=0.0
if (l==0) then
    kernel=0.0d0
else
!$omp parallel private(k, p)
!$omp do reduction(+:kernel)
    do k=1, n-1
        do p=max(1, abs(l-k)), l+k         
            if (p<=n-1) then
                kernel=kernel+((dble(k)*dble(p))/(dble(l)**5.))*(((dble(l)**2.+dble(k)**2.-dble(p)**2.)*ck(k)+(dble(l)**2.+dble(p)**2.-dble(k)**2.)*ck(p))**2.) &
                        *sk(l)*sk(k)*sk(p)*phi(k)*phi(p)
            elseif (p>n-1) then
                kernel=kernel+((dble(k)*dble(p))/(dble(l)**5.))*(((dble(l)**2.+dble(k)**2.-dble(p)**2.)*ck(k))**2.) &
                        *sk(l)*sk(k)*1.0*phi(k)*phi(n-1)
            endif
        enddo
    enddo
!$omp end do
!$omp end parallel
endif
kernel=((rho*(h**3.))/(32.0*(pi**2.)))*kernel

end function kernel

!--*******************************************************************--!

end module function_module
