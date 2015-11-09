module evolution_module
use parameter_module
use pattern_module    
        
contains
    ! subroutine initial_parameter{{{
subroutine initial_parameter()
integer            ::            i
real*8             ::            ra
    
    call basis_function()
    call conjugates_function()
    
    dr=0.05; dk=(2.0*pi)/(n*dr)
    r(0)=0.0; k(0)=0.0
    do i=1, n-1
        r(i)=dr*dble(i); k(i)=dk*dble(i)
    enddo
    
    fmayer(:)=mayer(r(:))

    ksi=rho*(4./3.)*pi*((di/2.)**3.)
    gamma_1=((1.0+2.0*ksi)**2.)/((1.0-ksi)**4.)
    gamma_2=-((2.0+ksi)**2.)/(4.0*(1.0-ksi)**4.)

    do i=0, n-1
        if (r(i)<1.0) then
            cr(i)=-gamma_1-6.0*ksi*gamma_2*r(i)-0.5*ksi*gamma_1*(r(i)**3.)
        else
            cr(i)=0.0
        endif
    enddo
    
    open(1, file='cr.txt')
        do i=0, n-1
            if (i==0) then
                write(1, '(i5, f18.8)') i, 0.0
            else
                write(1, '(i5, f18.8)') i, cr(i)
            endif
        enddo
    close(1)
    pause
    delta_2=5.0
    call initial_gamma()

    aa(:)=a_alpha(gammar(:))
    dgamma(:)=d_gamma(gammar(:), aa(:))    


end subroutine initial_parameter

    !}}}
!subroutine evolution{{{
subroutine evolution()

    
!--***************part_1******************--!
    gammar(:)=f_gamma(aa(:), dgamma(:))
!--***************************************--!    

!--***************part_2******************--!
    cr(:)=(1.0+gammar(:))*fmayer(:)
    ck(:)=fst(cr(:), 1)

    gammak(:)=(rho*(ck(:)**2.))/(1.0-rho*ck(:))
    pgammar(:)=fst(gammak(:), -1)
    
    call jacobian()
!--***************************************--!    

!--***************part_3******************--!
    
    paa(:)=a_alpha(pgammar(:))
    pdgamma(:)=d_gamma(pgammar(:), paa(:))


    do i=1, niu
        da(i)=aa(i)-paa(i)
    enddo

    call error_da(*100, *400)
!--***************************************--!
    
100    step=step+1
    call error_gamma(*200, *300)
    
200 open(1, file='g.txt')
    do i=1, n-1    
        write(1, '(i5, f18.8)') i, gammar(i)+cr(i)+1.0
    enddo
    close(1)

    open(1, file='cr.txt')
    do i=1, n-1
        write(1, '(i5, f18.8)') i, cr(i)
    enddo
    close(1)
    stop

300 dgamma(:)=pdgamma(:)

    return

400 call newton_raphson()

end subroutine evolution
!}}}
!subroutine error_gamma{{{
subroutine error_gamma(*, *)
integer            ::        i
real*8            ::        t, eta

t=0.0
do i=1, n-1
    t=t+(pgammar(i)-gammar(i))**2.    
enddo

eta=sqrt(dr*t)
print*, 'a', eta
if (eta<1.0e-5) then
    return 1
else
    return 2 
endif

end subroutine error_gamma
!}}}
!subroutine newton_raphson{{{
subroutine newton_raphson()
integer                   ::        i, j
real*8, dimension(niu)    ::        t
real*8, dimension(niu)    ::        naa
real*8                    ::        sde

ijac=.i.jac
t=0.0
do i=1, niu
    do j=1, niu
        t(i)=t(i)+ijac(i, j)*(da(j))
    enddo
enddo

sde=0.0
do i=1, niu
    sde=sde+t(i)**2.
enddo

if (sde<delta_2) then
    naa(:)=aa(:)-t(:)
    aa(:)=naa(:)
else
    do i=1, niu
        t(i)=(sqrt(delta_2)/sqrt(sde))*t(i)
    enddo
    naa(:)=aa(:)-t(:)
    aa(:)=naa(:)
endif

end subroutine newton_raphson
!}}}
!subroutine error_da{{{
subroutine error_da(*, *)
integer            ::        i
real*8             ::        sd

sd=0.0
do i=1, niu
    sd=sd+da(i)**2.
enddo
print*, sd
if (sqrt(sd)<1.0e-5) then
    return 1
else
    return 2
endif

end subroutine error_da
!}}}
!subroutine initial_gamma{{{
subroutine initial_gamma()

do i=0, n-1
    cr(i)=cx(r(i))
enddo

ck(:)=fst(cr(:), 1)

gammak(:)=rho*(ck(:)**2.)/(1.0-rho*ck(:))
gammar(:)=fst(gammak(:), -1)

end subroutine initial_gamma
!}}}
!function cx{{{
real*8 function cx(dis)
real*8            ::        dis

if (dis<1.0) then
    cx=-gamma_1-6.0*ksi*gamma_2*(dis/di)-0.5*ksi*gamma_1*(dis/di)**3.
else
    cx=0.0
endif
return
end function cx
!}}}
end module evolution_module
