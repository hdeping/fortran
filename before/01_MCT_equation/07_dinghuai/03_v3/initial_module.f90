module initial_module
use parameter_module
use fst_module


contains


!--******************************************--!

subroutine initial_parameter()
integer                    ::            i
character(len=10)        ::            filename

v0=0.0
di=1.0
packphi=0.512
rho=packphi/(pi/6.0)
mf=1.0
df=(v0**2.)/3.0

allocate(phie(0:n-1, 0:sn))
allocate(dphie(0:n-1, 0:sn))
allocate(ke_mee(0:n-1, 0:sn))

call calculate_oz(ck, sk)

open(1, file='sk.txt')
    do i=0, n-1
        write(1, '(2f25.16)') dble(i)*dk, sk(i)
    enddo
close(1)

do i=1, n-1    
    au(i)=(1.0+sk(i)*(v0**2.)/3.0)*(((dble(i)*dk)**2.)/sk(i))
    bu(i)=(((dble(i)*dk)**2.)/sk(i))/(1.0+sk(i)*(v0**2.)/3.0)
    aus(i)=(rho*((h)**3.))/(32.0*(pi**2.)*((dble(i))**5.))
enddo

au(0)=0.0d0; bu(0)=0.0d0
aus(0)=0.0d0

det_time=1.0e-6    

do i=1, n-1
    write(filename, '(i3.3, a4)') i, '.txt'
    open(unit=i, file='data/'//trim(filename))
enddo

end subroutine initial_parameter

!--******************************************--!

subroutine calculate_oz(tck, tsk)
integer                                ::                i
integer, parameter                    ::                on=1024, l=10
real*8                                ::                eta
real*8                                ::                lambta_1, lambta_2
real*8, dimension(0:on-1)            ::                r, k
real*8, dimension(0:on-1)            ::                cr, ccr, cck, chk, hk, sk, ck
real*8, dimension(0:n-1)            ::                tsk, tck

dr=0.01d0; dk=pi/(dble(on)*dr)
h=dk

do i=0, on-1
    r(i)=dble(i)*dr; k(i)=dble(i)*dk
enddo

eta=rho*pi/6.0

lambta_1=((1.0+2.0*eta)**2.)/((1.0-eta)**4.)
lambta_2=-((1.0+0.5*eta)**2.)/((1.0-eta)**4.)

do i=1, on-1
    if (r(i)<1.0d0) then
        cr(i)=-lambta_1-6.0*eta*lambta_2*r(i)-0.5*eta*lambta_1*((r(i))**3.)
    else
        cr(i)=0.0d0
    endif
enddo
cr(0)=0.0d0


ccr(:)=cr(:)*r(:)
cck(:)=fst(ccr(:), on, l, 1)
hk(:)=cck(:)/(k(:)-rho*cck(:))
hk(0)=0.0d0
sk(:)=1.0+rho*hk(:)
sk(0)=0.0d0
ck(:)=cck(:)/k(:)
ck(0)=0.0d0

tck(:)=ck(0:n-1); tsk(:)=sk(0:n-1)

end subroutine calculate_oz

!--*******************************************--!

end module initial_module
