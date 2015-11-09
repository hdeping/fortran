module function_module
use parameter_module

contains

!function delta{{{
real*8 function delta(i, j)
integer            ::            i, j

if (i==j) then
    delta=1.0
else
    delta=0.0
endif
return
end function delta
!}}}
!function fst{{{
function fst(indata, tt)
real*8, dimension(ib:n-1)            ::        indata
real*8, dimension(:), pointer        ::        fst
integer                              ::        tt

allocate(fst(ib:n-1))

fst=0.0
if (tt==1) then
    do j=1, n-1
        do i=1, n/2-1
            fst(j)=fst(j)+r(i)*indata(i)*sin(k(j)*r(i))
        enddo
        fst(j)=fst(j)*((4.0*pi*dr)/k(j))
    enddo
elseif (tt==-1) then
    do i=1, n-1
        do j=1, n/2-1
            fst(i)=fst(i)+k(j)*indata(j)*sin(k(j)*r(i))
        enddo
        fst(i)=fst(i)*(dk/(2.0*(pi**2.)*r(i)))
    enddo
endif

end function fst
!}}}
!function f_gamma{{{
function f_gamma(indata_fi, indata_se)
real*8, dimension(niu)                ::        indata_fi
real*8, dimension(ib:n-1)             ::        indata_se
real*8, dimension(:), pointer         ::        f_gamma

allocate(f_gamma(ib:n-1))
f_gamma(:)=0.0
do i=ib, n-1
    do j=1, niu
        f_gamma(i)=f_gamma(i)+indata_fi(j)*p(j, i)
    enddo
    f_gamma(i)=f_gamma(i)+indata_se(i)
enddo

end function f_gamma
!}}}
!function d_gamma{{{
function d_gamma(indata_fi, indata_se)
real*8, dimension(niu)               ::        indata_se
real*8, dimension(ib:n-1)            ::        indata_fi
real*8, dimension(:), pointer        ::        d_gamma
real*8                               ::        ta

allocate(d_gamma(ib:n-1))
d_gamma(:)=0.0
do i=ib, n-1
    ta=0.0
    do j=1, niu
        ta=ta+indata_se(j)*p(j, i)
    enddo
    d_gamma(i)=indata_fi(i)-ta
enddo

end function d_gamma
!}}}
!function a_alpha{{{
function a_alpha(indata)
real*8, dimension(ib:n-1)            ::        indata
real*8, dimension(:), pointer        ::        a_alpha

allocate(a_alpha(niu))
a_alpha(:)=0.0

do i=1, niu
    do j=ib, n-1
        a_alpha(i)=a_alpha(i)+q(i, j)*indata(j)
    enddo
enddo

end function a_alpha
!}}}
!function mayer{{{
function mayer(indata)
real*8, dimension(ib:n-1)            ::        indata
real*8, dimension(:), pointer        ::        mayer

allocate(mayer(ib:n-1))

do i=ib, n-1
    if (i==0) then
        mayer(i)=0.0
    else
        if (r(i)<di) then
            mayer(i)=-1.0
        else
            mayer(i)=0.0
        endif
    endif
!    mayer(i)=exp(-beta*vf(r(i)))-1.0
enddo

end function mayer
!}}}
!function gij{{{
real*8 function gij(i, j)
integer            ::            i, j

if (i==0) then
    if (abs(fmayer(j))<1.0e-5) then
        gij=0.0
    else
        gij=((2.0*dr*r(j))/pi)*(fmayer(j))*e(j)
    endif
else
    if (abs(fmayer(j))<1.0e-5) then
        gij=0.0
    else
        gij=((dr*r(j))/(pi*r(i)))*(fmayer(j))*(d(r(i)-r(j))-d(r(i)+r(j)))
    endif
endif
return
end function gij
!}}}
!function vf{{{
real*8 function vf(dis)
real*8            ::            dis

vf=4.0*((1.0/dis)**6.0)
return
end function vf
!}}}
!function d{{{
real*8 function d(dis)
real*8             ::            dis
integer            ::            m

    d=0.0
do m=0, n/2-1
    d=d+(((2.0*rho*ck(m))/(1.0-rho*ck(m)))+((rho*ck(m))/(1.0-rho*ck(m)))**2.)*cos(k(m)*dis)
enddo
    d=dk*d
return
end function d
!}}}
!function e{{{
real*8 function e(i)
integer              ::            i, m

e=0.0
do m=0, n/2-1
    e=e+k(m)*(((2.0*rho*ck(m))/(1.0-rho*ck(m)))+((rho*ck(m))/(1.0-rho*ck(m)))**2.)*sin(k(m)*r(i))
enddo
e=dk*e
return
end function e
!}}}



end module function_module
