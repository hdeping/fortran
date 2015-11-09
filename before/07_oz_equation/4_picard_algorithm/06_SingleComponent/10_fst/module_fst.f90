module module_fst
     use module_common
contains

!function fst{{{
function fst(id,ju)
integer                       :: ju
integer                       :: i, j
real*8, dimension(0:n-1)      :: id
real*8, dimension(0:n-1)      :: tbr, tbi, tor, toi
real*8, dimension(0:n-1)      :: fst
real*8                        :: th
!%%%%%%%%%%%%%%%%%%NOTICE%%%%%%%%%%%%%%%%%%%%!
!    dr should be equal to deltar
!%%%%%%%%%%%%%%%%%%NOTICE%%%%%%%%%%%%%%%%%%%%!

tbr(0)=0.0
tbi=0.0
do i=1, n-1
    tbr(i)=sin(pi*dble(i)/dble(n))*(id(i)+id(n-i))+0.5*(id(i)-id(n-i))
enddo

call fft(l,n,tbr, tbi, tor, toi)

fst(0)=0.0
fst(1)=0.5*tor(0)
do j=1, n/2-1
    i=j+j
    fst(i)=toi(j)
    fst(i+1)=fst(i-1)+tor(j)
enddo

if (ju==1) then
    th=4.0*deltar*pi
    fst(:)=th*fst(:)
elseif (ju==-1) then
    th=deltak/(2.0*(pi**2.))
    fst(:)=th*fst(:)
endif

end function fst
!}}}
!!original function fst{{{
!function fst(id, n, l, ju)
!integer                                ::            n, l, ju
!integer                                ::            i, j
!real*8, dimension(0:n-1)               ::            id
!real*8, dimension(0:n-1)               ::            tbr, tbi, tor, toi
!real*8, dimension(0:n-1)               ::            fst
!real*8, parameter                      ::            pi=3.1415926
!real*8                                 ::            dr, dk, th
!
!dr=0.05
!dk=(pi/dble(n))/dr
!
!tbr(0)=0.0
!tbi=0.0
!do i=1, n-1
!    tbr(i)=sin(pi*dble(i)/dble(n))*(id(i)+id(n-i))+0.5*(id(i)-id(n-i))
!enddo
!
!call fft(l, n, tbr, tbi, tor, toi)
!
!fst(0)=0.0
!fst(1)=0.5*tor(0)
!do j=1, n/2-1
!    i=j+j
!    fst(i)=toi(j)
!    fst(i+1)=fst(i-1)+tor(j)
!enddo
!
!if (ju==1) then
!    th=4.0*dr*pi
!    fst(:)=th*fst(:)
!elseif (ju==-1) then
!    th=dk/(2.0*(pi**2.))
!    fst(:)=th*fst(:)
!endif
!
!end function fst
!!}}}
!function may{{{
function may(d)
    real(8),intent(in)         :: d
    real(8)                    :: may(n)
    integer                    :: ii
    do ii = 1,n
        if(dr(ii) < d)then
             may(ii) = - 1
         else
             may(ii) = 0
         endif
    end do
end function may
!}}}
!function judge{{{
function judge(a,b)
    real(8),intent(in)      :: a(n)
    real(8),intent(in)      :: b(n)
    real(8)                 :: judge
    integer                 :: ii
    
    judge = 0
    do ii = 1,n 
        judge = judge + abs(a(ii) - b(ii))
    end do
    judge = judge/dble(n)
    
end function judge
!}}}
!subroutine evolution{{{
subroutine evolution()
    times = 0
    ckmm = 1.0
    do 
        test = ckmm
        gkmm(1) = 0.0
        do i = 2,n
            gkmm(i) = rhom*ckmm(i)**2.0/(dk(i) - rhom*ckmm(i))
        end do
        ! inverse fft to calculate grmm
        grmm = fst(gkmm,- 1)
        ! calculate crmm with PY closure
        do i = 1,n
            crmm(i) = (dr(i) + grmm(i))*maymm(i)
        end do
        ! fft to calculate ckmm
        ckmm = fst(crmm,1)
        ! judge the convergence
        ! get a mean value of ckmm
        do i = 1,n
            test1(i) = test(i)*gold + ckmm(i)*(1.0 - gold)
            test2(i) = test(i)*(1.0 - gold) + ckmm(i)*gold 
        end do
        
        
        lambda1 = judge(test,test1)
        lambda2 = judge(test,test2)
        if(lambda1 < lambda2)then
            ckmm   = test1
            lambda = lambda1
        else
            ckmm   = test2
            lambda = lambda2
        endif
        times = times + 1

        !if(times > int(1E4))exit
        if(mod(times,100) == 0)then
            print *,"lambda  = ",lambda
        endif
        if(lambda < error)exit
    end do
end subroutine evolution
!}}}
! fft����ӳ���,����:
!                       fft: һά���ٸ���Ҷ�任
!                      ifft: һά���ٸ���Ҷ��任
!                 butterfly: ���ٸ���Ҷ�任������Ҫ�ĵ����㷨
!********************************************************************
!subroutine fft{{{
subroutine fft(l,n,x_real,x_imag,fx_real,fx_imag)
! һά���ٸ���Ҷ�任, nΪ���ݸ�����n=2^l��xΪ�����ʱ������
integer       :: l        !
integer       :: n        ! n = 2^l
integer       :: flag     !
integer       :: ia       !
integer       :: ib       !
! fxΪ�����Ƶ�����У�real��imagΪ��Ӧ��ʵ�����鲿
real*8        :: x_real(n)
real*8        :: x_imag(n)
real*8        :: fx_real(n)
real*8        :: fx_imag(n)
real*8        :: x1(2)
real*8        :: x2(2)
integer       :: i
! ������������˳��Ϊfft����˳��
!         ����������λ��ת��Ϊ2��������ת��Ϊfft����˳��λ�õ�2������
!         ��Ҫע�⣬����ת��λ��˳���0��ʼ��n-1��fortran�洢λ��˳��
!         ���Ǵ�1��ʼ��n
do i=0,n-1
    flag=mod(i,2)*2**(l-1)
    do j=1,l-1
        if(i<2**j) exit
        flag=flag+mod(int(i/(2**j)),2)*(2**(l-j-1))
    end do !j
    !fx_real(n-flag)=x_real(i+1)
    !fx_imag(n-flag)=x_imag(i+1)
    fx_real(flag+1)=x_real(i+1)
    fx_imag(flag+1)=x_imag(i+1)
end do !i
do i=1,l
    flag=2**i
    do j=1,n/flag
        ia=(j-1)*flag+1
        ib=j*flag
        call butterfly(flag,fx_real(ia:ib),fx_imag(ia:ib))
    end do !j
end do !i
end subroutine fft 
!}}}
!subroutine getfftfreq{{{
subroutine getfftfreq(length,freq)
! ��fft���˳����ӦƵ��
real*8    freq(n)
do i=0,n-1
    flag=i
    if(i>n/2) flag=i-n    
    freq(i+1)=flag*2*pi/n
end do !i
end subroutine getfftfreq
!}}}
!subroutine ifft{{{
subroutine ifft(l,n,x_real,x_imag,fx_real,fx_imag)
! һά���ٸ���Ҷ�任����任
integer       ::    l,n 
integer       ::    flag,ia,ib
real*8        ::    x_real(n),x_imag(n),fx_real(n),fx_imag(n)
real*8        ::    y(n)
y=-x_imag
call fft(l,n,x_real,y,fx_real,fx_imag)
fx_real=fx_real/n
fx_imag=-fx_imag/n
end subroutine ifft 
!********************************************************************
!}}}
!subroutine butterfly{{{
subroutine butterfly(length,x_real,x_imag)
! ���ٸ���Ҷ�任��Ҫ�õ��ĵ����㷨
integer        ::    length
real*8        ::    x_real(length),x_imag(length),fx_real(length),fx_imag(length)
real*8        ::    pi=3.1415926535
fx_real=x_real;fx_imag=x_imag
do i=1,length/2
    x_real(i)=fx_real(i)+fx_real(i+length/2)*cos(2*pi*(i-1)/length)&
              -fx_imag(i+length/2)*sin(2*pi*(i-1)/length)
    x_imag(i)=fx_imag(i)+fx_imag(i+length/2)*cos(2*pi*(i-1)/length)&
              +fx_real(i+length/2)*sin(2*pi*(i-1)/length)
    x_real(i+length/2)=fx_real(i)-fx_real(i+length/2)*cos(2*pi*(i-1)&
              /length)+fx_imag(i+length/2)*sin(2*pi*(i-1)/length)
    x_imag(i+length/2)=fx_imag(i)-fx_imag(i+length/2)*cos(2*pi*(i-1)&
              /length)-fx_real(i+length/2)*sin(2*pi*(i-1)/length)
end do !i
!it can also be expressed by the following codes 
!complex x(n)
!complex y(n)
!complex w
!integer k
!do k = 1,n/2
!    w        = exp((0,2*pi*(k-1)/n))
!    x(k)     = y(k) + y(k+n/2)*w
!    x(k+n/2) = y(k) - y(k+n/2)*w
!end do
end subroutine butterfly
!********************************************************************
!}}}

end module module_fst
