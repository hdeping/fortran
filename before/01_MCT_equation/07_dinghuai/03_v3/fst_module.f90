module fst_module

contains

!--****************************************************************--!

function fst(id, n, l, ju)
integer                                ::            n, l, ju
integer                                ::            i, j
real*8, dimension(0:n-1)            ::            id
real*8, dimension(0:n-1)            ::            tbr, tbi, tor, toi
real*8, dimension(0:n-1)            ::            fst
real*8, parameter                    ::            pi=3.1415926
real*8                                ::            dr, dk, th

dr=0.01d0; dk=(pi/dble(n))/dr

tbr(0)=0.0; tbi=0.0
do i=1, n-1
    tbr(i)=sin(pi*dble(i)/dble(n))*(id(i)+id(n-i))+0.5*(id(i)-id(n-i))
enddo

call fft(l, n, tbr, tbi, tor, toi)

fst(0)=0.0; fst(1)=0.5*tor(0)
do j=1, n/2-1
    i=j+j
    fst(i)=toi(j)
    fst(i+1)=fst(i-1)+tor(j)
enddo

if (ju==1) then
    th=4.0*dr*pi
    fst(:)=th*fst(:)
elseif (ju==-1) then
    th=dk/(2.0*(pi**2.))
    fst(:)=th*fst(:)
endif

end function fst

!--****************************************************************--!


!********************************************************************
! α�׷�����׵�������ӳ���,����:
!                 afar1_1d: һάһ�׵���
!                 afar2_1d: һά���׵���
!                 afar4_1d: һά�Ľ׵���
!                 afar_1dgmam: 1d gmam���õ��ĸ��׵���
!                 afar_2dgmam: 2d gmam���õ��ĸ��׵���
!********************************************************************
subroutine afar1(n,x_real,x_imag,freq)
! n ����ɢ����ĸ����
! x ������ַͬ�㷨
!     ���룺ԭʼ���ݵ�Ƶ����ʽ����ԭʼ���ݾ�������Ҷ�任�Ľ��
!         real,imagΪ��Ӧ��ʵ�����鲿
!       �����һ�׵�����Ƶ����ʽ���辭������Ҷ��任Ϊһ�׵�����ʱ����ʽ
!     real,imagΪ��Ӧ��ʵ�����鲿
! freq: x��Ӧλ�õ�Ƶ��ֵ
integer        :: n
real*8        :: x_real(n),x_imag(n),freq(n),f_real,f_imag    
do    i=1,n
    f_real=x_real(i)
    f_imag=x_imag(i)
    x_real(i)=f_imag*freq(i)
    x_imag(i)=-f_real*freq(i)
end do !i
end subroutine afar1
!********************************************************************
subroutine afar2(n,x_real,x_imag,freq)
! n ����ɢ����ĸ����
! x ������ַͬ�㷨
!     ���룺ԭʼ���ݵ�Ƶ����ʽ����ԭʼ���ݾ�������Ҷ�任�Ľ��
!         real,imagΪ��Ӧ��ʵ�����鲿
!       ��������׵�����Ƶ����ʽ���辭������Ҷ��任Ϊ���׵�����ʱ����ʽ
!     real,imagΪ��Ӧ��ʵ�����鲿
! freq: x��Ӧλ�õ�Ƶ��ֵ
integer        :: n
real*8        :: x_real(n),x_imag(n),freq(n)    
do    i=1,n
    x_real(i)=-x_real(i)*freq(i)**2
    x_imag(i)=-x_imag(i)*freq(i)**2
end do !i
end subroutine afar2
!********************************************************************
subroutine afar4(n,x_real,x_imag,freq)
! n ����ɢ����ĸ����
! x ������ַͬ�㷨
!     ���룺ԭʼ���ݵ�Ƶ����ʽ����ԭʼ���ݾ�������Ҷ�任�Ľ��
!         real,imagΪ��Ӧ��ʵ�����鲿
!       ������Ľ׵�����Ƶ����ʽ���辭������Ҷ��任Ϊ�Ľ׵�����ʱ����ʽ
!     real,imagΪ��Ӧ��ʵ�����鲿
! freq: x��Ӧλ�õ�Ƶ��ֵ
integer        :: n
real*8        :: x_real(n),x_imag(n),freq(n)
do    i=1,n
    x_real(i)=x_real(i)*freq(i)**4
    x_imag(i)=x_imag(i)*freq(i)**4
end do !i
end subroutine afar4    
!********************************************************************
! fft����ӳ���,����:
!                 fft: һά���ٸ���Ҷ�任
!                   ifft: һά���ٸ���Ҷ��任
!                 butterfly: ���ٸ���Ҷ�任������Ҫ�ĵ����㷨
!********************************************************************
subroutine fft(l,n,x_real,x_imag,fx_real,fx_imag)
! һά���ٸ���Ҷ�任, nΪ���ݸ�����n=2^l��xΪ�����ʱ������
! fxΪ�����Ƶ�����У�real��imagΪ��Ӧ��ʵ�����鲿
integer        ::    l,n,flag,ia,ib
real*8        ::    x_real(n),x_imag(n),fx_real(n),fx_imag(n),x1(2),x2(2)
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
!********************************************************************
subroutine getfftfreq(n,freq)
! ��fft���˳����ӦƵ��
real*8    freq(n)
do i=0,n-1
    flag=i
    if(i>n/2) flag=i-n    
    freq(i+1)=flag*2*3.1415926535/n
end do !i
end subroutine getfftfreq
!********************************************************************
subroutine ifft(l,n,x_real,x_imag,fx_real,fx_imag)
! һά���ٸ���Ҷ�任����任
integer        ::    l,n,flag,ia,ib
real*8        ::    x_real(n),x_imag(n),fx_real(n),fx_imag(n)
real*8        ::    y(n)
y=-x_imag
call fft(l,n,x_real,y,fx_real,fx_imag)
fx_real=fx_real/n
fx_imag=-fx_imag/n
end subroutine ifft 
!********************************************************************
subroutine butterfly(n,x_real,x_imag)
! ���ٸ���Ҷ�任��Ҫ�õ��ĵ����㷨
integer        ::    n
real*8        ::    x_real(n),x_imag(n),fx_real(n),fx_imag(n)
real*8        ::    pi=3.1415926535
fx_real=x_real;fx_imag=x_imag
do i=1,n/2
    x_real(i)=fx_real(i)+fx_real(i+n/2)*cos(2*pi*(i-1)/n)-fx_imag(i+n/2)*sin(2*pi*(i-1)/n)
    x_imag(i)=fx_imag(i)+fx_imag(i+n/2)*cos(2*pi*(i-1)/n)+fx_real(i+n/2)*sin(2*pi*(i-1)/n)
    x_real(i+n/2)=fx_real(i)-fx_real(i+n/2)*cos(2*pi*(i-1)/n)+fx_imag(i+n/2)*sin(2*pi*(i-1)/n)
    x_imag(i+n/2)=fx_imag(i)-fx_imag(i+n/2)*cos(2*pi*(i-1)/n)-fx_real(i+n/2)*sin(2*pi*(i-1)/n)
end do !i
end subroutine butterfly
!********************************************************************

end module fst_module
