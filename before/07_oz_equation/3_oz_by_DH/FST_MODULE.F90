module FST_MODULE
contains

!--****************************************************************--!

FUNCTION FST(ID, N, L, JU)
INTEGER								::			N, L, JU
INTEGER								::			I, J
REAL*8, DIMENSION(0:N-1)			::			ID
REAL*8, DIMENSION(0:N-1)			::			TBR, TBI, TOR, TOI
REAL*8, DIMENSION(0:N-1)			::			FST
REAL*8, PARAMETER					::			PI=3.1415926
REAL*8								::			DR, DK, TH

DR=0.05; DK=(PI/DBLE(N))/DR

TBR(0)=0.0; TBI=0.0
DO I=1, N-1
	TBR(I)=SIN(PI*DBLE(I)/DBLE(N))*(ID(I)+ID(N-I))+0.5*(ID(I)-ID(N-I))
ENDDO

CALL FFT(L, N, TBR, TBI, TOR, TOI)

FST(0)=0.0; FST(1)=0.5*TOR(0)
DO J=1, N/2-1
	I=J+J
	FST(I)=TOI(J)
	FST(I+1)=FST(I-1)+TOR(J)
ENDDO

IF (JU==1) THEN
	TH=4.0*DR*PI
	FST(:)=TH*FST(:)
ELSEIF (JU==-1) THEN
	TH=DK/(2.0*(PI**2.))
	FST(:)=TH*FST(:)
ENDIF

END FUNCTION FST

!--****************************************************************--!


!********************************************************************
! 伪谱法求各阶导数相关子程序,包括:
! 				aFar1_1D: 一维一阶导数
! 				aFar2_1D: 一维二阶导数
! 				aFar4_1D: 一维四阶导数
! 				aFar_1DgMaM: 1D gMAM中用到的各阶导数
! 				aFar_2DgMaM: 2D gMAM中用到的各阶导数
!********************************************************************
subroutine aFar1(n,x_real,x_imag,freq)
! n ：离散化后的格点数
! x ：采用同址算法
!     输入：原始数据的频域形式，即原始数据经过富里叶变换的结果
! 		real,imag为相应的实部与虚部
! 	  输出：一阶导数的频域形式，需经过富里叶逆变换为一阶导数的时域形式
!     real,imag为相应的实部与虚部
! freq: x相应位置的频率值
integer		:: n
real*8		:: x_real(n),x_imag(n),freq(n),f_real,f_imag	
do	i=1,n
	f_real=x_real(i)
	f_imag=x_imag(i)
	x_real(i)=f_imag*freq(i)
	x_imag(i)=-f_real*freq(i)
end do !i
end subroutine aFar1
!********************************************************************
subroutine aFar2(n,x_real,x_imag,freq)
! n ：离散化后的格点数
! x ：采用同址算法
!     输入：原始数据的频域形式，即原始数据经过富里叶变换的结果
! 		real,imag为相应的实部与虚部
! 	  输出：二阶导数的频域形式，需经过富里叶逆变换为二阶导数的时域形式
!     real,imag为相应的实部与虚部
! freq: x相应位置的频率值
integer		:: n
real*8		:: x_real(n),x_imag(n),freq(n)	
do	i=1,n
	x_real(i)=-x_real(i)*freq(i)**2
	x_imag(i)=-x_imag(i)*freq(i)**2
end do !i
end subroutine aFar2
!********************************************************************
subroutine aFar4(n,x_real,x_imag,freq)
! n ：离散化后的格点数
! x ：采用同址算法
!     输入：原始数据的频域形式，即原始数据经过富里叶变换的结果
! 		real,imag为相应的实部与虚部
! 	  输出：四阶导数的频域形式，需经过富里叶逆变换为四阶导数的时域形式
!     real,imag为相应的实部与虚部
! freq: x相应位置的频率值
integer		:: n
real*8		:: x_real(n),x_imag(n),freq(n)
do	i=1,n
	x_real(i)=x_real(i)*freq(i)**4
	x_imag(i)=x_imag(i)*freq(i)**4
end do !i
end subroutine aFar4	
!********************************************************************
! FFT相关子程序,包括:
! 				fft: 一维快速富里叶变换
! 	  			ifft: 一维快速富里叶逆变换
! 				butterfly: 快速富里叶变换中所需要的蝶形算法
!********************************************************************
subroutine fft(l,n,x_real,x_imag,fx_real,fx_imag)
! 一维快速富里叶变换, n为数据个数，n=2^l，x为输入的时间序列
! fx为输出的频谱序列，real和imag为相应的实部和虚部
integer		::	l,n,flag,ia,ib
real*8		::	x_real(n),x_imag(n),fx_real(n),fx_imag(n),x1(2),x2(2)
! 重排输入数据顺序为fft输入顺序：
! 		将输入数据位置转换为2进制数翻转即为fft输入顺序位置的2进制数
! 		需要注意，上述转换位置顺序从0开始到n-1，fortran存储位置顺序
! 		则是从1开始到n
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
subroutine getFFTfreq(n,freq)
! 求FFT输出顺序相应频率
real*8	freq(n)
do i=0,n-1
	flag=i
	if(i>n/2) flag=i-n    
	freq(i+1)=flag*2*3.1415926535/n
end do !i
end subroutine getFFTfreq
!********************************************************************
subroutine ifft(l,n,x_real,x_imag,fx_real,fx_imag)
! 一维快速富里叶变换的逆变换
integer		::	l,n,flag,ia,ib
real*8		::	x_real(n),x_imag(n),fx_real(n),fx_imag(n)
real*8		::	y(n)
y=-x_imag
call fft(l,n,x_real,y,fx_real,fx_imag)
fx_real=fx_real/n
fx_imag=-fx_imag/n
end subroutine ifft 
!********************************************************************
subroutine butterfly(n,x_real,x_imag)
! 快速富里叶变换需要用到的蝶形算法
integer		::	n
real*8		::	x_real(n),x_imag(n),fx_real(n),fx_imag(n)
real*8		::	PI=3.1415926535
fx_real=x_real;fx_imag=x_imag
do i=1,n/2
	x_real(i)=fx_real(i)+fx_real(i+n/2)*cos(2*PI*(i-1)/n)-fx_imag(i+n/2)*sin(2*PI*(i-1)/n)
	x_imag(i)=fx_imag(i)+fx_imag(i+n/2)*cos(2*PI*(i-1)/n)+fx_real(i+n/2)*sin(2*PI*(i-1)/n)
	x_real(i+n/2)=fx_real(i)-fx_real(i+n/2)*cos(2*PI*(i-1)/n)+fx_imag(i+n/2)*sin(2*PI*(i-1)/n)
	x_imag(i+n/2)=fx_imag(i)-fx_imag(i+n/2)*cos(2*PI*(i-1)/n)-fx_real(i+n/2)*sin(2*PI*(i-1)/n)
end do !i
end subroutine butterfly
!********************************************************************

end module FST_MODULE
