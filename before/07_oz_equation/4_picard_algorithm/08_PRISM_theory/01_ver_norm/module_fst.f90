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
real*8                        :: dr, dk, th
!%%%%%%%%%%%%%%%%%%NOTICE%%%%%%%%%%%%%%%%%%%%!
!    dr should be equal to deltar
!%%%%%%%%%%%%%%%%%%NOTICE%%%%%%%%%%%%%%%%%%%%!
dr=deltar
dk=(pi/dble(n))/dr

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
    th=4.0*dr*pi
    fst(:)=th*fst(:)
elseif (ju==-1) then
    th=dk/(2.0*(pi**2.))
    fst(:)=th*fst(:)
endif

end function fst
!}}}
! fft相关子程序,包括:
!                       fft: 一维快速傅立叶变换
!                      ifft: 一维快速傅立叶逆变换
!                 butterfly: 快速傅立叶变换中所需要的蝶形算法
!********************************************************************
!subroutine fft{{{
subroutine fft(l,n,x_real,x_imag,fx_real,fx_imag)
! 一维快速傅立叶变换, n为数据个数，n=2^l，x为输入的时间序列
integer       :: l        !
integer       :: n        ! n = 2^l
integer       :: flag     !
integer       :: ia       !
integer       :: ib       !
! fx为输出的频谱序列，real和imag为相应的实部和虚部
real*8        :: x_real(n)
real*8        :: x_imag(n)
real*8        :: fx_real(n)
real*8        :: fx_imag(n)
real*8        :: x1(2)
real*8        :: x2(2)
integer       :: i
! 重排输入数据顺序为fft输入顺序：
!         将输入数据位置转换为2进制数翻转即为fft输入顺序位置的2进制数
!         需要注意，上述转换位置顺序从0开始到n-1，fortran存储位置顺序
!         则是从1开始到n
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
! 求fft输出顺序相应频率
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
! 一维快速傅立叶变换的逆变换
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
! 快速傅立叶变换需要用到的蝶形算法
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
