module module_fst 
     use module_common
contains

!--****************************************************************--!
!subroutine fst{{{
    ! calculate the fourier sine transform
subroutine fst(od,id,ju)
    ! id  :   input data
    ! od  :  output data
    ! ju     1 for forward transform
    ! ju    -1 for backward transform
real*8, parameter              ::  pi = 3.1415926
integer                        ::  ju
integer                        ::  i
integer                        ::  j
real(8), dimension(0:n-1)      ::  id
real(8), dimension(0:n-1)      ::  od
real(8), dimension(0:n-1)      ::  tbr
real(8), dimension(0:n-1)      ::  tbi
real(8), dimension(0:n-1)      ::  tor
real(8), dimension(0:n-1)      ::  toi
real(8)                        ::  th

tbr(0)=0.0
tbi=0.0
do i=1, n-1
    tbr(i)=sin(pi*dble(i)/dble(n))*(id(i)+id(n-i))+0.5*(id(i)-id(n-i))
enddo

call fft(l, n, tbr, tbi, tor, toi)

od(0)=0.0
od(1)=0.5*tor(0)
do j=1, n/2-1
    i=j+j
    od(i)=toi(j)
    od(i+1)=od(i-1)+tor(j)
enddo

if (ju==1) then
    th=4.0*deltar*pi
    od(:)=th*od(:)
elseif (ju==-1) then
    th=deltak/(2.0*(pi**2.))
    od(:)=th*od(:)
endif

end subroutine fst
!}}}
!subroutine may{{{
    subroutine may(mayer,d)
        real(8)         ::  mayer(n)
        real(8)         ::  d
        real(8)         ::  x
        integer         ::  itmp
        do itmp = 1,n
            x = itmp*deltar
            if(x <= d)then
                mayer(itmp) = - 1
            else
                mayer(itmp) = 0
            endif
        end do
    end subroutine may
!}}}
!subroutine converg{{{
    !  convergence judgement
    subroutine converg(lambda,a,b)
        real(8),intent(out)           :: lambda
        real(8),intent(in)            :: a(n) 
        real(8),intent(in)            :: b(n) 
        integer                       :: itmp
        real(8)                       :: tmp 
        
        tmp = 0
        do i = 1,n
           tmp = tmp + abs(a(i) - b(i))
        end do
        tmp = tmp/n
        lambda = tmp

    end subroutine converg
!}}}
!subroutine evolution{{{
    subroutine evolution()
    do   !stop while convergent 
            test = ckmm

         hkmm(1) = 0D0  
        hkffc(1) = 0D0 
         hkfm(1) = 0D0 
        hkffb(1) = 0D0 
         gkmm(1) = 0D0 
        gkffc(1) = 0D0 
         gkfm(1) = 0D0 
        gkffb(1) = 0D0 
        do i = 2,n
           !************myself formula**********************
            hkmm(i) = dk(i)*ckmm(i)/(dk(i) - rhom*ckmm(i))
            hkffc(i)= dk(i)*ckffc(i)/(dk(i) - rhom*ckffc(i))
            hkfm(i) = (dk(i)*ckfm(i) + rhom*ckfm(i)*hkmm(i))/&
                      (dk(i) - rhof*ckffc(i))
            hkffb(i)= (dk(i)*ckffb(i) + rhom*ckfm(i)*hkfm(i)  &
                      + rhof*ckffb(i)*hkffc(i))/(dk(i) - rhof*ckffc(i))

            gkmm(i) = hkmm(i)  - ckmm(i)   
            gkffc(i)= hkffc(i) - ckffc(i)
            gkfm(i) = hkfm(i)  - ckfm(i) 
            gkffb(i)= hkffb(i) - ckffb(i)
        end do   ! i
        !!print *,hkmm
        !!pause
        !  inverse fft for gamma
        call fst(grmm ,gkmm ,-1)
        call fst(grffc,gkffc,-1)
        call fst(grfm ,gkfm ,-1)
        call fst(grffb,gkffb,-1)
        !!print *,grmm
        !!pause
        !  for crffb
        do i = 2,n
          if(dr(i) < dmm) then
              crffb(i) =  dr(i) + crmm(i) + grmm(i)
          else
              crffb(i) = crmm(i)
          endif
        end do
         crmm(1) = 0 
         crfm(1) = 0 
        crffc(1) = 0
        crffb(1) = 0
        do i = 2,n
            crmm(i)  = (dr(i)+grmm(i))*maymm(i)
            crfm(i)  = (dr(i)+grfm(i))*mayfm(i)
            crffc(i) = (dr(i) + grffb(i) + grffc(i))&
                *mayff(i) - crffb(i)
        end do
        !  fft to calculate c(k)
        call fst(ckmm ,crmm ,1)
        call fst(ckffc,crffc,1)
        call fst(ckfm ,crfm ,1)
        call fst(ckffb,crffb,1)
        !  test the convergence
        call converg(lambda,test,ckmm)
        !  to abort the cycle
        times = times + 1
        if(times == fre)then
            call cpu_time(t2)
            print *,"the time is ",t2 - t1
        endif
        if(lambda < error) exit
        if(times > int(2E7))exit
        !if(mod(times,fre) == 0) print *,lambda
    enddo 
    end subroutine evolution
!}}}
!********************************************************************
! fft相关子程序,包括:
!                 fft: 一维快速富里叶变换
!                   ifft: 一维快速富里叶逆变换
!                 butterfly: 快速富里叶变换中所需要的蝶形算法
!********************************************************************
!subroutine fft{{{
subroutine fft(l,n,x_real,x_imag,fx_real,fx_imag)
! 一维快速富里叶变换, n为数据个数，n=2^l，x为输入的时间序列
! fx为输出的频谱序列，real和imag为相应的实部和虚部
integer        ::    l,n,flag,ia,ib
real*8        ::    x_real(n),x_imag(n),fx_real(n),fx_imag(n),x1(2),x2(2)
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
subroutine getfftfreq(n,freq)
! 求fft输出顺序相应频率
real*8    freq(n)
do i=0,n-1
    flag=i
    if(i>n/2) flag=i-n    
    freq(i+1)=flag*2*3.1415926535/n
end do !i
end subroutine getfftfreq
!}}}
!subroutine ifft{{{
subroutine ifft(l,n,x_real,x_imag,fx_real,fx_imag)
! 一维快速富里叶变换的逆变换
integer        ::    l,n,flag,ia,ib
real*8        ::    x_real(n),x_imag(n),fx_real(n),fx_imag(n)
real*8        ::    y(n)
y=-x_imag
call fft(l,n,x_real,y,fx_real,fx_imag)
fx_real=fx_real/n
fx_imag=-fx_imag/n
end subroutine ifft 
!}}}
!subroutine butterfly{{{
subroutine butterfly(n,x_real,x_imag)
! 快速富里叶变换需要用到的蝶形算法
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
!}}}

end module module_fst
