module module_mct
    use module_common
    use module_algebra

contains
!subroutine getmatrixes{{{
subroutine getmatrixes()
    !real(8)             :: tmpvalue(m,m)
    filename = "sk.txt"
    open(20,file = filename, status = "old", iostat = ierror)
    ! read Sk and Ck
    do i = 1,n
        read(20,"(3f18.9)",iostat = ierror)sk(1,1,i),sk(1,2,i),sk(2,2,i)
        sk(2,1,i) = sk(1,2,i)
    end do
    close(20)
    filename = "ck.txt"
    open(20,file = filename, status = "old",iostat = ierror)
    do i = 1,n
        read(20,"(3f18.9)",iostat = ierror)ck(1,1,i),ck(1,2,i),ck(2,2,i)
        ck(2,1,i) = ck(1,2,i)
    end do
    close(20)

    ! get new sk, ck
    do i = 1,m
        do j = i,m
            sk(i,j,:) = sk(i,j,:)/sqrt((xrate(i)*xrate(j))) 
            ck(i,j,:) = ck(i,j,:)/sqrt((xrate(i)*xrate(j)))
        end do
        do j = 1,i - 1
            sk(i,j,:) = sk(j,i,:)
            ck(i,j,:) = ck(j,i,:)
        end do
    end do
    ! get inver_sk(m,m,n)
    deltafun = 0.0
    do i = 1,n
        deltafun(i,i) = 1.0
    end do
    do i = 1,n
        inver_sk(:,:,i) = sol_mat(sk(:,:,i),deltafun,m)
    end do
    ! initial D
    diffu(1,:) = (/1.0,1.0/)
    diffu(2,:) = (/1.0,1.0/)
    !  get K and R
    do i = 1,n
        tmpvalue = 0.0
        do itmp = 1,m
            do jtmp = 1,m
            mat_K(itmp,jtmp,i) = i**2.0*diffu(itmp,jtmp)*v(jtmp)
            end do
        end do
        mat_R(:,:,i) = multi_mat(mat_K(:,:,i),inver_sk(:,:,i),m,m,m)
    end do

end subroutine getmatrixes
!}}}
!get MCT equation{{{
function getmct()
   real(8)             :: getmct
   integer             :: ii 
   integer             :: jj
   real(8)             :: tmpa(m,m)
   real(8)             :: tmpb(m,m)
   real(8)             :: tmpc(m,m)
   real(8)             :: tmpd(m,m)

   ! get f(t1) 
   do ii = 1,m
        do jj = 1,m
            f(ii,jj,:,1) = sk(ii,jj,:)
        end do
   end do
   ! get U
   do q = 1,n
        memory(:,:,q,1) = getmemory(q,1)
        mat_U(:,:,q) = deltafun(:,:) + multi_mat(mat_K(:,:,i),&
                       memory(:,:,i,1),m,m,m)
   end do
   do t = 1,tmnum/2
        do q = 1,n
            !  get memory kernel
            tmpa(ii,jj) = - dt*mat_R(ii,jj,q)
            do ii = 1,t - 1
               tmpa(:,:) = tmpa(:,:) - mat_K(:,:,q)*memory(:,:,q,t+1-ii)&
                           *par_f(:,:,q,ii) 
            end do
            mat_V(:,:,q) = tmpa(:,:)
            par_f(:,:,q,t) = sol_mat(mat_U(:,:,q),mat_V(:,:,q),m)
            f(:,:,q,t+1) = f(:,:,q,t) + par_f(:,:,q,t)
        end do
   end do
   getmct = 1
end function getmct
!}}}
!get memory(q,t){{{
function getmemory(q,t)
   integer,intent(in)  :: q
   integer,intent(in)  :: t
   real(8)             :: getmemory(m,m)
   real(8)             :: pre_factor
   real(8)             :: totalvalue
   integer             :: ii
   integer             :: jj

   ! get matrix mat_A, mat_B, mat_D
   do k = 1,n
        mat_B(:,:,k) = multi_mat(ck(:,:,k),f(:,:,k,t),m,m,m)
        mat_D(:,:,k) = multi_mat(f(:,:,k,t),ck(:,:,k),m,m,m)
        mat_A(:,:,k) = multi_mat(mat_B(:,:,k),ck(:,:,k),m,m,m)
   end do
   totalvalue = rhosum**2.0*h**3.0/(16.0*pi**2.0*q**5.0)
   do ii = 1,m
        do jj = 1,m
            pre_factor = totalvalue/sqrt(xrate(ii)&
                         *xrate(jj))
            tmp = 0
            do k = 1,300
                !do p = abs(q - k),q + k
                do p = abs(q - k),300   !  
                    if(p == 0)cycle
                    l1 = dble(q**2.0 + k**2.0 - p**2.0)
                    l2 = dble(p**2.0 + q**2.0 - k**2.0)
                    tmp = tmp + p*k*l1*(l1*mat_A(ii,jj,k)&
                          *f(ii,jj,p,t) + l2*mat_B(ii,jj,k)&
                          *mat_D(ii,jj,p))
                end do
            end do
            getmemory(ii,jj) = tmp*pre_factor
        end do
        !print *,"h = ",h
   end do
end function getmemory
!}}}
!**********************************************
!  fast sine transform and fourier transform **
!**********************************************
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
    dr=0.04
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
    integer       ::    length
    real*8        ::    x_real(length)
    real*8        ::    x_imag(length)
    real*8        ::    fx_real(length)
    real*8        ::    fx_imag(length)
    real*8        ::    pi = 3.1415926535
    fx_real=x_real;fx_imag=x_imag
    do i = 1,length/2
        x_real(i) = fx_real(i) + fx_real(i + length/2)&
                    *cos(2*pi*(i - 1)/length) - fx_imag(i + length/2)&
                    *sin(2*pi*(i - 1)/length)
        x_imag(i) = fx_imag(i) + fx_imag(i + length/2)&
                    *cos(2*pi*(i - 1)/length) + fx_real(i + length/2)&
                    *sin(2*pi*(i - 1)/length)
    
        x_real(i + length/2) = fx_real(i) - fx_real(i + length/2)&
                               *cos(2*pi*(i - 1)/length) + &
                               fx_imag(i + length/2)*sin(2*pi*(i - 1)/length)
        x_imag(i + length/2) = fx_imag(i) - fx_imag(i + length/2)&
                               *cos(2*pi*(i - 1)/length) - &
                               fx_real(i + length/2)*sin(2*pi*(i - 1)/length)
    end do !i
end subroutine butterfly
!********************************************************************
!}}}

end module module_mct
