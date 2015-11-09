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

call fft(tbr, tbi, tor, toi)

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
!subroutine inverse1{{{
subroutine inverse1(a,n,l)
    integer,intent(in)    :: l
    integer,intent(in)    :: n
    integer,intent(out)   :: a(n)
    integer               :: b(n)
    integer               :: c(n)
    integer               :: i

    b(1) = 2
    do i = 2,n
        b(i) = b(i-1)*2
    end do
    do i = 1,n-1
        c(i) = b(l-i)
    end do
    c(n) = 1
    a(1) = 1
    a(n) = b(l) 
    do i = 1,n/2 -1
        if(mod(i,2) == 1)then
            a(i+1) = a(i) + b(l-1)
        else
            num = i - 1
            do j = 1,l-1
               k = mod(num,2)
               if(k == 0) exit
               num = num/2
            end do
            j = j-1
            !print *,"j = ",j
            a(i+1) = a(i) - sum(c(1:j)) + c(j+1)
        endif
        j = i + n/2
        a(j) =a(i) + 1
    end do
end subroutine inverse1
!}}}
!******************************************************
!  fft for fast fourier transform
! n for the length of the transform
! where n is 2**l
!subroutine fft{{{
subroutine fft(x_re,x_im,y_re,y_im)
    real(8),intent(in)    :: x_re(n)
    real(8),intent(in)    :: x_im(n)
    real(8),intent(out)   :: y_re(n)
    real(8),intent(out)   :: y_im(n)
    integer               :: a(n)    ! original order
    integer               :: b(l)    ! 2 power 
    integer               :: c(l)    ! n/b(l) 
    integer               :: ia 
    integer               :: ib 
    integer               :: k 
    integer               :: num 
    integer               :: cy

    b(1) = 2
    do i = 1,l-1
       b(i+1) = b(i)*2 
    end do
    do i = 1,l-1
       c(i) = b(l - i)
    end do
    c(n) = 1
    ! calculate a(n)
    call inverse1(a,n,l)

    !  fft
    do i = 1,l
        num = b(i)
        do j = 1,c(i) 
            ia = (j-1)*num
            ib = (j-1)*num + num/2
            call butterfly(y_re(ia:ib),y_im(ia:ib),num)
        end do
    end do
end subroutine fft
!}}}
!  butterfly for butterfly algorithm
!subroutine butterfly{{{
subroutine butterfly(y_re,y_im,n)
    integer,intent(in)       :: n
    real(8),intent(inout)    :: y_re(n)
    real(8),intent(inout)    :: y_im(n)
    real(8)                  :: x_re(n)
    real(8)                  :: x_im(n)
    real(8)                  :: w(n)
    
    x_re = y_re
    x_im = y_im
    do i = 1,n/2
        j = i + n/2 
        w(i) = cos(2*pi*(i-1)/dble(n))
        w(j) = sin(2*pi*(i-1)/dble(n))
        y_re(i) = x_re(i) + w(i)*x_re(j) + w(j)*x_im(j)
        y_im(i) = x_im(i) + w(i)*x_im(j) - w(j)*x_re(j)
        y_re(i) = x_re(i) - w(i)*x_re(j) - w(j)*x_im(j)
        y_im(i) = x_im(i) - w(i)*x_im(j) + w(j)*x_re(j)
    end do
end subroutine butterfly
!}}}
end module module_fst
