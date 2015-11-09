module my_own_fft
    integer,parameter              :: l = 7
    integer,parameter              :: n = 2**l
    real(8),parameter              :: pi = 3.141592653
    contains
    integer function num(n)
            integer,intent(in)     :: n 
            integer                :: ntmp
            integer                :: itmp
             
            num  = 0 
            ntmp = n
            do
                itmp = mod(ntmp,2)
                ntmp = ntmp/2
                if(itmp == 0)exit
                num  = num + 1
            end do  ! itmp
          
        end function num
!subroutine invertbite{{{
!    calculate the inversebite number
    subroutine invertbite(a)
        integer,intent(out)        :: a(n)
        integer                    :: i
        integer                    :: array_l(l)
        integer                    :: tmp
        
        ! initial the array_l
        do i = 1,l
            array_l(i) = 2**(l - i)    
        end do

        a(1) = 0
        a(n) = n - 1
        do i = 2,n-1
            tmp = num(a(i - 1))
        end do
    end subroutine invertbite
!}}}
!subroutine invnum{{{
!   calculate the fft  serial number
    subroutine invnum(n,m)
        integer,intent(in)         :: m
        integer,intent(out)        :: n
        integer                    :: i
        integer                    :: tmpn
        integer                    :: tmpm
       
        tmpm = m
        tmpn = n
        do while(tmpm /= 0)
           i = mod(tmpm,2)
           tmpn = 2*tmpn + i
           tmpm = tmpm/2
        end do
        n = tmpn
    end subroutine invnum
!}}}
!subroutine bite{{{
!    calculate the inverse number
    subroutine bite(a)
        integer,intent(out)        :: a(n)
        integer                    :: i
        integer                    :: tmp
        do i = 1,n
            call invnum(a(i),i - 1)
            a(i) = a(i) + 1 
        end do
    end subroutine bite
!}}}
!subroutine fft{{{
!  calculate the fast fourier transformation
    subroutine fft(y,x)
        complex(8),intent(in)      :: x(n)
        complex(8),intent(out)     :: y(n) 
        complex(8)                 :: xtmp(n)
        complex(8)                 :: ytmp(n)
        complex(8)                 :: w(n)
        complex(8)                 :: comtmp 
        real                       :: new
        integer                    :: i
        integer                    :: j
        integer                    :: k
        integer                    :: a(n)
        integer                    :: ist    ! start number
        integer                    :: ien    ! end   number
        integer                    :: ntmp   ! temperary length 

        call bite(a)
        do i = 1,n
            xtmp(i) =  x(a(i)) 
        end do
        do i = 1,l
            new    =  - pi*2.0**(1 - 1.0*i)
            comtmp = (0,new)
            w(i) = exp(comtmp)
        end do
        !  fast fourier transform algorithm
        do i = 1,l 
            ntmp = 2**i
            do j = 1, n/ntmp
                ist  = ntmp*(j - 1)
                do k = ist + 1,ist + ntmp/2
                   ytmp(k) = xtmp(k) +  w(i)*xtmp(k+ntmp/2)
                   ytmp(k + ntmp/2) = xtmp(k) -  w(i)*xtmp(k+ntmp/2)
                end do
            end do
            if(i == l)exit
            xtmp = ytmp
        end do
       
        y = ytmp
        
        
    end subroutine fft
!}}}
!subroutine ifft{{{
!  calculate the inverse fast fourier transformation
    subroutine ifft(y,x)
        complex(8),intent(in)      :: x(n)
        complex(8),intent(out)     :: y(n) 

        call fft(y,x)
        y = y/n
        
    end subroutine ifft
!}}}
end module my_own_fft
