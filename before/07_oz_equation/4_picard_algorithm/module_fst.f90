module module_fst
    use mkl_dfti
    use module_common
    !  variables for fft
    integer                        :: status 
    type(dfti_descriptor), pointer :: my_desc1_handle
    type(dfti_descriptor), pointer :: my_desc2_handle
    
    contains
!subroutine fft{{{
    subroutine fft(x,y)
       complex(8),intent(in)       :: x(n)
       complex(8),intent(out)      :: y(n)

       y = x
        status = dfticreatedescriptor( my_desc1_handle, dfti_single, &
                 dfti_real, 1, n) 
        status = dfticommitdescriptor( my_desc1_handle ) 
        status = dfticomputeforward( my_desc1_handle, y) 
        status = dftifreedescriptor(my_desc1_handle) 
    end subroutine fft
!}}}
!subroutine ifft{{{
    subroutine ifft(x,y)
       complex(8),intent(in)       :: x(n)
       complex(8),intent(out)      :: y(n)

       y = x
        status = dfticreatedescriptor( my_desc1_handle, dfti_single, &
                 dfti_real, 1, n) 
        status = dfticommitdescriptor( my_desc1_handle ) 
        status = dfticomputebackward( my_desc1_handle, y) 
        status = dftifreedescriptor(my_desc1_handle) 
    end subroutine ifft
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
!!subroutine fst{{{
!subroutine fst(od,id,ju)
!    ! id  :   input data
!    ! od  :  output data
!    ! ju     1 for forward transform
!    ! ju    -1 for backward transform
!integer                        ::  ju
!integer                        ::  i
!integer                        ::  j
!real(8), dimension(0:n-1)      ::  id
!real(8), dimension(0:n-1)      ::  od
!real(8), dimension(0:n-1)      ::  tbr
!real(8), dimension(0:n-1)      ::  tbi
!real(8), dimension(0:n-1)      ::  tor
!real(8), dimension(0:n-1)      ::  toi
!real(8)                        ::  th
!tbr(0)=0.0
!tbi=0.0
!do i=1, n-1
!    tbr(i)=sin(pi*dble(i)/dble(n))*(id(i)+id(n-i))+0.5*(id(i)-id(n-i))
!enddo
!
!call fft(tbr, tbi, tor, toi)
!
!od(0)=0.0
!od(1)=0.5*tor(0)
!do j=1, n/2-1
!    i=j+j
!    od(i)=toi(j)
!    od(i+1)=od(i-1)+tor(j)
!enddo
!
!if (ju==1) then
!    th=4.0*deltar*pi
!    od(:)=th*od(:)
!elseif (ju==-1) then
!    th=deltak/(2.0*(pi**2.))
!    od(:)=th*od(:)
!endif
!
!end subroutine fst
!!}}}
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

end module module_fst
