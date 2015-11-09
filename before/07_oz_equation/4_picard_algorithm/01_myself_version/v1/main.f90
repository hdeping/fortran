program main
    !use module_fst 
    use module_common
    use fst_module
    integer  itmp

    real(8)    :: a(n)
    real(8)    :: b(n)
    real(8)    :: c(n)
    !print *,10
    

    filename = "data_g_rmm.txt"
    open(10,file = filename)
    filename = "data_g_rfm.txt"
    open(20,file = filename)
    filename = "data_g_rffb.txt"
    open(30,file = filename)
    filename = "data_g_rff.txt"
    open(40,file = filename)
    filename = "run.log"
    open(50,file = filename)
    ! initial the value of c(four arrays)

      ckmm = 1.0  
      ckfm = 1.0
     ckffb = 1.0
     ckff  = 1.0
    !call random_number(ckmm(:) )
    !call random_number(ckfm(:) )
    !call random_number(ckffb(:))
    !call random_number(ckffc(:))
    !!print *,ckmm
    !!pause

    j = 0
    times = 0

    !  initial the k and r
        dk(1) =  0.0
        dr(1) =  0.0
    do i = 2,n
        dk(i) = (i - 1)*deltak
        dr(i) = (i - 1)*deltar
    end do
    !  call the mayer function
    !call may(maymm,dmm)
    !call may(mayfm,dfm)
    !call may(mayff,dff)
!********************************************************** 

    !  solve the equation
    call cpu_time(t1)

    !call evolution()

    call cpu_time(t2)

    do i = 1,n
        a(i) = (i - 1)*1.0
    end do
    write(*,*)"before forward backward fft"
    b = fst(a,n,l,1)
    c = fst(b,n,l,-1)
    do i = 1,n
        write(*,"(i3,3f18.9)")i,a(i),b(i),c(i)
    end do
     



!********************************************************** 

!    print *,"lambda is ",lambda
!    !  write ck to the file
!    crff(:) = crffb(:) + crffc(:)
!    grff(:) = grffb(:) + grffc(:)
!    do i = 1,n
!       g_rmm(i) = (crmm(i)  +  grmm(i) )/dr(i)  + 1
!       g_rfm(i) = (crfm(i)  +  grfm(i) )/dr(i)  + 1
!       g_rff(i)=  (crff(i)  +  grff(i) )/dr(i)  + 1
!    end do   !  i
!    do i = 1,n
!        write(10,*)dr(i),g_rmm(i) 
!        write(20,*)dr(i),g_rfm(i) 
!        write(30,*)dr(i),g_rffb(i) 
!        write(40,*)dr(i),g_rff(i)
!    end do
!    write(50,*)"times  = ",times
!    write(50,*)"lambda = ",lambda

    close(50)
    close(40)
    close(30)
    close(20)
    close(10)
    !print *,hkmm(10)
end program main
