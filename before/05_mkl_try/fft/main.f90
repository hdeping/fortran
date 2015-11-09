 program main
     ! fortran example. 
     ! compile line: ! ifort -mkl
     ! /opt/intel/compiler/11.1/069/mkl/include/mkl_dfti.f90 -o3 tfft.f90 -o tfft
     ! 1d complex to complex, and real to conjugate-even 
     use mkl_dfti 
      integer,parameter              :: mysize  = 32
      integer,parameter              :: mysize2 = 514
      integer                        :: count_rate
      integer                        :: count_max 
      integer                        :: ic 
      integer                        :: status 
      complex(8)                     :: x(mysize)
      complex(8)                     :: y(mysize)
      real                           :: x_re(mysize)
      real                           :: x_im(mysize)
      real                           :: y_re(mysize)
      real                           :: y_im(mysize)
      real                           :: t1
      real                           :: t2
      type(dfti_descriptor), pointer :: my_desc1_handle
      type(dfti_descriptor), pointer :: my_desc2_handle
      character(20)                  :: filename

      filename = "data.txt"
      open(10,file = filename)
      !...put input data into x(1),...,x(mysize); y(1),...,y(mysize) 
      do i=1,mysize 
          x(i)=i
      end do
      ! perform a complex to complex transform 
      call cpu_time(t1)
      write(*,*)"before fft"
      do i=1,mysize
          write(*,*)x(i)
      end do 
      y = x
      status = dfticreatedescriptor( my_desc1_handle, dfti_single, &
               dfti_single, 1, mysize ) 
      status = dfticommitdescriptor( my_desc1_handle ) 
      status = dfticomputeforward( my_desc1_handle,y) 
      status = dftifreedescriptor(my_desc1_handle) 
      write(*,*)"forward fft"
      do i=1,mysize
          write(*,*)y(i)
      end do 
      !  perform a backforward fft
      status = dfticreatedescriptor( my_desc1_handle, dfti_single, &
               dfti_single, 1, mysize ) 
      status = dfticommitdescriptor( my_desc1_handle ) 
      status = dfticomputebackward( my_desc1_handle,y) 
      status = dftifreedescriptor(my_desc1_handle) 
      write(*,*)"backward fft"
      do i=1,mysize
          write(*,*)x_re(i),x_im(i)
      end do 
      call cpu_time(t2)
      print *,"the time is",t2 - t1
  end program main
