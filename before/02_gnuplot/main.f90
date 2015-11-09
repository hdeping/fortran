program main
    use  module_new
    use  module_common

   call random_seed()

   !do j = 1,n
   !    do i = 1,n
   !         x(i) = dble(i)
   !         call random_number(x1) 
   !         y(i) = 2.0*x(i) + x1*dble(j)
   !    end do
   !    write(output,"('new',i2.2)")j
   !    call data_process(x,y,output,n)
   !end do
   do i = 1,n
       x(i) = 0.01*dble(i)
       y(i) = besj0(x(i))
   end do
   output = "bessel"
   call data_process(x,y,output,n)
  
   !do i=1,20
   !    call choose_exmp(i)
   !end do   !  
    
    end program main
     !real*8  x(n),y(n)
   !integer  i,j
   !do i=1,n
      ! x(i)=i
      ! y(i)=i**2
   !end do   !i 
   !call data_line(x,y,n)
   !open(10,file="code.txt")
   !write(10,*)"select case(i)"
   !do i=2,21
      ! write(num,'(i2.2)')i
      ! if(i>10)then
         !  write(num1,"(i2.2,')')")i-1
         !  write(10,*)"case(",trim(num1)
      ! else
         !  write(num1,"(i1.1,')')")i-1
      !     write(10,*)"case(",trim(num1)
      ! end if
      ! write(10,*)"call exmp"//trim(num)//"()"
   !end do   !  i
   !write(10,*)"end select"
   !close(10)
