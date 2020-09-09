program main
    use module_graph
   
    
    implicit none
    integer,parameter              :: ux1 = ux+k*m,uy1=uy+k*n
    call initgraphwindow(mycolor(17))
    !call setgraphwindow(ux,uy,dx,dy)
    call RANDOM_SEED()
    do i = 1,n
    do j = 1,m
    call RANDOM_NUMBER(x1)
    state(i,j) = int(2*x1)
    end do !j
   end do !i
   x = 400
   y = 400
  call gettriarray(tri)
  call gethexaarray(hexa)
  call trigonal()
  x0 = tri(1,1,1)-r2
  y0 = tri(1,1,2)-r2-sqrt(3.)*r2
  call hexagonal(x0,y0)
 
  open(10,file = "data.txt")
  do i = 1,n-1
  do j = 1,m-1
  write(10,*)tri(i,j,:)
  end do !j
  end do !i
  close(10)
 
    end program main