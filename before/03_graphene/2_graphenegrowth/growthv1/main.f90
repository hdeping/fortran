program main
    use module_graph
    use module_lattice
    use module_growth
   
    implicit none
    integer   num,output(m)
    integer   array1(4),array2(8),array3(4,2)
    integer  times
    
    call RANDOM_SEED()
    call initgraphwindow(mycolor(2))
    call setviewport(ux,uy,dx,dy)
    ! initialize the interface
    state = 0
    array1 = (/1,2,m-1,m/)
    !do i = 1,4
    !    array2(i) = i
    !    array2(i+4) = n+1-i
    !    state(:,array1(i)) = 1
    !end do  !i
    !do i = 1,8
    !    state(array2(i),:) = 1
    !end do  !i
    do i = 4,n-10,20
        do j = 8,m-10,20
            array3(:,1) = (/i,i,i+10,i+10/)
            array3(:,2) = (/j,j+7,j-3,j+10/)
            do itmp = 1,4
                state(array3(itmp,1),array3(itmp,2)) = 1
                call neighbor2(nei2,array3(itmp,1),array3(itmp,2))
                do jtmp = 1,6
                    state(nei2(jtmp,1),nei2(jtmp,2)) = 1
                end do  ! jtmp
            end do   !itmp
        end do  !j
    end do   ! i
    call latticeinterface(n,m,state)
    
    
    end program main
    !call boundary()
    !jcolor = saveimage(filename,0,0,1440,900)
   !  initial the border array
   
   !do i = 1,25
      ! do j = 23+i,77-i
         !  state(j,i+2) = 1
         !  call drawcircle(1,j,i+2)
         !   !state(i+2,j) = 1
         !  !call drawcircle(1,i+2,j)
      ! end do  !j
   !end do   ! i
   
    !do  jtmp = 1,10
 !  call boundary()
 !  times = 0
 !  p = >head
 !  call RANDOM_NUMBER(x1)
 !  junit = int(1000*x1)+1
 !  do while(associated(p))
    !   times = times+1
    !   if(times>junit.and.times<junit+100)then
    !   i = p.i
    !   j = p.j
    !   state(i,j) = 1
    !   call drawcircle(state(i,j),i,j)
    !   endif
    !   p = >p.next
 !  end do  ! while
 !  end do  ! jtmp