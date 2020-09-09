module module_lattice 
    use module_common
    use module_graph
    
    contains
!*************************************************
  subroutine statenum(status,num)
    integer  num,tmp,status
    tmp = 0
     do i = 2,n-1
     do j = 2,m-1
     if (state(i,j) = =status)tmp=tmp+1
     end do !j
     end do !i
     num = tmp
     end subroutine statenum   
!*************************************************
   subroutine area_multi(s,n,x,y)
   integer  n
   real s, x(n), y(n)
   real  tmp,tmp1,xtmp(3),ytmp(3)
   tmp = 0
   xtmp(1) = x(1)
   ytmp(1) = y(1)
   do i = 2,n-1
    xtmp(2) = x(i)
   ytmp(2) = y(i)
    xtmp(3) = x(i+1)
   ytmp(3) = y(i+1)
   call area_tri(tmp1,xtmp,ytmp)
   tmp = tmp+tmp1
   end do  !i
   s = tmp
   end subroutine area_multi
!*************************************************
    subroutine area_tri(s,x,y)
    real s,x(3),y(3)
    real tmp(4)
    tmp(1) = x(2)-x(1)
    tmp(2) = y(2)-y(1)
    tmp(3) = x(3)-x(1)
    tmp(4) = y(3)-y(1)
    s = abs(tmp(1)*tmp(4)-tmp(2)*tmp(3))/2.
    end subroutine area_tri
!*************************************************
    ! avoid the most outside border
    subroutine neighbornum(num,i,j)
    integer  num,i,j,times,i1,i2,j2
    times = 0
    call  neighbor2(nei2,i,j)
    do i1 = 1,6
    i2 = nei2(i1,1)
    j2 = nei2(i1,2)
    !if(i2>0.and.i2<n+1.and.j2<m+1.and.j2>0)then
    if(state(i2,j2) = =1)times=times+1                !  the number of the occupied sites
    end do !i1
   num = times
    end subroutine neighbornum
!*************************************************
    subroutine neighbor2(nei2,i,j)
    integer i,j,nei2(6,2)
    nei2(:,1) = (/i+2,i+1,i-1,i-2,i-1,i+1/)
    nei2(:,2) = (/j,j+1,j+1,j,j-1,j-1/)
    end subroutine neighbor2
!*************************************************
    subroutine neinum(num,i,j)
    integer  num,i,j,times
     times = 0
    call  neighbor(nei1,i,j)
    do i1 = 1,3
    i2 = nei1(i1,1)
    j2 = nei1(i1,2)
    if(i2>0.and.i2<n+1.and.j2<m+1.and.j2>0)then
    if(state(i2,j2) = =0)times=times+1
    end if
    end do !i1
   num = times
    end subroutine neinum
!*************************************************
    ! nearest 
    subroutine neighbor(nei1,i,j)
    integer  nei1(4,2),i,j,k
    nei1(1,:) = (/i+1,j/)
    nei1(2,:) = (/i-1,j/)
     if(mod(i,2) = =0)then
     if(mod(j,2) = =0)then
    k = 4
    else
    k = 1
    endif
    else
     if(mod(j,2) = =0)then
    k = 3
    else
    k = 2
    endif
    endif
    select case(k)
    case (2,4)
    nei1(3,:) = (/i,j+1/)    
    nei1(4,:) = (/i,j-1/)  
    case default
    nei1(3,:) = (/i,j-1/)  
    nei1(4,:) = (/i,j+1/)    
    end select
    nei2(:,1) = (/i+2,i+1,i-1,i-2,i-1,i+1/)
    nei2(:,2) = (/j,j+1,j+1,j,j-1,j-1/)
    end subroutine neighbor
!*************************************************
    !   two   rectangle
    subroutine updateborder(x,y)
    integer  x(4),y(4)
    integer  a(3,2)
    integer  i,j,k
    a(:,1) = (/1,2,1/)
    a(:,2) = (/4,3,4/)
    do k = 1,3
    do i = x(k),x(k+1)
        do j = y(a(k,1)),y(a(k,2))
            
        end do  !j
    end do   !i
    end do  !k
    end subroutine updateborder
!*************************************************
    subroutine drawboundary()
    call boundary
     p = >head
   do while(associated(p))
       call drawcircle(2,p.i,p.j)
       p = >p.next
   end do  ! while
    end subroutine drawboundary
!*************************************************
    !  consider the neighbor site and the second neighbor site
    subroutine boundary()
    integer  a,i1,j1
    integer  si,sj,k,icycle
    integer  i_pre,j_pre
    integer  tmpnew(2,3)    ! for the rest neighbors
    integer  times,icount     !  count the cycles
    nullify(head)
    nullify(last)
    open(10,file = "tmp.txt")
    a = n/4
    do 
        times = 0
    do  i = 1,n
        if(state(a,i) = =0)then
            call neinum(num,a,i)
            if(num<3.and.num>0)times = 1
            exit
        endif
    end do  ! i
    if(times = =1)exit
    a = a-1
    end   do   !
    !print *,a
    !read*
    j = i
    i = a
    si = i
    sj = j
    call neighbor(nei1,i,j)
    itmp = 1
    i_pre = nei1(itmp,1)
    j_pre = nei1(itmp,2)
    do 
        call link(i,j)
        call  drawcircle(2,i,j)                       !  drawcircle
        call neighbor(nei1,i,j)
        times = 0
        tmpnew = 0
        do itmp = 1,3
            i1 = nei1(itmp,1)
            j1 = nei1(itmp,2)
            call stateborder(num,i1,j1)                           !  judge if it's a border site
            if((i1 = =i_pre.and.j1==j_pre).or.num==0)cycle
            times = times+1
            tmpnew(times,1:2) = nei1(itmp,:)
            tmpnew(times,3) = num
        end do
        ! record the point  and point to the pre-site
        k = sum(tmpnew(:,3))                     !   2  for  two border  ,1   for  one border   ,    0  for  no  neighbor  border
            if(k = =2)then
                do itmp = 1,2 
                    i1 = tmpnew(itmp,1)
                    j1 = tmpnew(itmp,2)
                    call neinum(num,i1,j1)
                    if(num = =1)exit
                end do  !i
                call link(i1,j1)                                     !  link to this site
                call  drawcircle(2,i1,j1)                       !  drawcircle
                call link(i,j)                                         !   link to the previous site
                if(itmp<3)then
                i1 = tmpnew(3-itmp,1)
                j1 = tmpnew(3-itmp,2)
                else
                    i1 = nei1(4,1)
                    j1 = nei1(4,2)
                endif
            elseif(k = =1)then
                do itmp = 1,2
                    if(tmpnew(itmp,3) = =1)exit
                end do  !itmp
                    i1 = tmpnew(itmp,1)
                    j1 = tmpnew(itmp,2)
            else      !  consider the second neighbor site
                    call neighbor2(nei2,i,j)
                    do icycle = 1,6
                        i1 = nei2(icycle,1)
                        j1 = nei2(icycle,2)
                        if(i1 = =i_pre.and.j1==j_pre)cycle
                        call stateborder(num,i,j)
                        if(state(i1,j1) = =0.and.num==1)exit
                    end do   ! icycle
            endif
        i_pre = i
        j_pre = j
        i = i1
        j = j1
        if(i = =si.and.j==sj)exit
        read*
        write(10,*)i,j
    end do  !
    close(10)
    end subroutine boundary
!*************************************************
    !  link the data to the pointer
    subroutine link(i,j) 
    integer  i,j
    allocate(p)
        p.i = i
        p.j = j
        p.next = >null()
        if(associated(last))then
            last.next = >p
            last = >p
        else
            head = >p
            last = >p
        endif
    end subroutine link
!*************************************************
    !  judge a site border or not  ,  0 for no and 1 for yes
    subroutine stateborder(status,i,j)
    integer status,i,j
    integer junum
    call neighbornum(junum,i,j)
    if(state(i,j) = =0.and.junum>1.and.junum<6)then      !  the number of second neighbor sites is larger than 1 but less than 6
        status = 1
    else
        status = 0
    endif
    end subroutine stateborder
!*************************************************
    end module module_lattice
