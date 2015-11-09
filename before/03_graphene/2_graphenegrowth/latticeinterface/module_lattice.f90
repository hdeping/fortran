module module_lattice 
    use module_common
    
    
    contains
!*************************************************
    subroutine neighbornum(num,i,j)
    integer  num1,num2,i,j,times1,times2,i1,i2,j2
    times1 = 0
    times2 = 0
    call  neighbor(nei,i,j)
    do i1 = 1,3
    i2 = nei(i1,1)
    j2 = nei(i1,2)
    if(i2>0.and.i2<n+1.and.j2>0.and.j2<m+1)then
    times1 = times1+1
    if(state(i2,j2) = =0)times2=times2+1
    end if
    end do !i1
    num1 =  times1
    num2 =  times2
    end subroutine neighbornum
!*************************************************
    ! nearest 
    subroutine neighbor(nei,i,j)
    integer  nei(3,2),i,j,k
    nei(1,:) = (/i+1,j/)
    nei(2,:) = (/i-1,j/)
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
    nei(3,:) = (/i,j+1/)    
    case default
    nei(3,:) = (/i,j-1/)  
    end select
    end subroutine neighbor
!*************************************************
    end module module_lattice

    
    
    
    
!    subroutine growthborder()
!     integer  i1,i2
!     open(10,file = "output.txt")
!    do i = 1,n
!    do j = 1,m
!    !call neighbornum(num1,num2,i,j)
!    !if(num2<3)then
!    !state(i,j) = 1
!    !else
!    !call neighbornum2(i1,i2,i,j)
!    !if(i2<6)state(i,j) = 1
!    !endif
!    call neighbornum2(i1,i2,i,j)
!    if(i2<6)state(i,j) = 1
!    !print *,i,j,i2
!    !read*
!    write(10,*)i,j,i1,i2
!    end do !j
!    end do !i
!    end subroutine growthborder
!!*************************************************
!    subroutine latticeborder()
!    integer  i1,i2
!    do i = 1,n
!    do j = 1,m
!    call neighbornum(num1,num2,i,j)
!    if(num1<3)then
!    state(i,j) = 1
!    else
!    call neighbornum2(i1,i2,i,j)
!    if(i1<6)state(i,j) = 1
!    !print *,i1
!    !read*
!    endif
!    end do !j
!    end do !i
!    end subroutine latticeborder
!!*************************************************
!    ! the second nearest number , num1 for the border, num2 for the growth
!    subroutine neighbornum2(num1,num2,i,j)
!     integer  num1,num2,i,j,times1,times2,i1,i2,j2
!    times1 = 0
!    times2 = 0
!    call  neighbor2(nei2,i,j)
!    do i1 = 1,6
!    i2 = nei2(i1,1)
!    j2 = nei2(i1,2)
!    if(i2>0.and.i2<n+1.and.j2>0.and.j2<m+1)then
!    times1 = times1+1
!    if(state(i2,j2) = =0)times2=times2+1
!    endif
!    end do !i1
!    num1 =  times1
!    num2 =  times2
!    end subroutine neighbornum2
!*************************************************
    ! nearest number, num1 for the border, num2 for the growth
    
    
    
    
!!*************************************************
!    ! the second nearest
!     subroutine neighbor2(nei2,i,j)
!    integer  nei2(6,2),i,j
!    nei2(1,:) = (/i+2,j/)
!    nei2(2,:) = (/i-2,j/)
!    nei2(3,:) = (/i+1,j-1/)
!    nei2(4,:) = (/i+1,j+1/)
!    nei2(5,:) = (/i-1,j+1/)
!    nei2(6,:) = (/i-1,j-1/)
!    end subroutine neighbor2
