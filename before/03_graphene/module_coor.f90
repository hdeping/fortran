module module_coor
    use module_common
     
    contains
!***********************************************************************
    subroutine gethexaarray(hexa)
    real*8 hexa(n,m,2),xtmp,ytmp
    integer   iitmp
    integer   xarray(3,2)
    x0 = 0
    y0 = 0
    do i = 1,n
    do j = 1,m
   call   coor_hexa(xtmp,ytmp,i,j,x0,y0)
   hexa(i,j,:) = (/xtmp,ytmp/)
    end do !j
    end do !i
    end subroutine gethexaarray
!***********************************************************************
    ! rotate anticlockwisely by the point (x1,y1)  (the center of the pattern)
    subroutine rot_coor(xnew,ynew,x,y,theta,x0,y0)
    real*8 xnew,ynew,x,y,theta,x0,y0
    real*8 angle,x1,y1,x2,y2,s1,s2
    angle = theta*pi/180
    call coor_hexa(x1,y1,1,1,x0,y0)
    call coor_hexa(x2,y2,n,m,x0,y0)
   print *,x1,y1,x2,y2
    x1 = (x1+x2)/2.
    y1 = (y1+y2)/2.
    y1 = 900.-y1
   print *,"theta  =  ",theta,"angle = ",angle
   print *,"the center"
   print *,x1,y1
    x2 = x
    y2 = 900.-y
   print *,"the initial coordinate"
   print *,x,y
   print *,"the transform coordinate"
   print *,x2,y2
    s1 = (x2-x1)*cos(theta)-(y2-y1)*sin(theta)
    s2 = (x2-x1)*sin(theta)+(y2-y1)*cos(theta)
    x2 = x1+s1
    y2 = y1+s2
   print *,"the rotation coordinate "
   print *,x2,y2
    xnew = x2
    ynew = 900.-y2
   print *,"the final coordinate"
   print *,xnew,ynew
    read*
    end subroutine rot_coor
!***********************************************************************
    subroutine coor_hexa(x,y,i,j,x0,y0)
    integer  i,j,k,s,s1
    integer  xarray(3,2)
    real*8 x,y,x0,y0
    
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
    s = (j-1)/2
    s1 = (i-1)/2
    select case(k)
   case(1)
   x = x0+s*6*r2+r2
   y = y0+s1*2*sqrt(3.)*r2+r2+r2*sqrt(3.)
   case(2)
   x = x0+s*6*r2+2*r2
   y = y0+s1*2*sqrt(3.)*r2+r2
   case (3)
    x = x0+s*6*r2+4*r2
   y = y0+s1*2*sqrt(3.)*r2+r2
   case default
   x = x0+s*6*r2+r2*5
   y = y0+s1*2*sqrt(3.)*r2+r2+r2*sqrt(3.)
   end select
    end subroutine coor_hexa
!***********************************************************************
    subroutine gettriarray(tri)
    real*8 tri(n-1,m-1,2),xtmp,ytmp,x,y
    integer   iitmp
    integer   xarray(3,2)
    do i = 1,n-1
    do j = 1,m-1
    xarray(:,1) = (/i,i+1,i/)
    xarray(:,2) = (/j,j,j+1/)
    xtmp = 0
    ytmp = 0
    do iitmp = 1,3
    call coor_tri(x,y,xarray(iitmp,1),xarray(iitmp,2))
    xtmp = xtmp+x
    ytmp = ytmp+y
    end do
    xtmp = xtmp/3.0
    ytmp = ytmp/3.0
    tri(i,j,:) = (/xtmp,ytmp/)
    end do !j
    end do !i
    end subroutine gettriarray
!***********************************************************************
    subroutine coor_tri(x,y,i,j)
    integer  i,j,k,s
    real*8  x,y,x0,y0
    x0 = 200
    y0 = 100
    k = mod(j-1,2)+1
    s = (j-1)/2
    if(k = =1)then
    x = x0+2*sqrt(3.)*r1*s+r1
    y = 2.*r1*(i-1)+r1+y0
    else
    x = x0+2*sqrt(3.)*r1*s+r1+sqrt(3.)*r1
    y = 2.*r1*(i-1)+2*r1+y0
    endif
    end subroutine coor_tri
!***********************************************************************
    end module module_coor