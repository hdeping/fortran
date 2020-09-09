        Module module_graph
        use module_common
                                     
    Contains
!***********************************************************************
    ! /* ------ Set the Screen Resolution and Color ------------------------ */ !
        Subroutine InitGraphWindow( jColor )
            Integer jColor
            integer (2)     result

            Logical StatusMode
            Integer oldColor, retInt

            theFrame.type  =  QWIN$MAX    !SET
            theFrame.x     =  100
            theFrame.y     =  50
            theFrame.h     =  650      ! This Number is in the Unit of Pixels
            theFrame.w     =  650
            retInt  =  SetWSizeQQ( QWIN$FRAMEWINDOW, theFrame )

            theChild.type  =  QWIN$MAX
             theChild.x  =  70          
             theChild.y  =  50
             theChild.h  =  100          ! Note: This Number is in the Unit of Character (Height,Width)
             theChild.w  =  100
            retInt  =  SetWSizeQQ( 0, theChild )
      
            thisScreen.numxPixels   =  -1 
            thisScreen.numyPixels   =  -1
            thisScreen.numTextCols  =  -1
            thisScreen.numTextRows  =  -1
            thisScreen.numColors    =  -1
            thisScreen.fontSize     =  -1

            StatusMode  =  SetWindowConfig( thisScreen )
            if( .NOT. StatusMode ) StatusMode  =  SetWindowConfig( thisScreen )

            g_jColor  =  jColor    ! Store the BackGround Color

            ! 设置子窗口的背景颜色
            result  =  SETBKCOLORRGB (jcolor)
            
            Call ClearScreen( $GCLEARSCREEN ) 

            Return
        End Subroutine InitGraphWindow   
!************************************************************************
    subroutine atomnum(numco,numo2,state)
    integer  numco,numo2,state(n,m)
    numco = 0
     numo2 = 0
     do i = 1,n
     do j = 1,m
     if(state(i,j) = =1)numco=numco+1
     if(state(i,j) = =2)numo2=numo2+1
     end do !j
     end do !i
     end subroutine atomnum
!***********************************************************************   
       subroutine linew(x1,y1,x2,y2)
        real *8     x1,y1,x2,y2
        call moveto_w(x1,y1,xy)
        jcolor = setcolorrgb(mycolor(10))
        jcolor = lineto_w(x2,y2)
        end subroutine linew
 !************************************************************************        
    subroutine latticeinterface()
         do i = 1,n
         do j = 1,m
         call drawcircle(state(i,j),i,j)
         end do !j
         end do !i
         
    end subroutine latticeinterface
!************************************************************************
    subroutine drawlines()
    do j = 1,m
        do i = 1,n-1
            call line(i,j,i+1,j)
        end do  !i
    end do   !  j
    ! 偶数列
    do j = 2,m,2
        do i = 1,n
            if(mod(i,2) = =1)then
                call line(i,j,i,j-1)
            else
                if(j<m)call line(i,j,i,j+1)
            endif
        end do   !  i
    end do   !  j
    
    end subroutine drawlines
!************************************************************************
    !  draw a circle in the site (i,j)
      subroutine drawcircle(status,i,j)
      ! status  0 for blank ,1 for occupied
      integer  i,j,status
      real*8  x,y,x1,x2,y1,y2
     call getcoor(x,y,i,j)
      x1 = x-r
      x2 = x+r
      y1 = y-r
      y2 = y+r
          
        if(status = =0)then
    jcolor = setcolorrgb(mycolor(1))
     jcolor = Ellipse_w($GBORDER,x1,y1,x2,y2)
    !elseif(status = =1)then
    !jcolor = setcolorrgb(mycolor(1))
    ! jcolor = Ellipse_w($GFILLINTERIOR,x1,y1,x2,y2)
    elseif(status = =1)then
     jcolor = setcolorrgb(mycolor(1))
     jcolor = Ellipse_w($GFILLINTERIOR,x1,y1,x2,y2)
     else
      jcolor = setcolorrgb(mycolor(4))
     jcolor = Ellipse_w($GBORDER,x1,y1,x2,y2)
    endif
      end subroutine drawcircle
!************************************************************************
      subroutine getcoor(x,y,i,j)
      real*8  x,y,x0,y0
      integer  i,j,s,s1,k
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
    !k = mod(j,4)
    s = (j-1)/2
    s1 = (i-1)/2
    x0 = 0
    y0 = 0
    select case(k)
   case(1)
   x = x0+s*3*r0+r
   y = y0+s1*sqrt(3.)*r0+r+r0/2.*sqrt(3.)
   case(2)
   x = x0+s*3*r0+r+r0/2.
   y = y0+s1*sqrt(3.)*r0+r
   case (3)
    x = x0+s*3*r0+3/2.*r0+r
   y = y0+s1*sqrt(3.)*r0+r
   case default
   x = x0+s*3*r0+2*r0+r
   y = y0+s1*sqrt(3.)*r0+r+r0/2.*sqrt(3.)
   end select
      end subroutine getcoor
!************************************************************************
    subroutine rects(status)
    integer  status(n,m)
    do i = 1,n
    do j = 1,m
    call rect(status(i,j),j,i)
    end do !j
    end do !i
    
    end subroutine rects
!***********************************************************************
    subroutine rect(status,i,j)
    integer   status,i,j
    real(8)    x1,x2,y1,y2
    
    x1 = (i-1)*k
    y1 = (j-1)*k
    x2 = i*k
    y2 = j*k
    if(status = =0)then
    jcolor = setcolorrgb(mycolor(17))
     jcolor = rectangle_w($GFILLINTERIOR,x1,y1,x2,y2)
    elseif(status = =1)then
    jcolor = setcolorrgb(mycolor(1))
     jcolor = rectangle_w($GFILLINTERIOR,x1,y1,x2,y2)
    else
     jcolor = setcolorrgb(mycolor(4))
     jcolor = rectangle_w($GFILLINTERIOR,x1,y1,x2,y2)
    endif
    
    end subroutine rect
!******************************************************
    !  draw a line between site(i1,j1)  and  site(i2,j2)
    subroutine  line(i1,j1,i2,j2)
    integer i1,j1,i2,j2
    real*8   x1,y1,x2,y2
    call getcoor(x1,y1,i1,j1)
    call getcoor(x2,y2,i2,j2)
    call moveto_w(x1,y1,xy)
    jcolor = setcolorrgb(mycolor(1))
    jcolor = lineto_w(x2,y2)
    end subroutine   line
!******************************************************
subroutine  rgbtohls(h,l,s,r,g,b)
    integer  r,g,b,itmp,jtmp    !  range 0 to 255
    integer h,l,s        ! hux,lightness , saturation ,range 0 to 240
   integer  a(6),array(3),i1,j1,judge,num
   do itmp = 1,6
   a(itmp) = 40*(itmp-1)
   end do !itmp
    itmp = min(r,g,b)
    jtmp = max(r,g,b)
    if(itmp/ = jtmp)then
    array = (/r,g,b/)
    do i1 = 1,3
    if(itmp = =array(i1))exit
    end do !i1
    do j1 = 1,3
    if(jtmp = =array(j1))exit
    end do !i1
    if(j1 = =1.and.i1==3)judge=1
    if(j1 = =2.and.i1==3)judge=2
    if(j1 = =2.and.i1==2)judge=3
    if(j1 = =3.and.i1==2)judge=4
    if(j1 = =3.and.i1==1)judge=5
    if(j1 = =1.and.i1==1)judge=6
    !  value for h
    select case(judge)
    case(6,2,4)
    number = jtmp-array(6-i1-j1)
    case default
    number = array(6-i1-j1)-itmp
    end select
    h = a(judge)+40*number/(jtmp-itmp)
    end if
    ! value for saturation
    if(itmp+jtmp<256)then
    s = 240*(jtmp-itmp)/(itmp+jtmp)
    else
    s = 240*(jtmp-itmp)/(510-itmp-jtmp)
    endif
    itmp = itmp-(itmp+8)/17
    jtmp = jtmp-(jtmp+8)/17
    l = (itmp+jtmp)/2            !  value for lightness
    
    end subroutine rgbtohls
!***************************************************************
 subroutine  hlstorgb(r,g,b,h,l,s)
    
    
    end subroutine hlstorgb
!***************************************************************
    End Module module_graph