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
!***********************************************************************   
       subroutine line(x1,y1,x2,y2)
        integer     x1,y1,x2,y2
        call moveto(x1,y1,intxy)
        jcolor = setcolorrgb(mycolor(10))
        jcolor = lineto(int2(x2),int2(y2))
        end subroutine line
!***********************************************************************   
       subroutine linew(x1,y1,x2,y2)
        real *8     x1,y1,x2,y2
        call moveto_w(x1,y1,xy)
        jcolor = setcolorrgb(mycolor(10))
        jcolor = lineto_w(x2,y2)
        end subroutine linew
 !************************************************************************        
    subroutine latticeinterface(n,m,state)
         integer   n,m
         integer   state(n,m)
         
         do i = 1,n
         do j = 1,m
         call drawcircle(state(i,j),i,j)
         end do !j
         end do !i
         
      end subroutine latticeinterface
!************************************************************************
      subroutine drawcircle(status,i,j)
      ! status  0 for blank ,1 for occupied
      integer  i,j,k,status,s,s1
      real*8  x,y,x1,x2,y1,y2,x0,y0
      
      !open(10,file = "data.txt")
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
      x1 = x-r
      x2 = x+r
      y1 = y-r
      y2 = y+r

        if(status = =0)then
     jcolor = setcolorrgb(mycolor(1))
     jcolor = Ellipse_w($GBORDER,x1,y1,x2,y2)
        elseif(status = =1)then
     jcolor = setcolorrgb(mycolor(1))
     jcolor = Ellipse_w($GFILLINTERIOR,x1,y1,x2,y2)
        else               !  for the border of graphene
      jcolor = setcolorrgb(mycolor(4))
     jcolor = Ellipse_w($GFILLINTERIOR,x1,y1,x2,y2)
        endif
    end subroutine drawcircle
!***************************************************************
    End Module module_graph