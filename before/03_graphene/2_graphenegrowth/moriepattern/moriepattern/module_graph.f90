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
        
!***********************************************************************  
   ! draw n lines ,length  =  600
   subroutine lines(n)
   integer  n,k
   k = n/2+1
   x1 = x0-length/2
   x2 = x0+length/2
   y1 = y0-k*height
   do i = 1,n
   y1 = y1+height
   y2 = y1
   call line(x1,y1,x2,y2)
   end do  !i
   end subroutine lines
!***********************************************************************  
        subroutine lines_rot(n,theta)
        integer  n,k
        real*8  theta,angle,x,y,length_new
        angle = theta*pi/180+pi/2. 
        length_new = length
         k = n/2
         do i = 1,n
         x = x0+(i-k-1)*height*cos(angle)
         y = y0+(i-k-1)*height*sin(angle)
          x1 = x-cos(angle+pi/2.)*length_new/2.
         x2 = x+cos(angle+pi/2.)*length_new/2.
         y1 = y-sin(angle+pi/2.)*length_new/2.
         y2 = y+sin(angle+pi/2.)*length_new/2.
        
         call line_new(x1,y1,x2,y2)
         end do  !i
        end subroutine lines_rot
!***********************************************************************  
        !!  theta  ::  Deg pattern(0~180)
        !subroutine rot_coor(xnew,ynew,theta,x,y)
        !real*8  xnew,ynew,x,y,theta,angle,s1,s2
        !angle = theta*pi/180.
        !!print *,x,y,cos(angle),sin(angle)
        !!print *,x-x0,y-y0
        !s1 = (x-x0)*cos(angle)-(y-y0)*sin(angle)
        !s2 = (x-x0)*sin(angle)+(y-y0)*cos(angle)
        !xnew = x0+s1
        !ynew = y0+s2
        !!print '(4f10.4)',x0,y0,xnew,ynew
        !!print '(2f10.4)',s1,s2
        !end subroutine rot_coor
!***********************************************************************
        subroutine line_new(x1,y1,x2,y2)
       real*8     x1,y1,x2,y2
        call moveto_w(x1,y1,xy)
        jcolor = setcolorrgb(mycolor(10))
        jcolor = lineto_w(x2,y2)
        end subroutine line_new
!***********************************************************************   
       subroutine line(x1,y1,x2,y2)
       real*8     x1,y1,x2,y2
        call moveto_w(x1,y1,xy)
        jcolor = setcolorrgb(mycolor(1))
        jcolor = lineto_w(x2,y2)
        end subroutine line
!***********************************************************************
        
        
!***********************************************************************
        
    End Module module_graph