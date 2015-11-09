        Module module_graph
        use module_common
        use module_coor
                                     
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
        subroutine setgraphwindow(ux,uy,dx,dy)
        integer ux,uy,dx,dy
        call setviewport(int2(ux),int2(uy),int2(dx),int2(dy))
        end subroutine setgraphwindow
!***********************************************************************
     subroutine  trigonal()
     real*8  x,y
     do i = 1,n
     do j = 1,m
    call coor_tri(x,y,i,j)
    call drawcircle(1,x,y,r1)
     end do !j
     end do !i
     end  subroutine trigonal
!***********************************************************************
    subroutine hexagonal(x0,y0)
    real*8 x0,y0,theta
    theta = 10
    do i = 1,int(n*2*r1/(r2*sqrt(3.)))
     do j = 1,int(m*r1/(r2*sqrt(3.)))
     call coor_hexa(x,y,i,j,x0,y0)
    !call rot_coor(x,y,x,y,theta,x0,y0)
     call drawcircle(2,x,y,r2)
     end do !j
     end do !i
    end subroutine hexagonal
!***********************************************************************
    subroutine drawcircle(jcolor,x,y,r)
    integer  jcolor
   real*8  x,y,r
   real*8 x1,y1,x2,y2
   x1 = x-r
   x2 = x+r
   y1 = y-r
   y2 = y+r
   if(jcolor = =1)    then
    icolor = setcolorrgb(mycolor(10))
    icolor = Ellipse_w($GFILLINTERIOR,x1,y1,x2,y2)
    else
    icolor = setcolorrgb(mycolor(4))
    icolor = Ellipse_w($GFILLINTERIOR,x1,y1,x2,y2)
    endif
    end subroutine drawcircle
!***********************************************************************
    End Module module_graph