!! This Module Defines Some Common-Used Customizer Graphic Modules
Module Graph_Module
Use MSFLIB
use common_module
use lattice_module

! Definitions of Globally Used Variables for Graphic Settings -------
  Integer g_xLeft, g_yLeft, g_xRight, g_yRight    ! g --> Global
  Integer g_jColor
  Integer g_xLeft1, g_yLeft1, g_xRight1, g_yRight1    ! g --> Global

  Type (WindowConfig) thisScreen
  Type (QWInfo)       theFrame, theChild 

! Define my Own Color Palette for 16 Colors --------------

  Integer :: myColor(16) = (/ #000000,   &        ! Black
                              #FFFFFF,   &        ! Bright White
                                 #000080,   &        ! DULL Red  
                              #0000FF,   &        ! Bright Red
                              #008000,   &        ! Dull Green
                              #00FF00,   &        ! Bright Green
                              #008080,   &        ! Dull Yellow
                              #00FFFF,   &        ! Bright Yellow
                              #800000,   &        ! Dull Blue
                              #FF0000,   &        ! Bright Blue
                              #800080,   &        ! Dull Magenta 
                              #FF00FF,   &        ! Bright Magenta
                              #808000,   &        ! Dull Turquoise
                              #FFFF00,   &        ! Bright Turquoise 
                              #808080,   &        ! Dark Gray
                              #C0C0C0   /)        ! Bright Gray
 
Contains

! Set the Screen Resolution and Color -------------------------------
Subroutine InitGraphWindow( jColor )
USE MSFLIB
Integer jColor                

Logical StatusMode
Integer oldColor, retInt

theFrame.type = QWIN$SET
theFrame.x    = 100
theFrame.y    = 50
theFrame.h    = 860      ! This Number is in the Unit of Pixels
theFrame.w    = 1400
retInt = SetWSizeQQ( QWIN$FRAMEWINDOW, theFrame )

theChild.type = QWIN$MAX
theChild.x = 700          
theChild.y = 50
theChild.h = 1000          ! Note: This Number is in the Unit of Character (Height,Width)
theChild.w = 1000
retInt = SetWSizeQQ( 0, theChild )
  
thisScreen.numxPixels  = -1 
thisScreen.numyPixels  = -1
thisScreen.numTextCols = -1
thisScreen.numTextRows = -1
thisScreen.numColors   = -1
thisScreen.fontSize    = -1

StatusMode = SetWindowConfig( thisScreen )
if( .NOT. StatusMode ) StatusMode = SetWindowConfig( thisScreen )

g_jColor = jColor    ! Store the BackGround Color

oldColor = SetBKColorRGB( mycolor(jColor) )
Call ClearScreen( $GCLEARSCREEN ) 
Return
End Subroutine InitGraphWindow       
!******************************************************************
! Set the Window Region for the Graphics ----------------------------
Subroutine SetGraphWindow( x1, y1, x2, y2 )
Use MSFLIB
Integer(2)    x1, y1, x2, y2

g_xLeft  = x1 ; g_yLeft  = y1 
g_xRight = x2 ; g_yRight = y2                                

Call SetViewPort( x1, y1, x2, y2 )
Call ClearScreen( $GVIEWPORT ) 

Return     
End Subroutine SetGraphWindow
!******************************************************************!
Subroutine Draw2DLattic_Hexagon(ii, jj, jColor,flag) 

Integer  jColor,ii,jj
integer     flag                ! 0为空心，其他为实心
Integer(2)  x1,x2,y1,y2,jRet,inte

x1 = g_xLeft ; x2 = g_xRight ; y1 = g_yLeft ; y2 = g_yRight
!! Draw the site  -------------------
if(mod(ii+jj,2)==0) then
    x1=g_xLeft+(ii-1)*len_grid+3
else
    x1=g_xLeft+(ii-1)*len_grid+int2(len_grid/2)
endif
y1 = hh_grid*jj
!inte = int(hh_grid/2)
inte=5
call getSiteType(ii,jj)
if(growthstatus(ii,jj)>0.and.type_site==1) then
    jRet = SetColorRGB( mycolor(Color_spc))
    jRet = Ellipse( $GFILLINTERIOR, x1-inte, y1-inte, x1+inte, y1+inte ) 
    inte=3
endif
jRet = SetColorRGB( mycolor(jColor))
if(flag==0) then
    jRet = Ellipse( $GBORDER, x1-inte, y1-inte, x1+inte, y1+inte ) 
else
    jRet = Ellipse( $GFILLINTERIOR, x1-inte, y1-inte, x1+inte, y1+inte ) 
endif
End Subroutine Draw2DLattic_Hexagon
!******************************************************************!
Subroutine Draw2DLattic_Hexagon_line(i,j,ii,jj,jColor ) 
!USE MSFLIB

Integer  jColor,i,j,ii,jj
Integer  x1,y1,jRet,inte,dummy
integer(2)  x2,y2
Type  (XYCOORD ) xyPos1

if(j==1.and.jj==ny_ltc) return
if(jj==1.and.j==ny_ltc) return
jRet = SetColorRGB( mycolor(jColor))
!! Draw the line  -------------------===
if(mod(i+j,2)==0) then
    x1=g_xLeft+(i-1)*len_grid+3
else
    x1=g_xLeft+(i-1)*len_grid+int2(len_grid/2)
endif
y1 = hh_grid*j
if(mod(ii+jj,2)==0) then
    x2=g_xLeft+(ii-1)*len_grid+3
else
    x2=g_xLeft+(ii-1)*len_grid+int2(len_grid/2)
endif
y2 = hh_grid*jj
!inte = int(hh_grid/2)
inte=4
Call MoveTo( x1, y1, xyPos1 )
dummy = LineTo( x2, y2 )

End Subroutine Draw2DLattic_Hexagon_line
!******************************************************************!
subroutine DrawSites()    
! 画出当前构型(需要表征当前构型的三个数组status，surface，active)
do i=1,nx_ltc
    do j=1,ny_ltc
        call getSiteClass(i,j)
        selectcase(Class_Site)
            case(3)
                if(showBody==1) call Draw2DLattic_Hexagon(i,j,color_body,1)
            case(2)
                if(showSurface==1) call Draw2DLattic_Hexagon(i,j,color_surface,1)    
            case(1)
                if(showActive==1) call Draw2DLattic_Hexagon(i,j,color_active,0)    
            case(0)
                if(showInactive==1)    call Draw2DLattic_Hexagon(i,j,color_inactive,0)
            case default
        endselect
!        call getSiteType(i,j)
!        if(type_site==1.and.showSPC==1) then
!            if(status(i,j)==0) then
!                call Draw2DLattic_Hexagon(i,j,color_spc,0)
!            else
!                call Draw2DLattic_Hexagon(i,j,color_spc,0)
!            endif
!        endif
    enddo !j
enddo !i 
end subroutine DrawSites
!******************************************************************!
subroutine DrawSites_growthStatus()    
integer        ::    icolor
! 画出当前构型(需要表征当前构型的三个数组status，surface，active)
do i=1,nx_ltc
    do j=1,ny_ltc
        call getSiteClass(i,j)
        selectcase(Class_Site)
            case(0)
                if(showInactive==1)    call Draw2DLattic_Hexagon(i,j,color_inactive,0)
            case default
                selectcase(growthstatus(i,j))
                    case(1) 
                        icolor=1
                    case(2)
                        icolor=6
                    case(3)
                        icolor=4
                    case(4)
                        icolor=8
                    case(5)
                        icolor=10
                    case(6)
                        icolor=12
                    case default
                        icolor=color_body
                endselect
                call Draw2DLattic_Hexagon(i,j,icolor,1)
        endselect
!        call getSiteType(i,j)
!        if(type_site==1) call Draw2DLattic_Hexagon(i,j,color_spc,0)
    enddo !j
enddo !i 
end subroutine DrawSites_growthStatus
!******************************************************************!
subroutine DrawBonds()
!画出当前CC键
do i=1,nx_ltc
    do j=1,ny_ltc
        if(status(i,j)==0) cycle
        call getNeighbors(i,j)
        do ii=2,3
            if(status(Neighbors(ii,1),Neighbors(ii,2))==0) cycle
            if(i==1.and.ii==3) cycle
            if(i==ny_ltc.and.ii==2) cycle
            call Draw2DLattic_Hexagon_line(i,j,Neighbors(ii,1),Neighbors(ii,2),color_bond )
        enddo !ii
    enddo !j
enddo !i 

end subroutine DrawBonds
!******************************************************************!
subroutine redraw(i_s,i_e)
integer i_s,i_e,nei(2),i_s1

i_s1=i_s-6
if(i_s1<=0) i_s1=1
do i=i_s1,i_e
    do j=1,ny_ltc
        call getNeighbors(i,j)
        do ii=1,3
            nei=Neighbors(ii,:)
            if(status(i,j)==1.and.status(nei(1),nei(2))==1) then
                call Draw2DLattic_Hexagon_line(i,j,nei(1),nei(2),color_bond) 
            else    
                call Draw2DLattic_Hexagon_line(i,j,nei(1),nei(2),color_background)
            endif 
            call getSiteClass(nei(1),nei(2))
            call Draw2DLattic_Hexagon(nei(1),nei(2),color_background,1)
            selectcase(Class_Site)
                case(3)
                    if(showBody==1) call Draw2DLattic_Hexagon(nei(1),nei(2),color_body,1)
                case(2)
                    if(showSurface==1) call Draw2DLattic_Hexagon(nei(1),nei(2),color_surface,1)    
                case(1)
                    if(showActive==1) call Draw2DLattic_Hexagon(nei(1),nei(2),color_active,0)    
                case(0)
!                    if(showInactive==1) call Draw2DLattic_Hexagon(nei(1),nei(2),color_inactive,0)
                    call Draw2DLattic_Hexagon(nei(1),nei(2),color_inactive,0)
                case default
            endselect
            call getSiteType(nei(1),nei(2))
            if(type_site==1) call Draw2DLattic_Hexagon(nei(1),nei(2),color_spc,0)
        enddo !ii
        call Draw2DLattic_Hexagon(i,j,color_background,1)
        call getSiteClass(i,j)
        selectcase(Class_Site)
            case(3)
                if(showBody==1) call Draw2DLattic_Hexagon(i,j,color_body,1)
            case(2)
                if(showSurface==1) call Draw2DLattic_Hexagon(i,j,color_surface,1)    
            case(1)
                if(showActive==1) call Draw2DLattic_Hexagon(i,j,color_active,0)    
            case(0)
!                if(showInactive==1) call Draw2DLattic_Hexagon(i,j,color_inactive,0)
                call Draw2DLattic_Hexagon(i,j,color_inactive,0)
            case default
        endselect
        call getSiteType(i,j)
        if(type_site==1) call Draw2DLattic_Hexagon(i,j,color_spc,0)
    enddo !j            
enddo !i
do i=i_s1,i_e+1
    do j=1,ny_ltc
        call getSiteType(i,j)
        if(type_site==2) call Draw2DLattic_Hexagon(i,j,color_spc,0)
    enddo !j            
enddo !i
end subroutine  redraw
!******************************************************************!
End Module Graph_Module
