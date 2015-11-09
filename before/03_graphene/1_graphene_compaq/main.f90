program main
use common_module
use Graph_module
use lattice_module
use KMCevent_module
use event_module

real*8                ::    rate
integer                ::    i_level
real                ::    dEH1(6),dEH2(6)
real                ::    dET1(6),dET2(6)
integer,parameter    ::    nT_tot = 1E6
real*8                ::    freq(6,4)
integer                ::    theProc,icolor

call random_seed()
! 初始化绘图窗口及记录鼠标事件的变量
jUnit  =  GetActiveQQ()
jfbk   =  RegisterMouseEvent( jUnit, Mouse$LBUTTONDOWN, PauseProgram )
jfbk2  =  RegisterMouseEvent( jUnit, Mouse$RBUTTONDOWN, PauseProgram2 )

open(11,file = 'EqDis.dat')
open(12,file = 'flux.dat')
do i = 1,191
    c_ci(1) = 0.007+0.0001*(i-1)
    call initial()
    write(11,"(I5,6G18.9)") i,c_ci
    write(12,"(25G18.9)") c_ci(1),k_ciH,k_ciH_D,k_ciT,k_ciT_D
enddo !i
close(11)
close(12)
icolor = settextcolor(int2(0))
call InitGraphWindow(color_background)
pause
dET1 = dEiT_B-dEiH_B
dET2 = dEiT_B2-dEiH_B2    

write(fName,"('growthCnt.dat')")
open(13,file = fName)
do irun = 0,0 !20
    nt = 0
    c_ci(1) = 0.01    
    !ii = 5
    dEiT_B = dEiH_B+ 0.05*real(20-irun)*dET1
    dEiT_B2 = dEiH_B2+ 0.05*real(20-irun)*dET2
    cnt_growth = 0.0
    call initial()
    call SetGraphWindow( x_scr, y_scr, x_scr+w_scr, y_scr+h_scr )
    Call SetTextWindow( int2(46), int2(5), int2(47), int2(860) )
    call drawBonds()
    call drawSites()    
!    call drawSites_growthStatus()    
    open(10,file = 'count.dat')
    write(fName,"('cnt',I2.2,'.dat')") 20-irun
    open(12,file = fName)
!    do while(nt<nT_tot)
    do while(sum(cnt_growth)<100000)
        !更新各参数
        call update()
        ! 单击鼠标左键暂停
        if( jfg_pause  = = 1 ) then
            open(11,file = 'initial.dat')
            do ii = 1,nx_ltc
                write(11,"(60I5)") status(ii,:)
            enddo !ii
            close(11)
            pause
            jfg_pause  = 0
        else if(jfg_pause  = = 2) then
            call drawBonds()
            call drawSites_growthStatus()
            pause
            jfg_pause  = 0    
        endif
        !表面生长
    !    call growth()
        call growth2(theProc)
        write(12,"(I5)") theProc
    
        ! 移动观察窗口
        if(pos_front>nx_ltc-1) then
            call moveLattice()
            Call ClearScreen($GCLEARSCREEN)
            call drawBonds()
            call drawSites()
!            call drawSites_growthStatus()    
!            call redraw(pos_graphene,pos_front)
        endif

        ! 隔一段时间画一次图
    !    if(dt>1E-1) then
    !        dt = .0
        if(mod(nT,200) = =1) then
            call redraw(pos_graphene,pos_front+1)
!            call drawSites_growthStatus()    
            rate = dble(0.142*1.5*(pos_graphene+xPos_ltc-2))/dble(time)
            print*,nt, time,rate,int(sum(cnt_growth))
            open(11,file = 'dt.dat')
            write(11,"(7G18.9)") 0.05*real(20-irun),cnt_growth/sum(cnt_growth)
            close(11)
        endif
!        write(13,"(7G18.9)") 0.05*real(20-irun),sum(cnt_growth)
    enddo
    close(10)
    close(12)
    write(13,"(7G18.9)") 0.05*real(20-irun),cnt_growth/sum(cnt_growth)
!    call freqCNT(irun,freq)
enddo !irun
close(13)
end program main
