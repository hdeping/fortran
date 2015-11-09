program main
    use graph_module
    real         :: t1 ,t2
    real         :: tv
    real         :: tvv
    integer      :: tvn
    integer      :: results
    open(10,file='mcenter.dat')
    open(11,file='test.dat')
    open(12,file='v_averger.dat')
    open(13,file='kaver.dat')
    open(14,file='v.dat')
call random_seed()
!    fric=4.0
!    goto 400
    call glutinitdisplaymode(ior(aux_double, aux_rgb))
    call glutinitposition(50, 10, int(graph_width), int(graph_height))
    results=fauxinitwindow("flow")
    call myinit()
    call initial_setup()
!    call saveconfig( 1 )
!    fspp=100.0d0

    !pause
    call glutidlefunc(loc(mydisplay))
    call glutreshapefunc(loc(myreshape))
    call keyregister()
    call glutmainloop(null)
    stop
400 print*,'aver_fric'
open(20,file='friction.dat')
print*,nsum/(lbox*lbox),'rou'
do tv = 0.0 , 5.0 , 0.1
    print*,tv
    call initial_setup()
    fric=tv
    tvv=0.0d0
    tvn=0
    do step = 1 , 2e4
        call viscek_move()
        if(step>1e4)then
            call cal_averv()
            tvv=tvv+averv
            tvn=tvn+1
        endif
        if(mod(step,1000)==0)print*,step
    enddo
    write(20,*)fric,tvv/tvn
enddo
end program main
