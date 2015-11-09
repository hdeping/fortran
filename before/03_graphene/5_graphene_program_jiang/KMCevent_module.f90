module KMCevent_module
use common_module
use lattice_module
!use Graph_module
real    :: p_add
contains

!******************************************************
subroutine growth()
real        :: tmp1
real        :: A0        ! 总趋势函数
real        :: drift(24)! 各过程趋势函数
integer        :: theCfg,theProc,flag,loc,iProc

drift(1:6) = k_ciH*n_cfgH
drift(7:12) = k_ciT*n_cfgT
drift(13:18) = k_ciH_D*n_cfgH_D
drift(19:24) = k_ciT_D*n_cfgT_D
A0 = sum(drift)
! 确定发生下一次反应的时间间隔
call random_number(tmp)
! 更新时间
time = time+log(1/tmp)/A0
dt = dt+log(1/tmp)/A0
nT = nT+1
! 确定当前时刻发生的过程
call random_number(tmp)
tmp = tmp*A0
do i = 1,24
    if(sum(drift(1:i))>tmp) then
        theProc = i
        exit
    endif
enddo !i
! 确定当前过程的具体构象并执行
call random_number(tmp1)
select case (theProc)
    case(:6)
        iProc = theProc
        theCfg = floor(tmp1*n_cfgH(iProc))+1
        do i = 1,iProc
            loc = cfgH(iProc,theCfg,2)+i-1
            if(loc>ny_ltc) loc = loc-ny_ltc
            status(cfgH(iProc,theCfg,1),loc) = 1
        enddo !i
    case(7:12)
        iProc = theProc-6
        theCfg = floor(tmp1*n_cfgT(iProc))+1
        do i = 1,iProc
            loc = cfgT(iProc,theCfg,2)+i-1
            if(loc>ny_ltc) loc = loc-ny_ltc
            status(cfgT(iProc,theCfg,1),loc) = 1
        enddo !i
    case(13:18)
        iProc = theProc-12
        theCfg = floor(tmp1*n_cfgH_D(iProc))+1
        do i = 1,iProc
            loc = cfgH_D(iProc,theCfg,2)+i-1
            if(loc>ny_ltc) loc = loc-ny_ltc
            status(cfgH_D(iProc,theCfg,1),loc) = 0
        enddo !i
    case(19:24)
        iProc = theProc-18
        theCfg = floor(tmp1*n_cfgT_D(iProc))+1
        do i = 1,iProc
            loc = cfgT_D(iProc,theCfg,2)+i-1
            if(loc>ny_ltc) loc = loc-ny_ltc
            status(cfgT_D(iProc,theCfg,1),loc) = 0
        enddo !i
    case default
endselect
end subroutine growth
!******************************************************
subroutine growth2(theProc)
integer                :: i_level
real                :: tmp1
real                :: A0        ! 总趋势函数
real,allocatable    :: drift(:)    ! 各过程趋势函数
real                :: drift1(24)
integer                :: np
integer                :: i_start
integer                :: theCfg,theProc,flag,loc,iProc
real                :: r_ciH(6),r_ciT(6),r_ciH_D(6),r_ciT_D(6)

r_ciH = k_ciH
r_ciT = k_ciT
r_ciH_D = k_ciH_D
r_ciT_D = k_ciT_D
i_level = 1
if(n_cfgH(1) = =0) then
    do i = 1,6
        if(r_ciT_D(i)>300.*r_ciT(i)) cycle !300.
!        if(n_cfgT(i)) cycle
        i_level = i
        exit
    enddo !i
endif
!if(nT-nT_stay>800)    then
!    i_level = 5
!    nT_stay = nT
!endif
if(sum(n_cfgH(i_level:6)+n_cfgT(i_level:6)) = =0) i_level=1
!i_level = 1
i_end = 7-i_level
np = i_end*4
allocate(drift(np))
drift(1:i_end) = r_ciH(i_level:)*n_cfgH(i_level:)
drift(i_end+1:2*i_end) = r_ciT(i_level:)*n_cfgT(i_level:)
drift(2*i_end+1:3*i_end) = r_ciH_D(i_level:)*n_cfgH_D(i_level:)
drift(3*i_end+1:4*i_end) = r_ciT_D(i_level:)*n_cfgT_D(i_level:)
A0 = sum(drift)
! 确定发生下一次反应的时间间隔
call random_number(tmp)
! 更新时间
write(11,*) nT,log(1./tmp)/A0/CE_time
time = time+log(1./tmp)/A0/CE_time
dt = dt+log(1./tmp)/A0/CE_time
nT = nT+1
! 确定当前时刻发生的过程
call random_number(tmp)
tmp = tmp*A0
do i = 1,np
    if(sum(drift(1:i))>tmp) then
        theProc = floor(real(i-1)/i_end)*6+mod(i-1,i_end)+i_level
        exit
    endif
enddo !i
! 确定当前过程的具体构象并执行
call random_number(tmp1)
i_tmp = 0;j_tmp=0;i_flag=0
select case (theProc)
    case(1:6)
        iProc = theProc
        theCfg = floor(tmp1*n_cfgH(iProc))+1
        do i = 1,iProc
            loc = cfgH(iProc,theCfg,2)+i-1
            if(loc>ny_ltc) loc = loc-ny_ltc
            status(cfgH(iProc,theCfg,1),loc) = 1
        enddo !i
    case(7:12)
        i_flag = 1
        iProc = theProc-6
        theCfg = floor(tmp1*n_cfgT(iProc))+1
        do i = 1,iProc
            loc = cfgT(iProc,theCfg,2)+i-1
            if(loc>ny_ltc) loc = loc-ny_ltc
            status(cfgT(iProc,theCfg,1),loc) = 1
            growthstatus(cfgT(iProc,theCfg,1),loc) = iProc
        enddo !i
    case(13:18)
        iProc = theProc-12
        theCfg = floor(tmp1*n_cfgH_D(iProc))+1
        do i = 1,iProc
            loc = cfgH_D(iProc,theCfg,2)+i-1
            if(loc>ny_ltc) loc = loc-ny_ltc
            status(cfgH_D(iProc,theCfg,1),loc) = 0
        enddo !i
    case(19:24)
        i_flag = 1
        iProc = theProc-18
        theCfg = floor(tmp1*n_cfgT_D(iProc))+1
        do i = 1,iProc
            loc = cfgT_D(iProc,theCfg,2)+i-1
            if(loc>ny_ltc) loc = loc-ny_ltc
            status(cfgT_D(iProc,theCfg,1),loc) = 0
            j_tmp = growthstatus(cfgT_D(iProc,theCfg,1),loc)
            if(j_tmp>i_tmp) i_tmp = j_tmp
            growthstatus(cfgT_D(iProc,theCfg,1),loc) = 0
        enddo !i
    case default
endselect
if(i_flag = =1) write(10,"(3I5)") i_level,iproc,theProc+100*i_tmp
deallocate(drift)
if(i_level>1)     nT_stay = nT
end subroutine growth2
!******************************************************
end module KMCevent_module
