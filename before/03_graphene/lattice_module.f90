module lattice_module
use common_module
!use Graph_module

contains
!*******************************************************************
subroutine update()
integer flag,loc,flag1,ii,value

do ii = pos_graphene,nx_ltc
    ! 更新完整的graphene界面位置
    if(sum(status(ii,:)) = =ny_ltc.and.ii==pos_graphene+1) then
        pos_graphene = ii !
        nT_stay = nT
        do jj = 1,ny_ltc
            if(growthStatus(ii,jj)>0) cnt_growth(growthStatus(ii,jj)) = cnt_growth(growthStatus(ii,jj))+1.
        enddo !jj
    endif
    do jj = 1,ny_ltc                                                 ! 更新x方向最大占据位置
        if(status(ii,jj) = =1.and.ii>pos_front) pos_front=ii
     enddo                                                       !jj
enddo                                                           !ii    ! 更新当前可被占据的Ci构型列表,构型限制在同一列，列表给出的位置为构型的第一个C位置，要求可用构型纵向坐标不能超过临列已占据的坐标
n_cfgH = 0
n_cfgT = 0
do ii = pos_graphene,nx_ltc                                    !pos_graphene+4 !pos_front+1,修改于2014-8-12,
    do jj = 1,ny_ltc                                                ! 这两重循环为遍历生长区域内所有点
        do j_c = 1,6                                                ! 这重循环为遍历以(ii,jj)点为起点的所有ci簇
            flag = 0
            do jj_c = 1,j_c                                       ! 这重循环遍历当前构型包含的所有位点
                loc = jj+jj_c-1
                if(loc>ny_ltc)then
                    loc = mod(loc,ny_ltc)                
                    if(loc = =0)loc=ny_ltc
                endif
                value = ii
                if(value>ny_ltc)then
                    value = mod(value,ny_ltc)
                    if(value = =0)value=ny_ltc
                endif
                if(status(ii,loc) = =1) goto 100                  ! 若该位置被占据，则当前ci构型不可用且所有i更大的ci构型都不可用
                if(status(ii-1,loc)/ = 1) goto 100                ! 检查当前点是否超过临列已占据的坐标
                call getSiteType(ii,loc)                        ! 记录该构型包含多少高能量位
                if(mod(ii+loc,2) = =0) flag=flag+type_site
                if(type_site = =1.and.j_c==1) then
                    call getSiteClass(ii,jj)
                    if(n_neighbor = =2) flag=0
                endif
            enddo !jj_c

            loc = jj-1
            if(loc = =0) loc=ny_ltc 
            if(flag<1.and.status(ii,loc)/ = 1.and.mod(ii+jj,2)==1) goto 200        ! 检查构型上端是否超过临列已占据的坐标
            flag1 = status(ii,loc)                                    ! 检查构型上端临位是否都被占据
            loc = jj+j_c
            if(loc>ny_ltc) loc = loc-ny_ltc
            if(flag<1.and.status(ii,loc)/ = 1.and.mod(ii+jj+j_c-1,2)==1) goto 200 ! 检查构型下端是否超过临列已占据的坐标
            flag1 = flag1+status(ii,loc)                                ! 检查构型下端临位是否都被占据
            if(flag1 = =2) then 
                if(j_c<5) flag = -7+flag                                ! 若上下端临位都被占据，对于i<4的构型无论是否包含高能量位，都按低能量位处理
                                                                    ! 否则在考虑能量问题时，仍需检查该构型是否包含高能量位
            endif
            if(flag<1) then
                n_cfgH(j_c) =  n_cfgH(j_c)+1                            ! 不包含高能量位或者上下端临位都被占据，则该构型为低能量构型
                cfgH(j_c,n_cfgH(j_c),:) = (/ii,jj/)
            else
                n_cfgT(j_c) =  n_cfgT(j_c)+1                            ! 否则该构型为高能量构型
                cfgT(j_c,n_cfgT(j_c),:) = (/ii,jj/)
                if(j_c>1.and.flag1 = =2) then
                    do itmp = j_c+1,6
                        n_cfgT(itmp) =  n_cfgT(itmp)+1                ! 特殊情形，若当前为c5构型且该构型恰好补全改列，则c6也可填充到该构型
                        cfgT(itmp,n_cfgT(itmp),:) = (/ii,jj/)
                    enddo !itmp
                endif
            endif
200            cycle                    
        enddo !j_c
100        cycle    
    enddo !jj
enddo !ii
n_cfgH_D = 0
n_cfgT_D = 0
do ii = pos_graphene,nx_ltc-1           ! 原为pos_front+1,修改于2014-8-12
    do jj = 1,ny_ltc                    ! 这两重循环为遍历生长区域内所有点
        do j_c = 1,6                    ! 这重循环为遍历以(ii,jj)点为起点的所有ci簇
            flag = 0
            do jj_c = 1,j_c            ! 这重循环遍历当前构型包含的所有位点
                loc = jj+jj_c-1
                if(loc>ny_ltc) loc = loc-ny_ltc
                if(status(ii,loc) = =0) goto 300                        ! 若给位置未占据，则当前ci构型不可用且所有i更大的ci构型都不可用
                if(ii<nx_ltc.and.status(ii+1,loc) = =1) goto 300        ! 检查当前点的右临列是否已占据，若占据则不可脱附
                call getSiteType(ii,loc)                            ! 记录该构型包含多少高能量位
                if(mod(ii+loc,2) = =0) flag=flag+type_site
            enddo !jj_c
            loc = jj-1
            if(loc = =0) loc=ny_ltc 
            if((j_c = =1.or.flag<1).and.status(ii,loc)==1.and.mod(ii+jj,2)==0) goto 300        ! 检查构型上端是否超过临列已占据的坐标
            flag1 = status(ii,loc)                                    ! 检查构型上端临位是否被占据
            loc = jj+j_c
            if(loc>ny_ltc) loc = loc-ny_ltc
            if((j_c = =1.or.flag<1).and.status(ii,loc)==1.and.mod(ii+jj+j_c-1,2)==0) goto 400 ! 检查构型下端是否超过临列已占据的坐标
            flag1 = flag1*status(ii,loc)                                ! 检查构型下端临位是否都被占据
            if(flag1 = =1) goto 400                                    ! 若上下端临位都被占据，则不可脱附
            if(flag<1) then
                n_cfgH_D(j_c) =  n_cfgH_D(j_c)+1                        ! 不包含高能量位或者上下端临位都被占据，则该构型为低能量构型
                cfgH_D(j_c,n_cfgH_D(j_c),:) = (/ii,jj/)
            else
                n_cfgT_D(j_c) =  n_cfgT_D(j_c)+1                        ! 否则该构型为高能量构型
                cfgT_D(j_c,n_cfgT_D(j_c),:) = (/ii,jj/)
            endif
400            cycle                    
        enddo !j_c
300        cycle    
    enddo !jj
enddo !ii
end subroutine update
!*******************************************************************
subroutine getSiteClass(i,j)
! 检查当前点的分类，class_site = 3: body; 2: surface; 1: active; 0: inactive
integer i,j,ss(3)
! neighbor 1
if(j = =1)then                                        !此处改为2013-12-3 18：59
    ss(1) = status(i,ny_ltc)
else
ss(1) = status(i,j-1)
endif
! neighbor 2
if(j = =ny_ltc)then
    ss(2) = status(i,1)
else
ss(2) = status(i,j+1)
endif
! neighbor 3
if(mod(i+j,2) = =0) then                              !此处改为2013-12-3 18：59
    if(i = =1)then
        ss(3) = 1
    else
    ss(3) = status(i-1,j)
    endif 
else
    if(i = =nx_ltc)then
        ss(3) = 0
    else
    ss(3) = status(i+1,j)
    endif
endif
class_site = sum(ss)
n_neighbor = sum(ss)
if(status(i,j) = =1.and.class_site<3) class_site=2
if(status(i,j) = =0.and.class_site>0) class_site=1
end subroutine getSiteClass
!*******************************************************************
subroutine getSiteType(i,j)
! 检查当前点发生生长过程的类型，1: atop; 0: hollow
integer :: i,j
integer :: i_pos,j_pos,i1

type_site = 0
i_pos = mod(i+xPos_ltc-2,nx_barrier)
if(i_pos = =0) i_pos=nx_barrier
j_pos = mod(j,ny_barrier)
if(j_pos = =0) j_pos=ny_barrier
select case (i_pos)
    case(2,4,15,17)
        if(j_pos = =4.or.j_pos==6) type_site=1
    case(3,16)
        if(j_pos = =3.or.j_pos==5.or.j_pos==7) type_site=1
    case(5,7,12,14)
        if(j_pos = =14.or.j_pos==16) type_site=1
    case(6,13)
        if(j_pos = =13.or.j_pos==15.or.j_pos==17) type_site=1
    case default
end select

end subroutine getSiteType
!*******************************************************************
subroutine getNeighbors(i,j)
! get neighbors of site(i,j)
integer    i,j
! neighbor 1
if(j = =1) then
    Neighbors(1,:) = (/i,ny_ltc/)                                  !修改于2014-8-14
else
Neighbors(1,:) = (/i,j-1/)
endif
! neighbor 2
if(j = =ny_ltc) then
    Neighbors(2,:) = (/i,1/)
else
Neighbors(2,:) = (/i,j+1/)
endif
! neighbor 3
if(mod(i+j,2) = =0) then                                          !此处于2013-12-3 18：52修改
    if(i = =1)then
     Neighbors(3,:) = (/nx_ltc,j/)   
    else
    Neighbors(3,:) = (/i-1,j/)
    endif
else
    if(i = =nx_ltc)then
    Neighbors(3,:) = (/nx_ltc,j/)    
    else
    Neighbors(3,:) = (/i+1,j/)
    endif
endif
endsubroutine getNeighbors
!*******************************************************************
subroutine moveLattice()
! 将当前观察窗口沿x正向平移偶数列
integer    dj,j_status(nx_ltc,ny_ltc)
dj = pos_graphene-1
if(mod(dj,2) = =1) dj=dj-1
if(dj<1) then
    print*,"Warning: the present configure is out of the observation window!"
    stop
endif
xPos_ltc = xPos_ltc+dj
pos_front = pos_front-dj
pos_graphene = pos_graphene-dj
j_status = status
status(1:pos_front,:) = j_status(dj+1:dj+pos_front,:)
status(1+pos_front:nx_ltc,:) = 0
j_status = growthstatus
growthstatus(1:pos_front,:) = j_status(dj+1:dj+pos_front,:)
growthstatus(1+pos_front:nx_ltc,:) = 0
!cfgH(:,:,1) = cfgH(:,:,1)-dj
!cfgT(:,:,1) = cfgT(:,:,1)-dj
end subroutine moveLattice
!*******************************************************************
end module lattice_module