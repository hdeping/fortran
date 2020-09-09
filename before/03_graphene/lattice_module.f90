module lattice_module
use common_module
!use Graph_module

contains
!*******************************************************************
subroutine update()
integer flag,loc,flag1,ii,value

do ii = pos_graphene,nx_ltc
    ! ����������graphene����λ��
    if(sum(status(ii,:)) = =ny_ltc.and.ii==pos_graphene+1) then
        pos_graphene = ii !
        nT_stay = nT
        do jj = 1,ny_ltc
            if(growthStatus(ii,jj)>0) cnt_growth(growthStatus(ii,jj)) = cnt_growth(growthStatus(ii,jj))+1.
        enddo !jj
    endif
    do jj = 1,ny_ltc                                                 ! ����x�������ռ��λ��
        if(status(ii,jj) = =1.and.ii>pos_front) pos_front=ii
     enddo                                                       !jj
enddo                                                           !ii    ! ���µ�ǰ�ɱ�ռ�ݵ�Ci�����б�,����������ͬһ�У��б������λ��Ϊ���͵ĵ�һ��Cλ�ã�Ҫ����ù����������겻�ܳ���������ռ�ݵ�����
n_cfgH = 0
n_cfgT = 0
do ii = pos_graphene,nx_ltc                                    !pos_graphene+4 !pos_front+1,�޸���2014-8-12,
    do jj = 1,ny_ltc                                                ! ������ѭ��Ϊ�����������������е�
        do j_c = 1,6                                                ! ����ѭ��Ϊ������(ii,jj)��Ϊ��������ci��
            flag = 0
            do jj_c = 1,j_c                                       ! ����ѭ��������ǰ���Ͱ���������λ��
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
                if(status(ii,loc) = =1) goto 100                  ! ����λ�ñ�ռ�ݣ���ǰci���Ͳ�����������i�����ci���Ͷ�������
                if(status(ii-1,loc)/ = 1) goto 100                ! ��鵱ǰ���Ƿ񳬹�������ռ�ݵ�����
                call getSiteType(ii,loc)                        ! ��¼�ù��Ͱ������ٸ�����λ
                if(mod(ii+loc,2) = =0) flag=flag+type_site
                if(type_site = =1.and.j_c==1) then
                    call getSiteClass(ii,jj)
                    if(n_neighbor = =2) flag=0
                endif
            enddo !jj_c

            loc = jj-1
            if(loc = =0) loc=ny_ltc 
            if(flag<1.and.status(ii,loc)/ = 1.and.mod(ii+jj,2)==1) goto 200        ! ��鹹���϶��Ƿ񳬹�������ռ�ݵ�����
            flag1 = status(ii,loc)                                    ! ��鹹���϶���λ�Ƿ񶼱�ռ��
            loc = jj+j_c
            if(loc>ny_ltc) loc = loc-ny_ltc
            if(flag<1.and.status(ii,loc)/ = 1.and.mod(ii+jj+j_c-1,2)==1) goto 200 ! ��鹹���¶��Ƿ񳬹�������ռ�ݵ�����
            flag1 = flag1+status(ii,loc)                                ! ��鹹���¶���λ�Ƿ񶼱�ռ��
            if(flag1 = =2) then 
                if(j_c<5) flag = -7+flag                                ! �����¶���λ����ռ�ݣ�����i<4�Ĺ��������Ƿ����������λ������������λ����
                                                                    ! �����ڿ�����������ʱ��������ù����Ƿ����������λ
            endif
            if(flag<1) then
                n_cfgH(j_c) =  n_cfgH(j_c)+1                            ! ������������λ�������¶���λ����ռ�ݣ���ù���Ϊ����������
                cfgH(j_c,n_cfgH(j_c),:) = (/ii,jj/)
            else
                n_cfgT(j_c) =  n_cfgT(j_c)+1                            ! ����ù���Ϊ����������
                cfgT(j_c,n_cfgT(j_c),:) = (/ii,jj/)
                if(j_c>1.and.flag1 = =2) then
                    do itmp = j_c+1,6
                        n_cfgT(itmp) =  n_cfgT(itmp)+1                ! �������Σ�����ǰΪc5�����Ҹù���ǡ�ò�ȫ���У���c6Ҳ����䵽�ù���
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
do ii = pos_graphene,nx_ltc-1           ! ԭΪpos_front+1,�޸���2014-8-12
    do jj = 1,ny_ltc                    ! ������ѭ��Ϊ�����������������е�
        do j_c = 1,6                    ! ����ѭ��Ϊ������(ii,jj)��Ϊ��������ci��
            flag = 0
            do jj_c = 1,j_c            ! ����ѭ��������ǰ���Ͱ���������λ��
                loc = jj+jj_c-1
                if(loc>ny_ltc) loc = loc-ny_ltc
                if(status(ii,loc) = =0) goto 300                        ! ����λ��δռ�ݣ���ǰci���Ͳ�����������i�����ci���Ͷ�������
                if(ii<nx_ltc.and.status(ii+1,loc) = =1) goto 300        ! ��鵱ǰ����������Ƿ���ռ�ݣ���ռ���򲻿��Ѹ�
                call getSiteType(ii,loc)                            ! ��¼�ù��Ͱ������ٸ�����λ
                if(mod(ii+loc,2) = =0) flag=flag+type_site
            enddo !jj_c
            loc = jj-1
            if(loc = =0) loc=ny_ltc 
            if((j_c = =1.or.flag<1).and.status(ii,loc)==1.and.mod(ii+jj,2)==0) goto 300        ! ��鹹���϶��Ƿ񳬹�������ռ�ݵ�����
            flag1 = status(ii,loc)                                    ! ��鹹���϶���λ�Ƿ�ռ��
            loc = jj+j_c
            if(loc>ny_ltc) loc = loc-ny_ltc
            if((j_c = =1.or.flag<1).and.status(ii,loc)==1.and.mod(ii+jj+j_c-1,2)==0) goto 400 ! ��鹹���¶��Ƿ񳬹�������ռ�ݵ�����
            flag1 = flag1*status(ii,loc)                                ! ��鹹���¶���λ�Ƿ񶼱�ռ��
            if(flag1 = =1) goto 400                                    ! �����¶���λ����ռ�ݣ��򲻿��Ѹ�
            if(flag<1) then
                n_cfgH_D(j_c) =  n_cfgH_D(j_c)+1                        ! ������������λ�������¶���λ����ռ�ݣ���ù���Ϊ����������
                cfgH_D(j_c,n_cfgH_D(j_c),:) = (/ii,jj/)
            else
                n_cfgT_D(j_c) =  n_cfgT_D(j_c)+1                        ! ����ù���Ϊ����������
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
! ��鵱ǰ��ķ��࣬class_site = 3: body; 2: surface; 1: active; 0: inactive
integer i,j,ss(3)
! neighbor 1
if(j = =1)then                                        !�˴���Ϊ2013-12-3 18��59
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
if(mod(i+j,2) = =0) then                              !�˴���Ϊ2013-12-3 18��59
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
! ��鵱ǰ�㷢���������̵����ͣ�1: atop; 0: hollow
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
    Neighbors(1,:) = (/i,ny_ltc/)                                  !�޸���2014-8-14
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
if(mod(i+j,2) = =0) then                                          !�˴���2013-12-3 18��52�޸�
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
! ����ǰ�۲촰����x����ƽ��ż����
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