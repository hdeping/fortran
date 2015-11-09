module Interface_module
use common_module
use math_module
use datalink_module
use dynamic_module

contains
!**********************************************************************************************
! �����·�Ӧ��������ӳ���,����:
!                newInterface: �����½���
!**********************************************************************************************
subroutine newInterface(jF,p_prev)
! determine L0
implicit    none
integer                ::    jF
real*8                ::    op_array(nT),ri,dOP
! ������ر���
type(IF_link),pointer::    p_prev
! ��������
integer                ::    j,iN

write(*,"('��ǰ����:',I5)") p_prev%i_IF+1 

write(*,"('        ���ڲ������ߡ�')") ! �����������
write(fName,"('curve_',I3.3,'.dat')") jF
open(10,file=fName)
in=Num*nx
if(f_MSSA==0) then
    call getCurve_1st(10,in,nT,op_array,jF)
    MSSA=sum(op_array)/dble(nT)
else
    call getCurve_ith(10,in,op_array,p_prev)
endif
close(10)

op_array=dble(op_drct)*(op_array-MSSA)
call getDis_R(nT,op_array,np_r,dis_r)

open(10,file='test.dat')
do j=1,np_r
    write(10,"(3G18.9)") dis_r(j,:)
enddo !i
close(10)

dOP=dis_r(2,1)-dis_r(1,1)
do j=1,np_r
    if(dis_r(j,3)>p_r) then
        if(j>1) then
            ri=dOP*((p_r-dis_r(j-1,3))/(dis_r(j,3)-dis_r(j-1,3))+0.5)+dis_r(j-1,1)
        else
            ri=dis_r(j,1)-dOP*(1.-p_r/dis_r(1,3))+0.5*dOP    
        endif
        exit
    endif
enddo !j

write(*,"('        ��¼������Ϣ��')")
call NewLinkUnit(p_IF)
p_IF%i_IF=p_prev%i_IF+1
if(f_MSSA==0) then
    p_IF%MSSA=MSSA
else
    p_IF%MSSA=p_prev%MSSA
endif
p_IF%r=ri
print*,"    r=",p_IF%r

write(*,"('        ��¼���湹�͡�')")
! �ҵ������ϵĹ��Ͳ�������������������
write(fName,"('curve_',I3.3,'.dat')") jF
open(10,file=fName)
in=Num*nx
if(f_MSSA==0) then
    call FindConfig_1st(10,in,nt,p_IF)
    f_MSSA=1
else
    call FindConfig_ith(10,in,len_cv,p_IF)    
endif
close(10)
write(fName,"('cfg/cfg_jF',I2.2,'_iF',I3.3,'.dat')") jF,p_IF%i_IF
open(11,file=fName)
in=3*Num
do j=1,n_cfg
    write(11,"(<in>G18.9)") config(j,:,1),config(j,:,2),config(j,:,3)
    write(11,"(<in>G18.9)") V_config(j,:,1),V_config(j,:,2),V_config(j,:,3)
    write(11,"(<in>G18.9)") F_config(j,:,1),F_config(j,:,2),F_config(j,:,3)
enddo !j
close(11)
write(*,"('��ǰ����:',I5,'��ɡ�')") p_prev%i_IF+1 
end subroutine newInterface
!**********************************************************************************************
! Ѱ���½��湹������ӳ���,����:
!                FindConfig_1st: ��һ������Ĺ������ļ���
!                FindConfig_ith: ��������Ĺ����������ʵļ���
!**********************************************************************************************
subroutine FindConfig_1st(ID,in,nt,p)
integer        ::    ID,in,nt
type(IF_link),pointer    ::    p,p1
integer        ::    it,f

open(601,file='cp/phi.dat')
! ǰһ���������ÿһ�����߱�Ȼ���ﵱǰ������߻ص�A����
cnt_CV=0;it=1
read(ID,"(<in>G18.9)") x(:,1),x(:,2),x(:,3)
read(ID,"(<in>G18.9)") v(:,1),v(:,2),v(:,3)
read(ID,"(<in>G18.9)") F_tot(:,1),F_tot(:,2),F_tot(:,3)
call getOP(op)
do while(cnt_cv(1)<n_cfg)
    f=InMSSA(op)
    it=it+1
    if(.not.EOF(ID)) then
        read(ID,"(<in>G18.9)") x(:,1),x(:,2),x(:,3)
        read(ID,"(<in>G18.9)") v(:,1),v(:,2),v(:,3)
        read(ID,"(<in>G18.9)") F_tot(:,1),F_tot(:,2),F_tot(:,3)
    else
        call RunDynamic()
    endif
    call getOP(op)
    do while(op>44.5)
        write(fName,"('initial/initialX_',I2.2,'.dat')") jF
        OPEN(100+jF, FILE=fName)
        write(fName,"('initial/initialV_',I2.2,'.dat')") jF
        OPEN(200+jF, FILE=fName)
        write(fName,"('initial/initialF_',I2.2,'.dat')") jF
        OPEN(300+jF, FILE=fName)
        do it=1,Num
            read(100+jF,'(3G18.9)') x(it,:)    
            read(200+jF,'(3G18.9)') V(it,:)
            read(300+jF,'(3G18.9)') F_tot(it,:)
        enddo !it    
        close(100+jF)
        close(200+jF)
        close(300+jF)
        call RunDynamic()
        call getOP(op)
    enddo
    if((.not.InMSSA(op)).and.f==1) then    
        cnt_cv(1)=cnt_cv(1)+1
        p%cp(1)=dble(cnt_cv(1))/dble(it*nt_jump*dt)
        p%cp(2)=dble(cnt_cv(3))/dble(it*nt_jump*dt)
        write(601,"(I5.5,2G18.9)") cnt_cv(1),dble(it*nt_jump*dt),dble(cnt_cv(1))/dble(it*nt_jump*dt)
        config(cnt_cv(1),:,:)=x
        V_config(cnt_cv(1),:,:)=v
        F_config(cnt_cv(1),:,:)=F_tot
    endif
enddo !while
close(601)
close(111)
end subroutine FindConfig_1st
!----------------------------------------------------------------------------
subroutine FindConfig_ith(ID,in,Ni_cv,p)
implicit    none
integer        ::    ID,in,Ni_cv(n_CV)
type(IF_link),pointer    ::    p
integer        ::    it,ii,i,i_c,flag,f
real*8        ::    x_tmp

write(fName,"('cp/cp_',I3.3,'.dat')") p%i_IF
open(601,file=fName)
! ��ǰ���������ÿһ�����߱�Ȼ������һ������߻ص�A����
cnt_cv=0
do i=1,n_CV
    it=0;flag=0;f=1
    do while(it<=Ni_cv(i).or.flag==0)
        it=it+1
        if(it>Ni_cv(i)) then
            call RunDynamic()
        else
            read(ID,"(I9,<in>G18.9)") i_c,x(:,1),x(:,2),x(:,3)
            read(ID,"(I9,<in>G18.9)") i_c,v(:,1),v(:,2),v(:,3)
            read(ID,"(I9,<in>G18.9)") i_c,F_tot(:,1),F_tot(:,2),F_tot(:,3)
        endif
        if(flag==1) cycle
        if(it>nT_try) then        ! ����3��δ��AҲδ����һ����
            cnt_cv(3)=cnt_cv(3)+1
            flag=1                ! �����ǰ����û�ж�ȡ�꣬���ȡ������
            cycle
        endif
        call getOP(op)
        if(InMSSA(op))  then    ! ����2���ص�A����
            cnt_cv(2)=cnt_cv(2)+1
            flag=1                
            cycle
        endif
        ! �ɹ���Խ���¼����
        if((.not.InIF(op,p)).and.f==1)  then    ! ����1��������һ����
            cnt_cv(1)=cnt_cv(1)+1
            flag=1                
            p%cp(1)=dble(cnt_cv(1))/dble(sum(cnt_cv))
            p%cp(2)=dble(cnt_cv(3))/dble(sum(cnt_cv))
            write(601,*) sum(cnt_cv),p%cp,sum(p%cp)
            if(cnt_cv(1)<=n_cfg) then
                config(cnt_cv(1),:,:)=x
                V_config(cnt_cv(1),:,:)=V
                F_config(cnt_cv(1),:,:)=F_tot
            endif
            if(cnt_cv(1)>=n_cfg) goto 100
            cycle
        endif
        f=InIF(op,p)
    enddo
enddo !i

do while(cnt_cv(1)<n_cfg)
    call random_number(x_tmp)        
    i=floor(dble(n_cfg)*x_tmp)+1
    if(i<1) i=1
    if(i>n_cfg) i=n_cfg
    x=myconfig(i,:,:)
    V=V_myconfig(i,:,:)
    F_tot=F_myconfig(i,:,:)
    it=0;f=1
    do while(1)
        it=it+1
        call RunDynamic()    
        if(it>nT_try) then        ! ����3��δ��AҲδ����һ����
            cnt_cv(3)=cnt_cv(3)+1
            exit
        endif
        call getOP(op)
        if(InMSSA(op))  then    ! ����2���ص�A����
            cnt_cv(2)=cnt_cv(2)+1
            exit
        endif
        ! �ɹ���Խ���¼����
        if((.not.InIF(op,p)).and.f==1)  then    ! ����1��������һ����
            cnt_cv(1)=cnt_cv(1)+1
            p%cp(1)=dble(cnt_cv(1))/dble(sum(cnt_cv))
            p%cp(2)=dble(cnt_cv(3))/dble(sum(cnt_cv))
            write(601,*) sum(cnt_cv),p%cp,dble(cnt_cv(3))/dble(sum(cnt_cv))
            if(cnt_cv(1)<=n_cfg) then
                config(cnt_cv(1),:,:)=x
                V_config(cnt_cv(1),:,:)=V
                F_config(cnt_cv(1),:,:)=F_tot
            endif
            exit
        endif
        f=InIF(op,p)
    enddo
enddo 
100 close(601)
end subroutine FindConfig_ith
!*********************************************************************************************
end module Interface_module
