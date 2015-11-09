module dynamic_module
use common_module
use math_module
use Go_dynamics
use Go_Structure

contains
!***************************************************************************************************************
subroutine initial(jF)
implicit none
integer    :: jF,it,i_IF,in,i,j
type(IF_link),pointer    ::    p1

call initial_GO()
call CreatLink(head_IF) ! 启用新界面串行结构
head_IF%i_IF=0
head_IF%i_IF=0
p_IF=>head_IF

call getF_tot(F_tot)
DO I=1, NUM
!    if(f_fixEnd==1.and.i==num) cycle
    do j=1,3
        F_tot(i,j)=F_tot(i,j)+SQRT(2.0*TEMP*gamma*kB/dt)*GaussianNoise(dble(0.),dble(1.))
    enddo !j
    IF (f_pullingN==1.and.I==1) F_tot(i,:)=F_tot(i,:)+F*A0    
    IF (f_pullingC==1.and.I==Num) F_tot(i,:)=F_tot(i,:)-F*A0    
ENDDO

if(f_restart==1) then
    print*,'输入重新开始的起始界面:'
    read(*,*) i_IF
    write(fName,"('cfg/cfg_',I3.3,'.dat')") i_IF
    open(11,file=fName)
    in=nx*Num
    do j=1,n_cfg
        read(11,"(<in>G18.9)") config(j,:,1),config(j,:,2),config(j,:,3)
        read(11,"(<in>G18.9)") V_config(j,:,1),V_config(j,:,2),V_config(j,:,3)
        read(11,"(<in>G18.9)") F_config(j,:,1),F_config(j,:,2),F_config(j,:,3)
    enddo !j
    close(11)
    write(fName,"('IF_',I2.2,'.dat')") jF
    open(10,file=fName)
    call InputLink(10,i_IF,p_IF)
    close(10)
    f_MSSA=1
    MSSA=p_IF%MSSA
    print*,'MSSA',MSSA
else
    f_MSSA=0
    if(f_newinitial==1) then
        write(fName,"('initial/initialX.dat')")
        OPEN(100+jF, FILE=fName)
        write(fName,"('initial/initialV.dat')")
        OPEN(200+jF, FILE=fName)
        write(fName,"('initial/initialF.dat')")
        OPEN(300+jF, FILE=fName)
        do it=1,Num
            read(100+jF,'(3G18.9)') x(it,:)    
            read(200+jF,'(3G18.9)') V(it,:)
            read(300+jF,'(3G18.9)') F_tot(it,:)
        enddo !it    
        close(100+jF)
        close(200+jF)
        close(300+jF)
        write(fName,"('ts_',I2.2,'.dat')") jF
        OPEN(10, FILE=fName)
        do it=1,nt_rlx
            call runDynamic()
            call getOP(op)
            if(mod(it,100)==0) then
                write(10,"(3G18.9)") dble(it*nt_jump)*dt,op
                print*,'    ',dble(it*dt*nt_jump),op
            endif
        enddo !it
        close(10)
        write(fName,"('initial/initialX_',I2.2,'.dat')") jF
        OPEN(100+jF, FILE=fName)
        write(fName,"('initial/initialV_',I2.2,'.dat')") jF
        OPEN(200+jF, FILE=fName)
        write(fName,"('initial/initialF_',I2.2,'.dat')") jF
        OPEN(300+jF, FILE=fName)
        do it=1,Num
            write(100+jF,'(3G18.9)') x(it,:)    
            write(200+jF,'(3G18.9)') V(it,:)
            write(300+jF,'(3G18.9)') F_tot(it,:)
        enddo !it    
        close(100+jF)
        close(200+jF)
        close(300+jF)
    endif
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
endif
endsubroutine initial
!----------------------------------------------------------------------------
subroutine getCurve_1st(ID,in,nT,y,jF)
implicit    none
integer        ::    ID,nT,i,in,j,jF
real*8        ::    y(nT),tmp
!open(110,file='op.dat')
do i=1,nT
    call RunDynamic()
    call getOP(tmp)
    do while(tmp>44.5)
        write(fName,"('initial/initialX_',I2.2,'.dat')") jF
        OPEN(400+jF, FILE=fName)
        write(fName,"('initial/initialV_',I2.2,'.dat')") jF
        OPEN(200+jF, FILE=fName)
        write(fName,"('initial/initialF_',I2.2,'.dat')") jF
        OPEN(300+jF, FILE=fName)
        do j=1,Num
            read(400+jF,'(3G18.9)') x(j,:)    
            read(200+jF,'(3G18.9)') V(j,:)
            read(300+jF,'(3G18.9)') F_tot(j,:)
        enddo !it    
        close(400+jF)
        close(200+jF)
        close(300+jF)
        call RunDynamic()
        call getOP(tmp)
    enddo
    y(i)=tmp
!    if(mod(i,100)==0) print*,dble(i*nt_jump)*dt,y(i)
!    write(110,*) dble(i*nt_jump)*dt,y(i)
    write(ID,"(<in>G18.9)") x(:,1),x(:,2),x(:,3)
    write(ID,"(<in>G18.9)") v(:,1),v(:,2),v(:,3)
    write(ID,"(<in>G18.9)") F_tot(:,1),F_tot(:,2),F_tot(:,3)
enddo !i
!close(110)

end subroutine getCurve_1st
!----------------------------------------------------------------------------
subroutine getCurve_ith(ID,in,y,p)
implicit    none
integer                ::    ID,in,i,ny,it,i_tmp,cnt
real*8                ::    y(nT),f(nt+1),tmp
type(IF_link),pointer    ::    p

myconfig=config
V_myconfig=V_config
F_myconfig=F_config
ny=0;n_CV=0;f=0;cnt=0
do while(ny<nT)
    n_cv=n_cv+1
    call random_number(tmp)        ! 每条轨线的初始值都从前一个界面的构型随机产生
    i_tmp=floor(dble(n_cfg)*tmp)+1
    if(i_tmp<1) i_tmp=1
    if(i_tmp>n_cfg) i_tmp=n_cfg
    x=myconfig(i_tmp,:,:)
    V=V_myconfig(i_tmp,:,:)
    F_tot=F_myconfig(i_tmp,:,:)
    cnt=cnt+1
    f(n_cv)=cnt
    write(ID,"(I9,<in>G18.9)") n_cv,x(:,1),x(:,2),x(:,3)
    write(ID,"(I9,<in>G18.9)") n_cv,v(:,1),v(:,2),v(:,3)
    write(ID,"(I9,<in>G18.9)") n_cv,F_tot(:,1),F_tot(:,2),F_tot(:,3)
    call getOP(op)
    if((.not.InMSSA(op)).and.ny<nT) then
        ny=ny+1
        y(ny)=op
    endif
    do it=2,niT
        call RunDynamic()
        call getOP(op)
        cnt=cnt+1
        write(ID,"(I9,<in>G18.9)") n_cv,x(:,1),x(:,2),x(:,3)
        write(ID,"(I9,<in>G18.9)") n_cv,v(:,1),v(:,2),v(:,3)
        write(ID,"(I9,<in>G18.9)") n_cv,F_tot(:,1),F_tot(:,2),F_tot(:,3)
        if(InMSSA(op)==1) exit    ! 如果轨线返回第0界面内，则轨线结束
        if(InIF(op,p)) cycle        ! 如果返回当前界面内侧，不计入统计 
        if(ny<nT) then
            ny=ny+1 ! 否则记录该点
            if(mod(ny,1000)==0) print*,'n_cv=',n_cv,'ny=',ny
            y(ny)=op
        endif 
    enddo !it
enddo !iCurve
f(n_cv+1)=cnt+1
close(10)
if(allocated(len_CV)) deallocate(len_cv)
allocate(len_cv(n_cv))
open(10,file='test_len.dat')
do i=1,n_cv
    len_cv(i)=f(i+1)-f(i)
    write(10,*) i,len_cv(i)
enddo !i
close(10)
!pause
return
end subroutine getCurve_ith
!**********************************************************************************************
end module dynamic_module
