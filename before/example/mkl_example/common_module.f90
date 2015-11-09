module common_module

real*8, PARAMETER                ::    PI=3.14159265
!-------------------------------------------------------------------------------------------------------
! controlling flag
integer                          ::    f_MSSA        ! ��־������1��Aλ����ȷ����0��δȷ��
integer                          ::    f_newinitial=1
integer                          ::    f_restart=0
integer                          ::    f_shake=1
integer                          ::    f_FDi_try=1
integer                          ::    f_fixend=1
integer                          ::    f_pullingN=1
integer                          ::    f_pullingC=1
integer                          ::    f_printV=0
!-------------------------------------------------------------------------------------------------------
! parameters and variables for Go model 
INTEGER, PARAMETER               ::    Num=67                    ! number of alpha C in protein I27
integer,parameter                ::    nx=3

real*8                           ::    x(num,nx)                ! positions of each bead (alpha C), 
real*8                           ::    v(num,nx)                ! velocities of each bead (alpha C), 
real*8                           ::    F_tot(num,3)
real*8                           ::    mass(num)                ! masses of each bead (alpha C), 
real*8                           ::    X0(num,nx)                ! native-state positions of each bead

real*8                           ::    A0(nx)                    ! end-to-end (N->1) vector in native state 
real*8, DIMENSION(NUM-1)         ::    R_bond,R0_bond            ! bond lengths and native-state bond lengths
real*8, DIMENSION(NUM-1,3)       ::    vct_bond                ! unit bond vectors
real*8, PARAMETER                ::    EPS_bond=378.0            ! ������ǿ�� 

real*8, DIMENSION(NUM-2)         ::    Theta,Theta0            ! angles and native-state angles of two adjacent bonds 
real*8, PARAMETER                ::    EPS_Theta=75.6            ! �Ƕ�Ťת����ǿ�� 

real*8, DIMENSION(NUM-3,4)       ::    DiTheta0                ! native-state dihedral angles of two next-nearest bonds
real*8, DIMENSION(NUM-3)         ::    DiTheta01                ! dihedral angles of two next-nearest bonds
real*8, DIMENSION(NUM-3)         ::    DiTheta                    ! dihedral angles of two next-nearest bonds
real*8, DIMENSION(NUM-3,4)       ::    EPS_DiTheta                ! �����Ťת����ǿ�� 

real*8, DIMENSION(NUM,NUM)       ::    Sigma_nonNei                ! repulsive radius or native state distance of particles i and j for abs(j-i)>2
real*8, DIMENSION(NUM,NUM)       ::    R0_nonNei                ! repulsive radius or native state distance of particles i and j for abs(j-i)>2
real*8, DIMENSION(NUM,NUM)       ::    R_nonNei                ! distance of particles i and j for abs(j-i)>2
real*8, DIMENSION(NUM,NUM)       ::    EPS_nonNei                ! interaction strength of particles i and j for abs(j-i)>2
INTEGER,DIMENSION(NUM,NUM)       ::    Matrix                    ! contact matrix, for abs(j-i)>2, bond= 0; nonbond=1; native contact=2 

integer                          ::    n_NC                    ! number of native contact of native state
real*8, parameter                ::    CoEff_break=1.5            ! coefficient of break distance (related to native distance) of a native contact 

CHARACTER                        ::    AF*22                    !�����γɹ淶PDB��ʽ���ַ�
!-------------------------------------------------------------------------------------------------------
! parameters and variables for dynamics
integer,parameter                ::    nt_rlx=1000
integer                          ::    nt_jump=10
real*8                           ::    op            ! һά�����
integer,parameter                ::    op_drct=1    ! ���������1����A��B��������ӣ�-1����С
real*8                           ::    MSSA        ! location of A
real*8                           ::    MSSB        ! location of B
real*8                           ::    temp=300.                ! temperature
real*8,parameter                 ::    kB=1./1000.                    
real*8,parameter                 ::    gamma=5.                ! friction constant
real*8, PARAMETER                ::    dt=1.E-2                ! time step
real*8                           ::    F=0.5                    ! pulling force (on bead 1)
integer,parameter                ::    nF=10                    
real*8                           ::    F_array(nF)= (/    1.428,1.,2.,3.,4.,&
                                                    5.,6.,7.,8.,9./)                ! pulling force array (on bead 1)
!---------------------------------------------------------------------------------------------------------
! ����ֲ�������ر���
integer,parameter     ::    np_r=50    
real*8                ::    dis_r(np_r,3)
real*8                ::    p_r=0.9                    ! ����ǰ�ۻ��ֲ���ֵ 
! ����ѧ������ز���
integer,parameter     ::    nT=3E4                    ! ��һ�����棺���߳���
integer               ::    n_CV                    ! �������棺��������
integer,allocatable   ::    len_cv(:)
integer,parameter     ::    niT=100                ! �������棺������󳤶ȣ�Ҫ������ϵ��ԥʱ���൱
integer,parameter     ::    nT_try=20*niT            ! ������nT_try����δ�ص�A̬Ҳδ������һ�����Ϊ�����ֹ���
integer               ::    cnt_CV(3)                ! ���ֹ��޵Ĺ��߼���,1:����һ���棬 2: ����A�� 3: ������    
! ������صĽṹ�塢ָ�롢�������Լ�����
integer,parameter     ::    n_cfg=100
real*8                ::    config(n_cfg,num,nx)
real*8                ::    myconfig(n_cfg,num,nx)
real*8                ::    V_config(n_cfg,num,nx)
real*8                ::    V_myconfig(n_cfg,num,nx)
real*8                ::    F_config(n_cfg,num,nx)
real*8                ::    F_myconfig(n_cfg,num,nx)
!-------------------------------------------------------------------------
type    ::    IF_link                                ! ����ṹ������IF_link
    integer                      :: i_IF          ! ��������
    real*8                       :: MSSA          ! MSSA���������С
    real*8                       :: r             ! ��MSSA�ľ���
    real*8                       :: cp(2)         ! cp ǰһ���浽��ǰ����öε��������ʣ��������������Ϊ�� 
    type(IF_link),pointer        :: prev
    type(IF_link),pointer        :: next
end type IF_link
type(IF_link),pointer            :: head_IF        ! ����ṹ���ͷָ��
type(IF_link),pointer            :: p_IF        ! ָʾ����ṹ�嵱ǰ��Ԫ��ָ��
!-------------------------------------------------------------------------
! other variables
character*50           ::    fName
integer                ::    ERR

contains
!---------------------------------------------------------------------------------------------------------
subroutine Save2xyz( fUnit )
    implicit none
    integer        fUnit, i

    write( fUnit, * ) Num
    write( fUnit, * )
        
    do i = 1, Num
        if(i==1) then
            write( fUnit, "(A8,3G12.6)") 'S    ',x(I,:)
        elseif(i==Num) then
            write( fUnit, "(A8,3G12.6)") 'S    ',x(I,:)
        else
            write( fUnit, "(A8,3G12.6)") 'C    ',x(I,:)
        endif
    end do
    return
end subroutine Save2xyz
subroutine  new()
!******************************************************************************************************************
end module common_module
