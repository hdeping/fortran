MODULE Go_Structure
use common_module
use math_module
CONTAINS 
!***********************************************************************************************************
subroutine getStructureVariable()

!!!!!!!!!!!!键长
call getBond()
!!!!!!!!!!!!键键夹角
call getAngle()
!!!!!!!!!!!!二面角
call getDiAngle()
!!!!!!!!!!!!自然连接和非键连接长度
call getNonNeighbor()
end subroutine getStructureVariable
!***********************************************************************************************************
subroutine initial_GO()
real*8    ::    r_tmp(num)
real*8,allocatable ::    RNat(:)

!!!!!!!!!!!!初始速度
V=0.
!!!!!!!!!!!!从原始文件中读取质量
CALL readMass()
!mass=mass/109.   ! 质量单位
mass=1.   ! 质量单位
!!!!!!!!!!!!从原始文件中读取个粒子坐标
CALL readStructure()                
X=X0
!!!!!!!!!!!!计算首尾相连的粒子的向量
A0=X0(1,:)-X0(NUM,:)
A0=A0/SQRT(A0(1)**2+A0(2)**2+A0(3)**2)
!!!!!!!!!!!!读取自然状态键长
call readBond()
!!!!!!!!!!!!读取自然状态键键夹角
call readAngle()
!!!!!!!!!!!!读取自然状态二面角
call readDiAngle()
!!!!!!!!!!!!读取自然状态自然连接和非键连接矩阵、自然连接长度和排斥半径、相互作用强度
call readNativeContact()
call readNonBond()
!if(f_FDi_try==1) then
call getDiAngle()
DiTheta01=DiTheta
!endif
call getStructureVariable()
allocate(RNat(n_nc))
R0_bond=R_bond
Theta0=Theta    
call getR_Native_Repulsion(n_nc,RNat,r_tmp)
where(matrix==1) Sigma_nonNei=R_nonNei*1.12246
end subroutine initial_GO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read the mass of partical from mass.txt
SUBROUTINE readMass()    
CHARACTER    ::        FR*18     

OPEN(301, FILE="GO_NativeState/mass.txt")
DO I=1, NUM
    READ(301, '(A18, F10.6)') FR, mass(i)
ENDDO 

CLOSE(301)
END SUBROUTINE readMass
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read the partical position of native state from ORG.txt
SUBROUTINE readStructure()    
CHARACTER    ::        FR*31     

OPEN(301, FILE="GO_NativeState/ORG.txt")
DO I=1, NUM
    READ(301, '(A31, F7.3, 2F8.3, A22 )') FR, X0(I,:), AF
ENDDO 
CLOSE(301)
END SUBROUTINE readStructure
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read the native bond lengths from bond.txt
SUBROUTINE readBond()    
CHARACTER    ::        FR*28     

OPEN(301, FILE="GO_NativeState/bond.txt")
DO I=1, NUM-1
    READ(301, '(A28, F10.6)') FR, R0_bond(i)
ENDDO 
CLOSE(301)
END SUBROUTINE readBond
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! get the bond lengths from particle positions
SUBROUTINE getBond()    

do i=1,num-1
    vct_bond(i,:)=x(i+1,:)-x(i,:)
    R_bond(i)=SQRT(vct_bond(i,1)**2+vct_bond(i,2)**2+vct_bond(i,3)**2)
    vct_bond(i,:)=vct_bond(i,:)/R_bond(i)
enddo !i
END SUBROUTINE getBond
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read the native angles from angle.txt
SUBROUTINE readAngle()    
CHARACTER    ::    FR*35     

OPEN(301, FILE="GO_NativeState/angle.txt")
DO I=1, NUM-2
    READ(301, '(A35, F12.6)') FR, Theta0(i)
ENDDO 
CLOSE(301)
Theta0=Theta0*PI/180.
END SUBROUTINE readAngle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! get the angles from particle positions
SUBROUTINE getAngle()    
real*8    ::    A(3), B(3)
real*8    ::    TT

DO I=1, NUM-2
    A=X(I,:)-X(I+1,:)
    B=X(I+2,:)-X(I+1,:)
    TT=sum(A*B)/(SQRT(A(1)**2+A(2)**2+A(3)**2)*SQRT(B(1)**2+B(2)**2+B(3)**2))
    IF (ABS(TT)>1.0) then
        if(TT<0.0) THEN
            THETA(I)=PI
        ELSEIF (TT>0.0) THEN
            THETA(I)=0.0
        endif
    ELSE
        THETA(I)=ACOS(TT)
    ENDIF
ENDDO
END SUBROUTINE getAngle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read the native dihedral angles and contact strength from DiAngle.txt
SUBROUTINE readDiAngle()    
CHARACTER    ::    FR*20,FS*3     

OPEN(301, FILE="GO_NativeState/DiAngle.txt")
DO I=1, NUM-3
    do j=1,4
        READ(301, '(A22,F10.6,A3,F12.6)') FR, EPS_DiTheta(i,j),FS,DiTheta0(i,j)
    enddo 
ENDDO 
DiTheta0=DiTheta0*PI/180.
!if(f_FDi_try==1) EPS_DiTheta=0.4*0.025*EPS_Theta
CLOSE(301)
END SUBROUTINE readDiAngle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! get the dihedral angles from particle positions
SUBROUTINE getDiAngle()    
real*8    ::    A(3),B(3),C(3)
real*8    ::    N(3), NN(3),SB,xt1,xt2

DO I=1, NUM-3
    A=X(I+1,:)-X(I,:)
    B=X(I+2,:)-X(I+1,:)
    C=X(I+3,:)-X(I+2,:)

    N=CrossProduct(A,B)
     NN=CrossProduct(B,C)

    SB=sqrt(B(1)**2+B(2)**2+B(3)**2)
    
    xt1=SB*sum(A*NN)
    xt2=sum(N*NN)
    DiTheta(i)=atan2(xt1,xt2)
    if(DiTheta(i)<0) DiTheta(i)=2.*PI+DiTheta(i)
ENDDO
END SUBROUTINE getDiAngle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read the native contact radius and contact strength from NativeContact.txt
SUBROUTINE readNativeContact()
integer        ::    i1,i2,i_tmp
integer,allocatable    :: j1(:),j2(:)
real*8        ::    eps_tmp,R0_tmp
CHARACTER    ::    FR1,FR2*4 
i_tmp=0
OPEN(302, FILE="GO_NativeState/Qlist.txt")
do while(.not.EOF(302))
    READ(302,*) 
    i_tmp=i_tmp+1
enddo
close(302)
n_nc=i_tmp
allocate(j1(n_nc),j2(n_nc))
OPEN(302, FILE="GO_NativeState/Qlist.txt")
do i=1,n_nc
    READ(302,*) j1(i),j2(i)
enddo
close(302)
Matrix=0;R0_nonNei=0.;EPS_nonNei=0.
OPEN(301, FILE="GO_NativeState/NativeContact.txt")
OPEN(303, FILE="test.dat")
i=1
DO while(.not.EOF(301))
    READ(301,'(A1,I4,A4,I4,F15.6,F12.6)') FR1,i1,FR2,i2,EPS_tmp,R0_tmp
    if(i2>num) i2=i2-num+1
    if(i2>num) i2=i2-num+4
    if(i1/=j1(i).or.i2/=j2(i)) cycle
    i=i+1
    Matrix(i1,i2)=2
    Matrix(i2,i1)=2
    EPS_nonNei(i1,i2)=-EPS_tmp*1.265414
    EPS_nonNei(i2,i1)=EPS_nonNei(i1,i2)
    R0_tmp=getABS(3,x(i1,:)-x(i2,:))
    Sigma_nonNei(i1,i2)=R0_tmp
    Sigma_nonNei(i2,i1)=R0_tmp
    write(303,*) i1,i2,EPS_nonNei(i1,i2),R0_tmp
ENDDO 
CLOSE(301)
CLOSE(303)
!pause
END SUBROUTINE readNativeContact
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read the native repulsive radius and contact strength from nonBond.txt
SUBROUTINE readNonBond()    
real*8        ::    R_half_repulsion(num)
CHARACTER    ::    FR*23     
real*8    ::    xtmp(3)

OPEN(301, FILE="GO_NativeState/nonbond.txt")
DO I=1, NUM
    READ(301, '(A23, F10.6)') FR, R_half_repulsion(i)
ENDDO 
CLOSE(301)
R0_nonNei=0.
do i=1,num-3
    do j=i+3,num
        xtmp=x(i,:)-x(j,:)
        R0_nonNei(i,j)=sqrt(xtmp(1)**2+xtmp(2)**2+xtmp(3)**2)    
        R0_nonNei(j,i)=R0_nonNei(i,j)    
        if(Matrix(i,j)==2) cycle
        Matrix(i,j)=1
        Matrix(j,i)=1
        EPS_nonNei(i,j)=0.000132
        EPS_nonNei(j,i)=0.000132
        Sigma_nonNei(i,j)=R_half_repulsion(i)+R_half_repulsion(j)
        Sigma_nonNei(j,i)=Sigma_nonNei(i,j)
    enddo
enddo !i
END SUBROUTINE readNonBond
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! get the nonNeighbor distance from particle positions. before using this subroutine, 
! one needs call readNativeContact and readNonBond first.
SUBROUTINE getNonNeighbor()    
real*8    ::    xtmp(3)

R_nonNei=0.
do i=1,num
    do j=i,num
        if(Matrix(i,j)==0) cycle
        xtmp=x(i,:)-x(j,:)
        R_nonNei(i,j)=sqrt(xtmp(1)**2+xtmp(2)**2+xtmp(3)**2)    
        R_nonNei(j,i)=R_nonNei(i,j)    
    enddo !j
enddo !i
END SUBROUTINE getNonNeighbor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine getR_Native_Repulsion(n_nc,RNative,RRP)
integer ::    N_NC
real*8    ::    RNative(N_NC),RRP(Num)

i_nc=0
do i=1,num
    RRP(i)=maxval(R_nonNei(i,:))
    do j=1,num
        if(matrix(i,j)==2.and.j>i) then
            i_nc=i_nc+1
            RNative(i_nc)=R_nonNei(i,j)
        endif
        if(matrix(i,j)==1.and.RRP(i)>R_nonNei(i,j)) RRP(i)=R_nonNei(i,j)
    enddo !j
enddo !i
RRP=RRP/2.
end subroutine getR_Native_Repulsion
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END module Go_Structure
