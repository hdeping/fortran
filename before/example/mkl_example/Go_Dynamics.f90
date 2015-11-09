module Go_dynamics
use common_module
use math_module
use Go_Structure
use Go_Force

contains
!*****************************************************************************************
SUBROUTINE RunDynamic()
! Dynamic evolution for one time step
! using the velocity form of the Verlet algorithm (J. Chem. Phys. 76， 637)
! and SHAKE algorithm to insure bond length
real*8    ::    F0_tot(num,3),vct0(num-1,3)
real*8    ::    Z(num*3,num*3),A(num*3,num*3)
real*8    ::    dx0(num*3),dx1(num*3),dx2(num*3)
real*8    ::    noise(num*3)

do itry=1,nt_jump
    call getStructureVariable()
    ! 保存前一步键向量用来进行shake
    do i=1,3
        vct0(:,i)=vct_bond(:,i)*R_bond(i)
    enddo !i
    ! 更新位置
    DO I=1, NUM
        X(I,:)=X(I,:)+V(i,:)*dt+F_tot(i,:)*dt*dt/(2.*mass(i))           
    ENDDO
    F0_tot=F_tot !!保存前一时刻的受力
    ! 采用shake算法保证键长
    if(f_shake==1) CALL SHAKE(x,vct0)
    if(f_fixend==1) then
        do i=1,3
            x(:,i)=x(:,i)-x(1,i)+x0(1,i)
        enddo !i        
    else
        do i=1,3
            x(:,i)=x(:,i)-sum(x(:,i))/dble(num)+sum(x0(:,i))/dble(num)
        enddo !i
    endif
    call rotation(num,-A0,x)
    ! 更新速度
    call getStructureVariable()
    call getF_tot(F_tot)
    call getTensor(gamma,x,Z,A)
    do i=1,num
        do j=1,3
            noise(i*3-3+j)=GaussianNoise(dble(0.),dble(1.))/sqrt(dt)
        enddo !j
        IF (f_pullingN==1.and.I==1) F_tot(i,:)=F_tot(i,:)+F*A0    
        IF (f_pullingC==1.and.I==Num) F_tot(i,:)=F_tot(i,:)-F*A0    
    enddo !i
    do i=1,num
        do j=1,3
            F_tot(i,j)=F_tot(i,j)+sum(A(i*3-3+j,:)*noise)
            dx0(i*3-3+j)=v(i,j)+dt*(F0_tot(i,j)+F_tot(i,j))/(2.*mass(i))
        enddo
    enddo 
    do i=1,3*num
        dx1(i)=-sum(Z(i,:)*dx0)*dt/2.
    enddo 
    do i=1,3*num
        dx2(i)=sum(Z(i,:)*dx1)*dt*dt/4.
    enddo
    do i=1,num
        do j=1,3
            V(i,j)=dx0(i*3-3+j)+dx1(i*3-3+j)/mass(i)+dx2(i*3-3+j)/mass(i)/mass(i)
            noise(i*3-3+j)=v(i,j) ! 暂存
        enddo !j
    enddo !i
    do i=1,num
        do j=1,3
            F_tot(i,j)=F_tot(i,j)-sum(Z(i*3-3+j,:)*noise)
        enddo
    enddo 

enddo
END SUBROUTINE RunDynamic
!----------------------------------------------------------------------------
SUBROUTINE RunDynamic_noHI()
! Dynamic evolution for one time step
! using the velocity form of the Verlet algorithm (J. Chem. Phys. 76， 637)
! and SHAKE algorithm to insure bond length
real*8    ::    F0_tot(num,3),vct0(num-1,3)
real*8    ::    noise0(num,3)
real*8    ::    AA,BB,CC

do itry=1,nt_jump
    call getStructureVariable()
    ! 保存前一步键向量用来进行shake
    do i=1,3
        vct0(:,i)=vct_bond(:,i)*R_bond(i)
    enddo !i
    ! 更新位置
    DO I=1, NUM
        X(I,:)=X(I,:)+V(i,:)*dt+(F_tot(i,:)-gamma*V(i,:))*dt*dt/(2.*mass(i))           
    ENDDO
    F0_tot=F_tot !!保存前一时刻的受力
    ! 采用shake算法保证键长
    if(f_shake==1) CALL SHAKE(x,vct0)
    if(f_fixend==1) then
        do i=1,3
            x(:,i)=x(:,i)-x(1,i)+x0(1,i)
        enddo !i        
    else
        do i=1,3
            x(:,i)=x(:,i)-sum(x(:,i))/dble(num)+sum(x0(:,i))/dble(num)
        enddo !i
    endif
    call rotation(num,-A0,x)
    ! 更新速度
    call getStructureVariable()
    call getF_tot(F_tot)
    AA=0.5*dt/mass(i)
    BB=1.-AA*gamma
    CC=BB+(AA*gamma)**2
    do i=1,num
!        if(f_fixEnd==1.and.i==num) cycle
        do j=1,3
            F_tot(i,j)=F_tot(i,j)+SQRT(2.0*TEMP*gamma*kB/dt)*GaussianNoise(dble(0.),dble(1.))
        enddo !j
        IF (f_pullingN==1.and.I==1) F_tot(i,:)=F_tot(i,:)+F*A0    
        IF (f_pullingC==1.and.I==Num) F_tot(i,:)=F_tot(i,:)-F*A0    
        v(i,:)=BB*CC*v(i,:)+AA*CC*(F0_tot(i,:)+F_tot(i,:))
    enddo !i
enddo
END SUBROUTINE RunDynamic_noHI
!----------------------------------------------------------------------------
subroutine getOP(op)
real*8        ::    op
op=sqrt(sum((x(1,:)-x(num,:))**2))
end subroutine getOP
!----------------------------------------------------------------------------
SUBROUTINE SHAKE(x, vct0)
! Using SHAKE alogrithm to enforce the bond length constrains
! x, 输入：当前时刻未约束坐标，输出：约束后坐标
! vct，输入：前一时刻的键向量
real*8    ::    x(num,3),vct0(num-1,3)
real*8    ::    dR(num-1)
real*8    ::  L(NUM-1)
call getBond()
dR=abs(R_bond**2-R0_bond**2)/(2.*R0_bond**2)
DO WHILE (maxval(dR)>1.E-4)
    DO I=1, NUM-1        
        L(I)=0.25*(R_bond(I)**2-R0_bond(I)**2)/sum(vct0(i,:)*vct_bond(i,:)*R_bond(i))/(1./mass(i)+1./mass(i+1)) *1.E-1
    ENDDO
!    if(f_pullingC==1) x(num,:)=x(num,:)-mass(num)*L(num-1)*vct0(num-1,:)
!    if(f_pullingN==1) x(1,:)=x(1,:)+mass(1)*L(1)*vct0(1,:)
    x(num,:)=x(num,:)-mass(num)*L(num-1)*vct0(num-1,:)
    x(1,:)=x(1,:)+mass(1)*L(1)*vct0(1,:)
    do i=2,num-1
        x(i,:)=x(i,:)+mass(i)*(L(i)*vct0(i,:)-L(i-1)*vct0(i-1,:))
    enddo !i
    call getBond()
    dR=abs(R_bond**2-R0_bond**2)/(2.*R0_bond**2)
ENDDO
END SUBROUTINE SHAKE

!----------------------------------------------------------------------------
subroutine getTensor(CoE_friction,x,Z,A)
!use imsl
use lapack95
real*8    ::    CoE_friction,x(num,3)
real*8    ::    M(num*3,num*3),Z(num*3,num*3),A(num*3,num*3)
real*8    ::    e(3),m1(3,3),aij,rij,mu
real*8    ::    r_bead=0.5
real*8    ::    Un(3,3)=(/(/1.,0.,0./),(/0.,1.,0./),(/0.,0.,1./)/)
integer    ::    ix,iy

mu=1./CoE_friction
Z=0.0;A=0.0
if(f_HI==0) then
    do i=1,3*NUM
        Z(i,i)=CoE_friction
        A(i,i)=sqrt(2.*kB*temp*CoE_friction)
    enddo
    return
endif
do i = 1, Num
    ix=3*i-2
    do j = i, Num
        iy=3*j-2
        if(i==j) then
            M(ix:ix+2,iy:iy+2)=mu*Un
        else
            ! 求Mobility矩阵的矩阵元
            e=x(j,:)-x(i,:)
            rij=getabs(3,e)
            e=e/rij
            m1(1,1)=e(1)*e(1);m1(1,2)=e(1)*e(2);m1(1,3)=e(1)*e(3)
            m1(2,1)=e(2)*e(1);m1(2,2)=e(2)*e(2);m1(2,3)=e(2)*e(3)
            m1(3,1)=e(3)*e(1);m1(3,2)=e(3)*e(2);m1(3,3)=e(3)*e(3)
            if(rij<2.*r_bead) then
                M(ix:ix+2,iy:iy+2)=mu*((1.-9.*rij/(32.*r_bead))*Un+(3.*r_bead/(32.*rij))*m1)
                M(iy,ix)=M(ix,iy);M(iy,ix+1)=M(ix+1,iy);M(iy,ix+2)=M(ix+2,iy)
                M(iy+1,ix)=M(ix,iy+1);M(iy+1,ix+1)=M(ix+1,iy+1);M(iy+1,ix+2)=M(ix+2,iy+1)
                M(iy+2,ix)=M(ix,iy+2);M(iy+2,ix+1)=M(ix+1,iy+2);M(iy+2,ix+2)=M(ix+2,iy+2)
            else
                M(ix:ix+2,iy:iy+2)=mu*(0.75*(r_bead/rij)*(Un+m1)+0.5*((r_bead/rij)**3)*(Un-3.*m1))
                M(iy,ix)=M(ix,iy);M(iy,ix+1)=M(ix+1,iy);M(iy,ix+2)=M(ix+2,iy)
                M(iy+1,ix)=M(ix,iy+1);M(iy+1,ix+1)=M(ix+1,iy+1);M(iy+1,ix+2)=M(ix+2,iy+1)
                M(iy+2,ix)=M(ix,iy+2);M(iy+2,ix+1)=M(ix+1,iy+2);M(iy+2,ix+2)=M(ix+2,iy+2)
            endif
        endif
    enddo !ij
enddo !j

A=M
call potrf(A,uplo='L')
do i=1,3*num-1
    A(i,i+1:)=0.
enddo
Z=A
call potri(Z,uplo='L')
do i=1,3*num-1
    Z(i,i+1:)=Z(i+1:,i)
enddo
A=Z
call potrf(A,uplo='L')
do i=1,3*num-1
    A(i,i+1:)=0.
enddo
end subroutine getTensor
!*************************************************************************************************************************
end module Go_dynamics
