MODULE Go_Force
use common_module
use math_module
USE Go_Structure

CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! get the total force
subroutine getF_tot(F)
real*8    :: F(num,3)
real*8    :: F_bond(num,3),F_theta(num,3)
real*8    :: F_ditheta(num,3),F_nonNei(num,3)

call getStructureVariable()
call getF_bond(num,F_bond)
call getF_Theta(num,F_theta)
call getF_DiTheta(num,F_ditheta)
call getF_nonNeighbor(num,F_nonNei)

if(f_printV==1) then
    print*,'F_bond:',maxval(F_bond)
    print*,'F_theta:',maxval(F_theta)
    print*,'F_DiTheta:',maxval(F_DiTheta)
    print*,'F_nonNei:',maxval(F_nonNei)
!    write(111,"(6G18.9)")Theta(1),Theta0(1), F_theta(1,1),DiTheta(1),DiTheta01(1),F_DiTheta(1,1)
!    pause
endif 

F=F_bond+F_theta+F_ditheta+F_nonNei
end subroutine getF_tot
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Bond force on each bead
SUBROUTINE getF_bond(n,F_bond)
integer    ::    n
real*8    ::    F_bond(n,3),F(n-1)

F=-0.5*EPS_bond*(R_bond-R0_bond)

F_bond(1,:)=-F(1)*vct_bond(1,:)
do i=2,n-1
    F_bond(i,:)=F(i-1)*vct_bond(i-1,:)-F(i)*vct_bond(i,:)
enddo 
F_bond(n,:)=F(n-1)*vct_bond(n-1,:)
END SUBROUTINE getF_bond
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Bending force
SUBROUTINE getF_Theta(n,F_theta)
real*8    ::    F_theta(n,3),F(n-2)
real*8    ::    p(3),d(3)

F=-EPS_Theta*(Theta-Theta0)

F_theta=0.
do i=1,n-2
    p=crossproduct(vct_bond(i,:),vct_bond(i+1,:))
    d=crossproduct(p,vct_bond(i,:))
    F_theta(i,:)=F_theta(i,:)-F(i)*d/sqrt(d(1)**2+d(2)**2+d(3)**2)
    F_theta(i+1,:)=F_theta(i+2,:)+F(i)*d/sqrt(d(1)**2+d(2)**2+d(3)**2)

    d=crossproduct(vct_bond(i+1,:),-p)
    F_theta(i+2,:)=F_theta(i+2,:)-F(i)*d/sqrt(d(1)**2+d(2)**2+d(3)**2)
    F_theta(i+1,:)=F_theta(i+1,:)+F(i)*d/sqrt(d(1)**2+d(2)**2+d(3)**2)
enddo !i
END SUBROUTINE getF_Theta
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Dihedral potential
SUBROUTINE getF_DiTheta(n,F_DiTheta)
real*8    ::    F_DiTheta(n,3),F(n-3)
real*8    ::    d(3)

F=0.
DO i=1,4
!    F=F+dble(i)*EPS_DiTheta(:,i)*(sin(dble(i)*DiTheta-DiTheta0(:,i))) 
!    F=F-EPS_DiTheta(:,i)*(sin(dble(i)*DiTheta-DiTheta0(:,i))) 
!    F=F-EPS_DiTheta(:,i)*(sin(dble(i)*(DiTheta-DiTheta0(:,i)))) 
    if(f_FDi_try==1) then
        F=-3.*EPS_DiTheta(:,2)*(sin(3.*(DiTheta-DiTheta01)))
        exit
    endif        
ENDDO

F_DiTheta=0.
do i=1,n-3
    d=crossproduct(vct_bond(i+1,:),vct_bond(i,:))    
    F_DiTheta(i,:)=F_DiTheta(i,:)+F(i)*d/sqrt(d(1)**2+d(2)**2+d(3)**2)
    F_DiTheta(i+1,:)=F_DiTheta(i+1,:)-F(i)*d/sqrt(d(1)**2+d(2)**2+d(3)**2)

    d=crossproduct(vct_bond(i+1,:),vct_bond(i+2,:))    
    F_DiTheta(i+3,:)=F_DiTheta(i+3,:)+F(i)*d/sqrt(d(1)**2+d(2)**2+d(3)**2)
    F_DiTheta(i+2,:)=F_DiTheta(i+2,:)-F(i)*d/sqrt(d(1)**2+d(2)**2+d(3)**2)
enddo !i
END SUBROUTINE getF_DiTheta
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Native and repulsion potential
SUBROUTINE getF_nonNeighbor(n,F_nonNei)
real*8    ::    F_nonNei(n,3),F(n,n)
real*8    ::    d(3)

F=0.0
where(Matrix==1.and.R_nonNei<=Sigma_nonNei)
!    F=12.*EPS_nonNei/R_nonNei*(Sigma_nonNei/R_nonNei/1.12246)**12
    F=EPS_nonNei/R_nonNei*(12.*(Sigma_nonNei/R_nonNei/1.12246)**12- &
        6.*(Sigma_nonNei/R_nonNei/1.12246)**6)
endwhere
where(Matrix==2) 
    F=EPS_nonNei/R_nonNei*(156.*(Sigma_nonNei/R_nonNei)**12- &
        180.*(Sigma_nonNei/R_nonNei)**10+24.*(Sigma_nonNei/R_nonNei)**6)
endwhere

F_nonNei=0.
do i=1,n-3
    do j=i+1,n
        d=x(j,:)-x(i,:)
!        if(i==1.and.j==4) then
        F_nonNei(i,:)=F_nonNei(i,:)-F(i,j)*d/sqrt(d(1)**2+d(2)**2+d(3)**2)
        F_nonNei(j,:)=F_nonNei(j,:)+F(i,j)*d/sqrt(d(1)**2+d(2)**2+d(3)**2)
!        endif
    enddo !j
enddo !i

END SUBROUTINE getF_nonNeighbor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


END Module Go_Force
