module potential_module

use common_module
use Go_Structure

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! get the total potential
subroutine getU_tot(v)
real*8    :: v
real*8    :: v_bond(num-1),v_theta(num-2)
real*8    :: v_ditheta(num-3),v_nonNei(num,num)

call getStructureVariable()
call getU_bond(v_bond)
call getU_Theta(v_theta)
call getU_DiTheta(v_ditheta)
call getU_nonNeighbor(v_nonNei)

v=sum(v_bond)+sum(v_theta)+sum(v_ditheta)+sum(v_nonNei)
end subroutine getU_tot
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Bond potential
SUBROUTINE getU_bond(u)
real*8    ::    u(num-1)

u=EPS_bond*(R_bond-R0_bond)**2

END SUBROUTINE getU_bond
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Bending potential
SUBROUTINE getU_Theta(u)
real*8    ::    U(num-2)

U=0.5*EPS_Theta*(Theta-Theta0)**2

END SUBROUTINE getU_Theta
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Dihedral potential
SUBROUTINE getU_DiTheta(u)
real*8    ::    U(num-3)

U=0.0
DO i=1,4
!    U=U+EPS_DiTheta(:,i)/dble(i)*(1.-cos(dble(i)*DiTheta-DiTheta0(:,i)))
    U=U+EPS_DiTheta(:,i)/dble(i)*(1.+cos(dble(i)*(DiTheta-DiTheta0(:,i))))
    if(f_FDi_try==1) then
        U=0.4*0.025*EPS_Theta*(1.-cos(DiTheta-DiTheta01))
        exit
    endif        
ENDDO

END SUBROUTINE getU_DiTheta
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Native and repulsion potential
SUBROUTINE getU_nonNeighbor(u)
real*8    ::    U(num,num)

u=0.0
where(Matrix==1.and.R_nonNei<=Sigma_nonNei)
    u=EPS_nonNei*(Sigma_nonNei/R_nonNei/1.12246)**12
elsewhere(Matrix==2) 
    u=EPS_nonNei*(13.*(Sigma_nonNei/R_nonNei)**12- &
        18.*(Sigma_nonNei/R_nonNei)**10+4.*(Sigma_nonNei/R_nonNei)**6)
endwhere
END SUBROUTINE getU_nonNeighbor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module potential_module
