module math_module

use common_module
use datalink_module

contains
!----------------------------------------------------------------------------
subroutine rotation(num,v0,x)
implicit    none
integer        ::    num,i,j
real(8)        ::    x(num,3),v0(3),r(num-1),v(num-1,3)
real(8)        ::    axle(3),angle(2),radius,rot(3,3)
do i = 1,num-1
    v(i,:) = x(i+1,:)-x(1,:)
    r(i) = getABS(3,v(i,:))
    v(i,:) = v(i,:)/r(i)
enddo !i
axle = CrossProduct(v(num-1,:),v0)
radius = getABS(3,axle)
if(radius<1.E-4) return
axle = axle/radius
angle(2) = radius
angle(1) = sum(v0*v(num-1,:))
call RotationMatrix(axle,angle,rot)
do i = 1,num-1
    axle = v(i,:)
    do j = 1,3
        v(i,j) = sum(rot(j,:)*axle)
    enddo
enddo !i
do i = 1,3
    x(2:num,i) = r*v(:,i)+x(1,i)
enddo
end subroutine
!----------------------------------------------------------------------------
Integer function InMSSA(op)
real*8    op
InMSSA = 0
if(dble(op_drct)*(op-MSSA)<head_IF%next%r) InMSSA = 1
return
end function
!----------------------------------------------------------------------------
Integer function InIF(op,p)
real*8    op
type(IF_link),pointer    ::    p
InIF = 0
if(dble(op_drct)*(op-MSSA)<p%r) InIF = 1
return
end function
!**********************************************************************************************
subroutine getDis_R(nT,r,np_r,dis_r)
! 计算r中存储的数据的分布函数和累计分布函数
implicit    none
integer                ::    nT,np_r    
real*8                ::    r(nT),rMax,rMin
real*8                ::    dis_r(np_r,3)
integer                ::    i,i_tmp

dis_r = .0
rMax = maxval(r)
rMin = minval(r)
if(rMax <=  rMin) goto 300
if(nT<np_r) return !数据过少，分布取0
do i = 1,nt
    i_tmp = ceiling(dble(np_r)*(r(i)-rMin)/(rMax-rMin))
    if(i_tmp<1) i_tmp = 1
    if(i_tmp>np_r) i_tmp = np_r
    dis_r(i_tmp,2) = dis_r(i_tmp,2)+1.
end do !i
dis_r(:,2) = dis_r(:,2)/dble(nT)
300 do i = 1,np_r
    dis_r(i,1) = (dble(i)-0.5)*(rMax-rMin)/dble(np_r)+rMin
    dis_r(i,3) = sum(dis_r(1:i,2))
enddo !i
end subroutine getDis_R
!**********************************************************************************************
! 其他数学子程序,包括:
!                getABS: 求d维空间向量的模
!                GaussianNoise: Function to Return a Gaussian White Noise 
!                CrossProduct: function to return cross product of two vector
!**********************************************************************************************
real*8 function getABS(d,x)
! Get abs(x), x is a d-dimention vector
integer        ::    d
real*8        ::    x(d),a
a = 0.0
do i = 1,d
    a = a+x(i)**2
enddo !i
if(a<.0) then
    print*,'here!',a
    pause
endif
getABS = sqrt(a)
return
end function getABS
!-------------------------------------------------------------------------
real*8 Function GaussianNoise(u,g)
! Function to Return a Gaussian Noise with Mean u and Variance g
implicit    none
real*8        :: u,g 
Integer        :: i,n
real*8        :: b,temp

n = 12 ; b = 0.0 
do i = 1,n
   call random_number(temp)
   b = b+temp
end do
GaussianNoise = u+g*(b-n/2)
Return
End Function GaussianNoise
!-------------------------------------------------------------------------
! function to return cross product of v1 and v2
function CrossProduct(v1,v2)
REAL*8, DIMENSION(3) :: CrossProduct
real*8    v1(3),v2(3)

CrossProduct(1) = v1(2) * v2(3) - v1(3) * v2(2)
CrossProduct(2) = v1(3) * v2(1) - v1(1) * v2(3)
CrossProduct(3) = v1(1) * v2(2) - v1(2) * v2(1)
return
end function CrossProduct
!-------------------------------------------------------------------------
! 给定旋转轴和旋转角度求旋转矩阵
subroutine RotationMatrix(axle,angle,rot)
implicit    none
real*8        :: axle(3),angle(2),rot(3,3) 

Rot(1,1) = angle(1) + axle(1)*axle(1)*(1.0-angle(1))
Rot(2,1) = axle(2)*axle(1)*(1.0-angle(1)) + axle(3)*angle(2)
Rot(3,1) = axle(3)*axle(1)*(1.0-angle(1)) - axle(2)*angle(2)
        
Rot(1,2) = axle(1)*axle(2)*(1.0-angle(1)) - axle(3)*angle(2)
Rot(2,2) = angle(1) + axle(2)*axle(2)*(1.0-angle(1))
Rot(3,2) = axle(3)*axle(2)*(1.0-angle(1)) + axle(1)*angle(2)
        
Rot(1,3) = axle(1)*axle(3)*(1.0-angle(1)) + axle(2)*angle(2)
Rot(2,3) = axle(2)*axle(3)*(1.0-angle(1)) - axle(1)*angle(2)
Rot(3,3) = angle(1) + axle(3)*axle(3)*(1.0-angle(1))
return
end subroutine RotationMatrix
!**********************************************************************************************
end module math_module
