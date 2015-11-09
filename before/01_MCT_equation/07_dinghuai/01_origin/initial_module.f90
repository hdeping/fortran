module initial_module
use parameter_module
use function_module
use communicate_module

contains

!--********************************************--!

subroutine initial_parameter()
integer                    ::                i, j, k
real*8                    ::                cutr

!beta=1.0
v0=0.0
!cutr=12.5663
!dr=cutr/dble(n); dk=(pi/dble(n))/dr 
dr=0.1; dk=(pi/dble(n))/dr
h=dk; rho=1.1
dt=1.0; dro=3.0*dt

!--***************************--!

do i=0, node-1
    inode(i, 1)=nnode*i
    inode(i, 2)=nnode*(i+1)-1
enddo

allocate(dphi(inode(senode, 1):inode(senode, 2), 0:sn))
allocate(ke_me(inode(senode, 1):inode(senode, 2), 0:sn))

if (senode==0) then
    open(1, file='sk.txt', action='read')
        do i=0, n-1
            read(1, '(i5, f18.8)') k, sk(i)
        enddo
    close(1)

    open(1, file='ck.txt', action='read')
        do i=0, n-1
            read(1, '(i5, f18.8)') k, ck(i)    
        enddo
    close(1)
    
    do i=0, n-1
        gammak(i)=dt*((h*dble(i))**2.)/sk(i)
        ak(i)=1.0+sk(i)*(v0**2.)/(dt*dro)
    enddo
    gammak(0)=0.0

endif

call comm_initial_parameter()

det_time=10.0d0**-5.

end subroutine initial_parameter

!--********************************************--!

end module initial_module
