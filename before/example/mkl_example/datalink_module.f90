module datalink_module
use common_module

contains
!*************************************************************************************
subroutine Outputlink(ID,p)
! ID为输出文件接口，p为指向界面结构体头指针的指针
integer                   ::    ID
type(IF_link),pointer     ::    p
real*8                    ::    k,k1
k=1.
k1=1.
do while(associated(p%next))
    p=>p%next
    k1=k1*p%cp(1)
    k=k*sum(p%cp)
    write(ID,"(I5,6G18.9)") p%i_IF,p%MSSA,p%r,p%cp,k1,k
enddo
end subroutine Outputlink
!-------------------------------------------------------------------------
subroutine Inputlink(ID,i_IF,p)
! ID为输出文件接口，p为指向界面结构体头指针的指针
implicit    none
integer                   ::    ID,i_IF,i
type(IF_link),pointer     ::    p
real*8                    ::    k

do i=1,i_IF
    call NewLinkUnit(p)
    read(ID,"(I5,6G18.9)") p%i_IF,p%MSSA,p%r,p%cp,k
enddo
end subroutine Inputlink
!-------------------------------------------------------------------------
subroutine Creatlink(head)
type(IF_link),pointer    :: head
if(associated(head)) call Deletelink(head)
allocate(head,stat=err)
if(err/=0) then
    print*,'Out of memory!'
    stop
endif
nullify(head%prev)
nullify(head%next)
end subroutine Creatlink
!-------------------------------------------------------------------------
subroutine Deletelink(list)
type(IF_link),pointer        :: next,list
do while(associated(list))
    next=>list%next
    deallocate(list)
    list=>next
enddo
end subroutine Deletelink
!-------------------------------------------------------------------------
subroutine NewlinkUnit(list)
type(IF_link),pointer    :: list
allocate(list%next,stat=err)
if(err/=0) then
    print*,'Out of memory!'
stop
endif
list%next%prev=>list
list=>list%next
nullify(list%next)
end subroutine NewlinkUnit
!*************************************************************************************
end module datalink_module
