program mpi_random 
use MPI 
integer :: ProcessID, numprocs, ierr, status(Mpi_Status_Size) 

integer :: l 
integer,allocatable :: iseed(:) 
real :: r 
real,allocatable :: s(:) 

call MPI_INIT(ierr) 
call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr) 
call MPI_COMM_RANK(MPI_COMM_WORLD, ProcessID, ierr) 

allocate(s(numprocs)) 
if(ProcessID == 0) then 
call random_seed() 
call random_number(s) 
s = s*(2**30) 
end if 
call MPI_BCAST(s, numprocs, MPI_REAL, 0, MPI_COMM_WORLD, ierr) 
call random_seed(SIZE = l) 
allocate(iseed(l)) 
iseed = floor(s(ProcessID+1)) 
call random_seed(PUT = iseed) 

call main(ProcessID)

deallocate(s) 
deallocate(iseed) 
call MPI_FINALIZE(ierr) 

end program mpi_random


subroutine main(myid)
use common_module
use Interface_module
use datalink_module
use dynamic_module
use Go_Structure

integer                ::    myid,jF,n_prev,nIF
real*8                ::    prob,M
real                ::    FF
type(IF_link),pointer :: p_prev,p1


! 初始化
jF=myid+1
F=dble(jF)*0.1
!F=0.1
call initial(jF)
print*,'Initial completes'
do while(p_IF%r<8.5)
    print*,"-----------------------------------------------------------------------"
    p_prev=>p_IF
    call newInterface(jF,p_prev)    ! 第0个界面
    p1=>head_IF
    write(fName,"('IF_',I2.2,'.dat')") jF
    open(10,file=fName)
    call outputLink(10,p1)
    close(10)
enddo
end subroutine main
