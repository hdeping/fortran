program main
    use module_algebra
    implicit none 

    integer  n
    real,allocatable     ::  a(:,:)
    real,allocatable     ::  b(:,:)
    real,allocatable     ::  re(:)
    real                 ::  deter
    filename="data.txt"
    open(10,file=filename,status="old",iostat=ierror)
    !  first get the number of rank  ( the first line of the input)
    read(10,"(i3)",iostat=ierror)n
    allocate(a(n,n))
    allocate(re(n))
    do i=1,n
        read(10,*,iostat=ierror)a(i,:)
        if(ierror/=0)exit
    end do
    close(10)
    call eigen(re,a,n)
    write(*,*)re
    call det(deter,a,n)
    write(*,"('the determinant is ==>',f10.5)")deter
    deallocate(a)
    deallocate(re)
    end program main

