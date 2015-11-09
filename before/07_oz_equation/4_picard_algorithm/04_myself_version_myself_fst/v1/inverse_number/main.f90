program main
    implicit none
    integer,parameter         :: l = 10
    integer,parameter         :: n = 2**l
    integer                   :: i
    integer                   :: j
    integer                   :: k
    integer                   :: num  
    integer                   :: tmp
    integer                   :: a(n)
    real                      :: t1
    real                      :: t2
    character(10)             :: filename 
    filename = "time.txt"
    open(10,file = filename)

    call cpu_time(t1)
    call inverse1(a,n,l)
    call cpu_time(t2)
    print *,"time cost is ==>",t2 - t1
    i = 0
    do i = 1,10
        write(*,"(2i10,2a<l+1>)")i,a(i),serial(i - 1),serial(a(i) - 1)
    end do
    close(10)
    pause
    contains
!subroutine inverse1{{{
subroutine inverse1(a,n,l)
    integer,intent(in)    :: l
    integer,intent(in)    :: n
    integer,intent(out)   :: a(n)
    integer               :: b(n)
    integer               :: c(n)
    integer               :: i

    b(1) = 2
    do i = 2,n
        b(i) = b(i-1)*2
    end do
    do i = 1,n-1
        c(i) = b(l-i)
    end do
    c(n) = 1
    a(1) = 1
    a(n) = b(l) 
    do i = 1,n/2 -1
        if(mod(i,2) == 1)then
            a(i+1) = a(i) + b(l-1)
        else
            num = i - 1
            do j = 1,l-1
               k = mod(num,2)
               if(k == 0) exit
               num = num/2
            end do
            j = j-1
            !print *,"j = ",j
            a(i+1) = a(i) - sum(c(1:j)) + c(j+1)
        endif
        j = i + n/2
        a(j) =a(i) + 1
    end do
end subroutine inverse1
!}}}
!subroutine inverse2{{{
subroutine inverse2(a,n,l)
    integer,intent(in)    :: l
    integer,intent(in)    :: n
    integer,intent(out)   :: a(n)
    integer               :: b(n)
    integer               :: c(n)
    
    a(1) = 1
    a(n) = n
    do i = 2,n/2
        tmp = i - 1
        a(i) = 0
        do while(tmp /= 0) 
            num = mod(tmp,2) 
            a(i) = 2*a(i) + num
            tmp = tmp/2
        end do
        a(i) = a(i) + 1
        j = i + n/2
        a(j) = a(i) + 1
    end do
    
end subroutine inverse2
!}}}
!function serial{{{
function serial(k)
    integer,intent(in)     :: k
    character(l)           :: serial
    integer                :: ii
    integer                :: jj
    num = k
    do  ii = 1,l
         tmp = mod(num,2) 
         num = num/2
         jj = l + 1 - ii
         if(tmp == 1)then
             serial(jj:jj) = '1'
         else
             serial(jj:jj) = '0'
         endif
    end do
end function serial
!}}}
end program main
