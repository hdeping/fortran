module module_common
    implicit none
    real   ,parameter            :: alpha = 1.0
    real   ,parameter            :: l     = 1.0
    integer,parameter            :: n = 1000
    integer,parameter            :: m = 150
    integer                      :: i
    integer                      :: j
    integer                      :: k
    character(10)                :: filename 
    contains
!function getMeanDegree{{{
function getMeanDegree(beta)
    real,intent(in)     :: beta
    real                :: getMeanDegree
    ! status
    integer             :: a(n,n)
    integer             :: s(m,2)
    integer             :: summary
    ! call random_number 
    real                :: x1,x2,x3,p
    real                :: dis  ! distance

    ! get m random sites
    call random_seed() 
    do i = 1,m
        call random_number (x1)
        j = int(n*x1) + 1
        s(i,1) = j
        call random_number (x1)
        j = int(n*x1) + 1
        s(i,2) = j
    enddo !cycle ends
    
    ! get random graph

    do   i = 1,10000000
         call random_number (x1)
         j = int(m*x1) + 1
         call random_number (x1)
         k = int(m*x1) + 1
         dis = sqrt((s(j,1) - s(k,1))**2.0 + &
                    (s(j,2) - s(k,2))**2.0) 
         call random_number (x1)

         !if ( a(j,k) == 1 )then
         !    cycle
         !endif ! if ends
         
         p = alpha*exp(-dis/(beta*l))
         
         if ( x1 < p )then
             a(j,k) = 1
         else
             a(j,k) = 0
         endif ! if ends

    enddo !cycle ends
     
    summary = 0
    do i = 1,n
        do j = 1,n
            summary = summary + a(i,j)
        enddo !cycle ends
         
    enddo !cycle ends

    getMeanDegree = summary*1.0 / m
 
end function getMeanDegree
!}}}
        
end module module_common
