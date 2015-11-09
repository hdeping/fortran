program main
    use blas95
    implicit none
    include 'omp_lib.h'


    integer,parameter          :: n = int(1E7)
    real(8)                    :: a(n)
    real(8)                    :: t1
    real(8)                    :: t2
    real(8)                    :: tpara  = 0.0
    real(8)                    :: tsing  = 0.0
    real(8)                    :: resp   = 0.0
    real(8)                    :: ress   = 0.0
    real(8)                    :: res
    integer                    :: i 
    integer                    :: j 
    integer                    :: si
    
    call random_seed()

       !$omp parallel num_threads(4)
       !$omp do
    do j = 1,100
       call cpu_time(t1)
       !$omp do
       do i = 1, n
        call random_number(a(i))
       enddo
       res=asum(a)/dble(n)
       !$omp enddo
       call cpu_time(t2)
       resp  = resp + res
       tpara = tpara + t2 - t1

       call cpu_time(t1)
       do i = 1,n
           call random_number(a(i))
       end do
       call cpu_time(t2)
       res   = asum(a)/dble(n)
       ress  = ress + res
       tsing = tsing + t2 - t1
    end do
       !$omp enddo
       !$omp end parallel


    print *,"     res parallel  ==> ",resp
    print *,"     res single    ==> ",ress
    print *,"time cost parallel ==> ",tpara
    print *,"time cost single   ==> ",tsing
end program main
