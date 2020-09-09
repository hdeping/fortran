program main
    use module_common
    use omp_lib
    implicit none
    real(8)             a(n,n)
    real(8)             b(n,n)
    real(8)             c(n,n)

    !$omp parallel
    do i = 1,10
        j = omp_get_num_threads()
        write(*,*)"Hello,world"
        write(*,*)"j = ",j
        call sleep(10)
    enddo !cycle ends
    !$omp end parallel

end program main
