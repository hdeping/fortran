program main
    use module_common
    implicit none

    do ii = 1,n
        x(ii) = dble(ii)
    end do
    val = ave_sum(x)
    print *,val
    val = ave_mul(x)
    print *,val
end program main
