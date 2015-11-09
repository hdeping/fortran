program main
    use module_common
    implicit none


    do i = 1,n
        a(i) = i
    end do

    print *,a(200000:200010)

    print *,sign(200)
    
end program main
