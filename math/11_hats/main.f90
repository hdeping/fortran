program main
    use module_common
    implicit none
    integer         :: a(n)
    integer         :: times

    a = 0
    do i = 1,30
        do j = n,1,-1
            if(a(j) == 0)then
                times    = j
                a(times) = 1
                exit
            endif
        end do
        times = times + 1
        do j = times,n
            a(j) = 0
        end do
        print "(<n>i1)",a(:)
    end do
end program main
