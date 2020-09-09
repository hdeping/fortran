program main
    use module_mct
    integer               :: itmp
    integer               :: jtmp
    real(8)               :: mat_final(ncut)
    
    filename = "final.txt"
    open(10,file = filename)


    call getmatrixes()
    mat_final = getfinal()
    do itmp = 1,ncut
        write(10,"(2f18.9)")dk(itmp),mat_final(itmp)
    end do
    close(10)


end program main
