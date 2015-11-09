program test
    use module_algebra
    use module_mct
    real(8)               :: tmpa(n)
    real(8)               :: tmpb(n)
    real(8)               :: tmpc(n)
    
    do i = 1,n
        if(i < 200)then
            tmpa(i) = i*deltar
        else
            tmpa(i) = 0.0
        endif
    end do
    ! get fst
    tmpb = fst(tmpa,1)
    tmpc = fst(tmpb,- 1)
    do i = 1,n
        print "(3f12.6)",tmpa(i),tmpb(i),tmpc(i)
        if(mod(i,10) == 0)pause
    end do
end program test
!code before {{{
    !call getmatrixes()
    !print *,"sk =           , cr = "
    !do i = 290,300
    !    print "(8f12.6)",sk(1,:,i),sk(2,:,i),cr(1,:,i),cr(2,:,i)
    !end do


!}}}
