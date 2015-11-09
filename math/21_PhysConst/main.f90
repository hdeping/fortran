program main
    use module_common
    implicit none
    real(8)     :: s
    real(8)     :: tmp
    filename = "data.txt"
    open(10,file = filename)

    r = 5.3D-11
    v = sqrt(ev**2.0/(4.0*pi*epsi0*emass*r))
    print *,"v = ",v
    s = ev**2.0/(4.0*pi*epsi0*emass*r*c**2.0)
    print *,"s = ",s
    do i = 1,118000
        tmp = dble(i)*s
        v   = sqrt(tmp)
        tmp = tmp**2.0
        tmp = (- tmp + sqrt(tmp**2.0 + 4*tmp))/2.0
        x   = sqrt(tmp)
        if ( mod(i,1000) == 0 )then
            write(10,*)i,v,x
        endif ! if ends
        !print *,"i = ",i,"x = ",x
    end do
    close(10)
    !s = g*emass**2.0/(r**2.0)
    !print *,"force is ",s
    !v = sqrt(13.6*ev*r/emass)
    !print *,"v = ",v
end program main
