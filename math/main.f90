program main
    use module_common
    implicit none

    real                h
    real                d
    real                theta
    real                r
    real                s

    write(*,*)"input h = "
    read(*,*)h
    write(*,*)"input d = "
    read(*,*)d

    r = (d**2.0 + h**2.0)/(2.0*h)
    theta = asin(d/r)

    s = theta*r**2.0 - d*(d**2.0 - h**2.0)/(2.0*h)

    write(*,*)"r is ",r
    write(*,*)"theta is ",theta
    write(*,*)"area is ",s


     
end program main
