program main
    use module_common
    use module_algebra
    implicit none

    a(1,:) = (/1.0,2.0,3.0,4.0/)
    a(2,:) = (/5.0,7.0,9.0,8.0/)
    a(3,:) = (/3.0,3.0,7.0,4.0/)
    a(4,:) = (/5.0,5.0,0.0,6.0/)
    b      = (/1.0,1.0,- 8.0, 1.0/)
    
    c = sol_equ(a,b,n)
    print *,c


end program main
