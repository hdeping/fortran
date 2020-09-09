program main
    use module_common
    implicit none
    real(8)           :: force(3)
    real(8)           :: a(3) 
    real(8)           :: b(3) 
    real(8)           :: t1
    real(8)           :: t2
    real(8)           :: force_new


    a = (/1.0,3.0,4.0/)
    b = (/2.0,6.0,9.0/)

    call cpu_time(t1)
    force  = getintegral(a,b)
    call cpu_time(t2)
    print *,"time cost is ",t2 - t1
    print *,"force = ",force
    force_new = sqrt((force(1))**2.0 + (force(2))**2.0 + (force(1))**2.0)
    print *,"scalar force is ",force_new

end program main
