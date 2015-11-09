module module_growth
    use  module_common
    use module_lattice
    
    contains
!******************************************************
    
!******************************************************
    subroutine getbordersite()
    integer x,y,count
    
    do i = 1,n
    do j = 1,m
    
    end do  !j
    end do  !i
    
    end subroutine getbordersite
!******************************************************
    !  j for the multiscale  ,1 for the whole type , and 6 only for C6
    subroutine getfluxtype(i,j)
    integer  i,j
    real   proba(6),a(6),tmp
    integer  k
    a(1) = 1.
    do k = 1,5
    a(k+1) = a(k)/100.
    end do !j
    proba(j:6) = a/sum(a(j:6))
    call random_number(x1)
    tmp = 0
    do k = j,6
    tmp = tmp+proba(j)
    if(tmp>x1)exit
    end do   !j
    i = k
    end subroutine getfluxtype
!******************************************************
    end module module_growth