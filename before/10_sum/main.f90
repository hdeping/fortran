program main
    implicit none
    integer  i
    integer  j
    integer  k
    integer  m
    integer  n
    integer  s
    
    do n= 1,10
        s=0
        do i=1,n
            do j=i+1,n
                do k=j+1,n
                    do m=k+1,n
                        s= s + i + j + k + m
                    end do  !  m
                end do  !  k
            end do  !  j 
        end do  !  i
        print *,n,s,fun(n+1,5),fun(n+2,5),fun(n,5)
    end do  !  n
    
    contains
    !****************************************************
    real function fun(n,m)
    integer  n,m
    integer  i
    real  tmp
    tmp=1
    do i=n-m+1,n
        tmp=tmp*i
    end do  !  i
    fun=tmp
    end function fun
    !****************************************************
end program main
