module module_common
    implicit none
    integer,parameter        :: n     = int(1E6)
    integer,parameter        :: fre   = 10
    real(8),parameter        :: pi    = 3.141592653
    real(8),parameter        :: error = 1E-4
    real(8)                  :: x(2)
    real(8)                  :: a    
    real(8)                  :: b    
    real(8)                  :: x1  ! for random number
    real(8)                  :: x2  ! for random number
    real(8)                  :: t1  
    real(8)                  :: t2  
    real(8)                  :: tmp
    real(8)                  :: lambda

    contains
!function getzero{{{
function getzero(a,b)
    real(8),intent(in)       :: a
    real(8),intent(in)       :: b
    real(8)                  :: getzero(2)
    real(8)                  :: atmp(2,2)
    real(8)                  :: btmp(2)
    real(8)                  :: xa(2)
    real(8)                  :: xb(2)
    integer                  :: k
    integer                  :: times
    
    xa = (/a,b/)
    times = 0
    do 
        ! get atmp
        atmp(1,1) = diffa(xa(1),xa(2))
        atmp(1,2) = diffb(xa(1),xa(2))
        atmp(2,1) = - atmp(1,2)
        atmp(1,1) = atmp(2,2)
        ! get btmp
        btmp(1)   = fab(xa(1),xa(2))
        btmp(2)   = gab(xa(1),xa(2))
        ! get xb
        xb = sol_equ(atmp,btmp,2)
        ! judge convergence 
        lambda = 0.0
        do k = 1,2
            lambda = lambda + abs(xb(k))
        end do
        times = times + 1
        !print *,times,"lambda = ",lambda
        ! get a new xa
        if(lambda < error)exit
        do k = 1,2
            xa(k) = xa(k) - 0.5*xb(k)
        end do
        if(lambda > 1000)exit
    end do
    getzero = xa
    print *,"times = ",times,"lambda = ",lambda

end function getzero
!}}}
!function fab{{{
function fab(a,b)
    real(8),intent(in)       :: a
    real(8),intent(in)       :: b
    real(8)                  :: fab
    integer                  :: k

    fab  = 0.0
    do k = 1,n
        fab = fab + exp(- a*log(dble(k)))&
              *cos(b*log(dble(k)))
    end do
end function fab
!}}}
!function gab{{{
function gab(a,b)
    real(8),intent(in)       :: a
    real(8),intent(in)       :: b
    real(8)                  :: gab
    integer                  :: k

    gab  = 0.0
    do k = 1,n
        gab = gab - exp(- a*log(dble(k)))&
              *sin(b*log(dble(k)))
    end do
end function gab
!}}}
!function diffa{{{
function diffa(a,b)
    real(8),intent(in)       :: a
    real(8),intent(in)       :: b
    real(8)                  :: diffa
    integer                  :: k

    diffa  = 0.0
    do k = 1,n
        diffa = diffa - log(dble(k))*exp(- a*log(dble(k)))&
              *cos(b*log(dble(k)))
    end do
end function diffa
!}}}
!function diffb{{{
function diffb(a,b)
    real(8),intent(in)       :: a
    real(8),intent(in)       :: b
    real(8)                  :: diffb
    integer                  :: k

    diffb  = 0.0
    do k = 1,n
        diffb = diffb - log(dble(k))*exp(- a*log(dble(k)))&
              *sin(b*log(dble(k)))
    end do
end function diffb
!}}}
!function sol_equ{{{
function sol_equ(a,b,n)
    ! in and out dummies
    integer,intent(in)    :: n
    real(8),intent(in)    :: a(n,n)
    real(8),intent(in)    :: b(n)
    real(8)               :: sol_equ(n)

    real(8)               :: c(n,n+1)
    real(8)               :: tmp
    real(8)               :: rate
    integer               :: s
    integer               :: i
    integer               :: j
    integer               :: k
    s = 0
    ! initial the new array c(n,n+1)
    do i = 1,n
        do j = 1,n
            c(i,j) = a(i,j)
        end do
        c(i,n+1) = b(i)
    end do
    !do i = 1,n
    !    print "(<n+1>f8.3)",c(i,:)
    !end do
    
     do  i = 1,n
         if(c(i,i) == 0)then
            do j = i+1,n
                if(c(j,i) /= 0)exit
            end do     !  i
            if(j == n+1)then
                s = 1
                cycle
            else
                do k = i,n+1          !  the ith column to  (n+1)th  column
                      tmp = c(i,k)
                   c(i,k) = c(j,k)
                   c(j,k) = tmp
               end do  !k
            endif
         endif
         if(s == 0)then
            do j = 1,n
                if(j == i)cycle
                if(c(j,i) == 0)cycle
                rate = c(j,i)/c(i,i)
                do k = i,n+1         !  the ith column to  (n+1)th  column
                    c(j,k) = c(j,k)-rate*c(i,k)
                end do   !k
            end do   !j
         end if     
     end do   !i
     do i = 1,n
         sol_equ(i) = c(i,n+1)/c(i,i)
     end do  !i
end function sol_equ
!}}}

end module module_common
