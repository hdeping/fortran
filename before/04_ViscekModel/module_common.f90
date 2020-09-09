module module_common
    implicit none
    integer,parameter            :: n   = 100000      ! times of evolution
    integer,parameter            :: num = 40          ! number of particles
    integer,parameter            :: fre = 100         ! frequence of printing
    real(8),parameter            :: pi  = 3.141592653 
    real(8),parameter            :: l   = 7D0         ! size of square
    real(8),parameter            :: v   = 3D-2        ! velocity
    real(8),parameter            :: r   = 1D0         ! interation radius
    real(8)                      :: eta               ! noise
    real(8)                      :: coor(num,2)       ! position
    real(8)                      :: angle(num)        ! noise
    real(8)                      :: vMean             ! mean velocity
    real(8)                      :: angleMean(num)    ! mean velocity
    real(8)                      :: t1
    real(8)                      :: t2

    integer                      :: i
    integer                      :: j
    integer                      :: k
    character(10)                :: filename 
    contains
!subroutine init_config{{{
! initialize the configuration
! position and the direction 
! of the velocities
subroutine init_config()
    integer                      :: ii
    integer                      :: jj
    real(8)                      :: x1
    ! initial the position
    do ii = 1,num
        do jj = 1,2
            call random_number(x1) 
            coor(ii,jj) = l*x1
        enddo !cycle ends
    enddo !cycle ends
    ! initialize the directions (0,2*pi)
    do ii = 1,num
        call random_number(x1)
        angle(ii) = 2*pi*x1
    enddo !cycle ends
end subroutine init_config
!}}}
!subroutine update{{{
subroutine update()
    integer                      :: ii
    integer                      :: jj
    real(8)                      :: delta_theta
    real(8)                      :: x1
    real(8)                      :: meanAngle(num)

    ! update the position
    do ii = 1,num
        coor(ii,1) = coor(ii,1) + v*cos(angle(ii))
        coor(ii,2) = coor(ii,2) + v*sin(angle(ii))
    enddo !cycle ends
    coor(:,1) = ch_range(coor(:,1),num,l)
    coor(:,2) = ch_range(coor(:,2),num,l)
    ! update the direction
    ! get mean angle
    meanAngle = getmean()
    do ii = 1,num
        call random_number(x1)
        delta_theta = eta*(x1 - 0.5)
        angle(ii)   = meanAngle(ii) + delta_theta
    enddo !cycle ends
    angle = ch_range(angle,num,2.0*pi)
end subroutine update
!}}}
!function getmean{{{
function getmean()
    integer                      :: ii
    integer                      :: jj
    integer                      :: kk
    integer                      :: mark      ! mark all the condition
    integer                      :: statu(num)! statu of particle
    integer                      :: times(num)! surrounding number with r 
    real(8)                      :: distance
    real(8)                      :: getmean(num)

    ! get statu of particles
    do ii = 1,num
       if ( coor(ii,1) < r )then
           jj = 1
       elseif ( coor(ii,1) < l - r )then
           jj = 2
       else
           jj = 3
       endif ! if ends
       if ( coor(ii,2) < r )then
           kk = 0
       elseif ( coor(ii,2) < l - r )then
           kk = 1
       else
           kk = 2
       endif ! if ends
       statu(ii) = 3*kk + jj
    enddo !cycle ends
    ! initialize the getmean and times
    getmean = angle
    times   = 1
    ! get the distance between two particles 
    ! get the meanAngle
!{{{
    do ii = 1,num - 1
        do jj = ii + 1,num 
            ! sort the condition 
            if     ( statu(ii) == 1 .and. statu(jj) == 7)then
                mark = 1
            elseif ( statu(ii) == 2 .and. statu(jj) == 8)then
                mark = 2
            elseif ( statu(ii) == 3 .and. statu(jj) == 9)then
                mark = 3
            elseif ( statu(ii) == 7 .and. statu(jj) == 1)then
                mark = 4
            elseif ( statu(ii) == 8 .and. statu(jj) == 2)then
                mark = 5
            elseif ( statu(ii) == 9 .and. statu(jj) == 3)then
                mark = 6
            elseif ( statu(ii) == 7 .and. statu(jj) == 9)then
                mark = 7
            elseif ( statu(ii) == 4 .and. statu(jj) == 6)then
                mark = 8
            elseif ( statu(ii) == 1 .and. statu(jj) == 3)then
                mark = 9
            elseif ( statu(ii) == 9 .and. statu(jj) == 7)then
                mark = 10
            elseif ( statu(ii) == 6 .and. statu(jj) == 4)then
                mark = 11
            elseif ( statu(ii) == 3 .and. statu(jj) == 1)then
                mark = 12
            elseif ( statu(ii) == 3 .and. statu(jj) == 7)then
                mark = 13
            elseif ( statu(ii) == 7 .and. statu(jj) == 3)then
                mark = 14
            elseif ( statu(ii) == 1 .and. statu(jj) == 9)then
                mark = 15
            elseif ( statu(ii) == 9 .and. statu(jj) == 1)then
                mark = 16
            else
                mark = 17
            endif ! if ends
            ! get distances
            select case(mark)
            case(1,2,3)
                distance = sqrt((coor(ii,1) - coor(jj,1))**2.0&
                                + (coor(ii,2) - coor(jj,2) + l)**2.0)
            case(4,5,6)
                distance = sqrt((coor(ii,1) - coor(jj,1))**2.0&
                                + (coor(ii,2) - coor(jj,2) - l)**2.0)
            case(7,8,9)
                distance = sqrt((coor(ii,1) - coor(jj,1))**2.0&
                                + (coor(ii,2) - coor(jj,2) + l)**2.0)
            case(10,11,12)
                distance = sqrt((coor(ii,1) - coor(jj,1) + l)**2.0&
                                + (coor(ii,2) - coor(jj,2) )**2.0)
            case(13)
                distance = sqrt((coor(ii,1) - coor(jj,1) - l)**2.0&
                                + (coor(ii,2) - coor(jj,2) + l )**2.0)
            case(14)
                distance = sqrt((coor(ii,1) - coor(jj,1) + l)**2.0&
                                + (coor(ii,2) - coor(jj,2) - l)**2.0)
            case(15)
                distance = sqrt((coor(ii,1) - coor(jj,1) + l)**2.0&
                                + (coor(ii,2) - coor(jj,2) + l)**2.0)
            case(16)
                distance = sqrt((coor(ii,1) - coor(jj,1) - l)**2.0&
                                + (coor(ii,2) - coor(jj,2) - l)**2.0)
            case default
                distance = sqrt((coor(ii,1) - coor(jj,1))**2.0&
                                + (coor(ii,2) - coor(jj,2))**2.0)
            end select ! select ends
            
            ! get the closed particles within r (r = 1 here)
            ! of every particle
            if ( distance < r )then
                times(ii)   = times(ii) + 1
                times(jj)   = times(jj) + 1
                getmean(ii) = (getmean(ii)*(times(ii) - 1) + &
                              angle(jj))/dble(times(ii))
                getmean(jj) = (getmean(jj)*(times(jj) - 1) + &
                              angle(ii))/dble(times(jj))
                !print *,ii,jj,"distance = ",distance
                !print "(4f12.4)",coor(ii,1),coor(ii,2),coor(jj,1),coor(jj,2)
                !print *,statu(ii),statu(jj),mark
                !pause
            endif ! if ends

        enddo !cycle ends
    enddo !cycle ends
!}}}
     
end function getmean
!}}}
!function ch_range{{{
function ch_range(a,n,inter)
    integer,intent(in)           :: n
    real(8),intent(in)           :: inter
    real(8),intent(in)           :: a(n)
    integer                      :: ii
    real(8)                      :: ch_range(n)
    ch_range = a

    do ii = 1,n
        if ( ch_range(ii) < 0 )then
            do  while(ch_range(ii) < 0)
                ch_range(ii) = ch_range(ii) + inter
            enddo !cycle ends
        elseif ( ch_range(ii) > inter )then
            do  while(ch_range(ii) > inter)
                ch_range(ii) = ch_range(ii) - inter
            enddo !cycle ends
        endif ! if ends
    enddo !cycle ends
end function ch_range
!}}}
!function getv{{{
! get the mean velocity
! a for angle
! n for the number of particles
function getv(a,n)
    integer,intent(in)           :: n
    real(8),intent(in)           :: a(n)
    integer                      :: ii
    real(8)                      :: getv
    real(8)                      :: xtmp
    real(8)                      :: ytmp

    xtmp  = 0D0
    ytmp  = 0D0
    do ii = 1,n
        xtmp = xtmp + cos(a(ii))
        ytmp = ytmp + sin(a(ii))
    enddo !cycle ends
    getv = sqrt(xtmp**2.0 + ytmp**2.0)/dble(n) 
end function getv
!}}}
end module module_common
