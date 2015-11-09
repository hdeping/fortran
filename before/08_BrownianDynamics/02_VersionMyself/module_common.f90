module module_common
    implicit none
    integer,parameter            :: n       = 100 ! n = nnum**2
    integer,parameter            :: nnum    = 10
    integer,parameter            :: m       = 2 ! denotes dimension
    real(8),parameter            :: limitBig= 1E10
    real(8),parameter            :: sigma   = 0.1  ! variance
    real(8),parameter            :: length(m) = (/11.0,11.0/)
    real(8),parameter            :: veloMax = 6.0 ! max velocity
    real(8),parameter            :: pi      = 3.1415927
    integer                      :: i
    integer                      :: j
    integer                      :: k
    integer                      :: partner(n) 
    real(8)                      :: coordinate(n,m)
    real(8)                      :: velocity(n,m)
    real(8)                      :: diameter(n)
    real(8)                      :: mass(n)
    ! record the colision time for every sphere
    real(8)                      :: collideTime(n) 
    real(8)                      :: timeIJ
    character(20)                :: filename 

    contains
!subroutine getInitialParameter{{{
subroutine getInitialParameter()
    integer              :: ii
    integer              :: jj
    integer              :: inum
    integer              :: jnum
    integer              :: coorIj(m)
    real(8)              :: tmp(m)


    ! initialize the diameter and the mass
    ! of the spheres

    mass     = 1.0
    diameter = 1.0


    ! initialize the velocity
    tmp = 0.0
    do ii = 1,m
        do jj = 1,n/2
            velocity(jj,ii)       = gauss()
            velocity(jj + n/2,ii) = - velocity(jj,ii)
        enddo !cycle ends
        tmp(ii) = sum(velocity(:,ii))
        print *,"the summation of velocity is ",tmp(ii)
    enddo !cycle ends
    ! initialize the coordinate
    do ii = 1,n
        inum   = mod(ii - 1,nnum) 
        jnum   = (ii - 1)/nnum
        coorIj = (/inum,jnum/)
        do jj  = 1,m
            coordinate(ii,jj) = length(jj)/2.0/nnum*(2*coorIj(jj) + 1)
        enddo !cycle ends
    enddo !cycle ends

    ! test the coordinatedo 
    filename = "coor_test.txt"
    open(10,file = filename)
    do ii = 1,n
        write(10,*)coordinate(ii,:)
    enddo !cycle ends
    close(10) 
     
end subroutine getInitialParameter
!}}}
!subroutine evolution{{{
subroutine evolution()
    integer,parameter      :: evoTimes = 100
    integer                :: ii 
    integer                :: jj
    integer                :: ij
    real(8)                :: timeLeast
    !integer                :: num
    integer,dimension(1)   :: it

    do ij =  1,evoTimes
     
        ! get the least time and its partner
        call getLeastTime()
        it(:) = minloc(collideTime)
        print *,"the minimum index is ",it(1)
        ii = it(1)
        jj = partner(ii)
        timeLeast = collideTime(ii)
        ! update the time and coordinate
        do ii = 1,n
            do jj = 1,m
                coordinate(ii,jj) = coordinate(ii,jj) + timeLeast*velocity(ii,jj)
                coordinate(ii,jj) = mod(coordinate(ii,jj),length(jj))
            enddo !cycle ends
             
        enddo !cycle ends
        ! update the velocity of the collide particles 
        call collideVelo(ii,jj)
   
    enddo !cycle ends


end subroutine evolution
!}}}
!!function gauss{{{
! get a gauss random number
!function gauss()
!    real(8)                    :: v
!    real(8)                    :: gauss
!    real(8)                    :: x1
!    real(8)                    :: x2
!    real(8)                    :: probability
!
!    x2          = 1.0
!    probability = 0.0
!    do while( x2 >  probability )
!        call random_number(x1) 
!        v = veloMax*(2.0*x1 - 1.0)
!        call random_number(x2)
!        probability = exp(- v**2.0/2.0/sigma)/sqrt(2.0*pi*sigma)
!    enddo !cycle ends
!    gauss = v
!end function gauss
!}}}
!function gauss{{{
!get a gauss random number
function gauss()
    real(8)                    :: v
    real(8)                    :: gauss
    real(8)                    :: x1
    integer                    :: ii


    v  = 0.0
    do ii = 1,6
        call random_number(x1)
        v = v + 2.0*x1 - 1.0
    enddo !cycle ends
    gauss = v
end function gauss
!}}}
!subroutine getLeastTime{{{
subroutine getLeastTime()
    integer                  :: ii
    integer                  :: jj
    real(8)                  :: tmp

    do ii = 1,n
        tmp = limitBig
        do jj = 1,n
            if(ii == jj)cycle
            if ( tmp > getCollisionTime(ii,jj) )then
                tmp = getCollisionTime(ii,jj)
                partner(ii) = jj
            endif ! if ends
        enddo !cycle ends
        collideTime(ii) = tmp
    enddo !cycle ends

end subroutine getLeastTime
!}}}
!function getCollisionTime{{{
function getCollisionTime(i,j)
    integer,intent(in)   :: i
    integer,intent(in)   :: j
    integer              :: ii
    !integer              :: jj
    real(8)              :: getCollisionTime
    real(8)              :: bij
    real(8)              :: rij
    real(8)              :: vij
    real(8)              :: relaVij(m)
    real(8)              :: relaRij(m)
    real(8)              :: a
    real(8)              :: deltaJudge

    ! get some variables
    a = (diameter(i) + diameter(j))/2.0
    bij = 0.0
    vij = 0.0
    rij = 0.0

    do ii = 1,m
        relaVij(ii) = velocity(i,ii) - velocity(j,ii)
        relaRij(ii) = coordinate(i,ii) - coordinate(j,ii)
    enddo !cycle ends

    do ii = 1,m
        vij = vij + relaVij(ii)**2.0
        rij = rij + relaRij(ii)**2.0
        bij = bij + relaRij(ii)*relaVij(ii)
    enddo !cycle ends

    deltaJudge = bij**2.0 - vij**2.0*(rij**2.0 - a**2.0)
    getCollisionTime = limitBig
    if ( bij < 0 )then
        if ( deltaJudge > 0 )then
            getCollisionTime = (- bij - sqrt(deltaJudge))/vij**2.0
        endif ! if ends
    endif ! if ends

end function getCollisionTime
!}}}
!subroutine collideVelo{{{
subroutine collideVelo(ii,jj)
    integer,intent(in)        :: ii
    integer,intent(in)        :: jj
    integer                   :: ij
    !integer                   :: k
    real(8)                   :: centerCycle(m)
    real(8)                   :: vertiCycle(m)
    real(8)                   :: components(m,m)  ! coordinate components
    !real(8)                   :: tmp

    ! transfer the velocity of the vertical part 
    ! to the direction between the center of circle
    do ij = 1,m
        centerCycle(ij) = coordinate(ii,ij) - coordinate(jj,ij)
        ! when m = 2, the vertiCycle can be expressed below
        vertiCycle(ij)  = centerCycle(m+1-ij)*(2*mod(ij,2) - 1)
    enddo !cycle ends
    
    ! get four coordinate components
    ! notice the order
    components(1,1) = dotProduct(velocity(jj,:),centerCycle,m)
    components(1,2) = dotProduct(velocity(ii,:),vertiCycle,m)
    components(2,1) = dotProduct(velocity(ii,:),centerCycle,m)
    components(2,2) = dotProduct(velocity(jj,:),vertiCycle,m)
    print *,dotProduct(centerCycle,vertiCycle,m)
    pause
    
    ! get new velocity
!    velocity(ii,:) = 0.0
!    velocity(jj,:) = 0.0

    print *, velocity(ii,:)
    print *, velocity(jj,:)
    do ij = 1,m
        velocity(ii,ij) = components(1,1)*centerCycle(ij) + &
                          components(1,2)*vertiCycle(ij)
        velocity(jj,ij) = components(2,1)*centerCycle(ij) + &
                          components(2,2)*vertiCycle(ij)
    enddo !cycle ends
    print *, velocity(ii,:)
    print *, velocity(jj,:)
     
   ! do ij = 1,m
   !     do k = 1,m
   !         velocity(ii,ij) = velocity(ii,ij) + centerCycle
   !     enddo !cycle ends
   ! enddo !cycle ends
     

end subroutine collideVelo
!}}}
!function dotProduct{{{
function dotProduct(a,b,n)
    integer,intent(in)          :: n
    real(8),intent(in)          :: a(n)
    real(8),intent(in)          :: b(n)
    real(8)                     :: dotProduct
    integer                     :: ii
    
    dotProduct = 0.0
    do ii = 1,n
        dotProduct = dotProduct + a(ii)*b(ii)
    enddo !cycle ends
end function dotProduct
!}}}
end module module_common
