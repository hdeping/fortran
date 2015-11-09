module common_module

! parameters for lattice
integer,parameter        ::    nx_ltc = 60                ! number of grids in x direction 
integer,parameter        ::    ny_ltc = 60                ! number of grids in y direction
integer,parameter        ::    nx_barrier = 20 !20            ! the spatial period for attachment barrier in x direction, 
integer,parameter        ::    ny_barrier = 20 !20            ! the spatial period for attachment barrier in y direction, 

! variables for lattice
integer                    ::    n_neighbor
integer                    ::    status(nx_ltc,ny_ltc)    ! the present status of the lacttice
integer                    ::    growthstatus(nx_ltc,ny_ltc)    ! the present status of the lacttice
integer                    ::    class_site                ! the class of site (i,j).  3: body; 2: surface; 1: active; 0: nonactive
integer                    ::    type_site                ! the type of site (i,j) for attachment. 1: hollow 2: top 
integer                    ::    Neighbors(3,2)            ! the neighbors of site (i,j)
integer                    ::    n_cfgH(6)                ! number of hollow config can be attached by ci
integer                    ::    cfgH(6,ny_ltc,2)        ! hollow config can be attached by ci
integer                    ::    n_cfgT(6)                ! number of top config can be attached by ci
integer                    ::    cfgT(6,ny_ltc,2)        ! top config can be attached by ci
integer                    ::    n_cfgH_D(6)                ! number of hollow config can be detached via ci
integer                    ::    cfgH_D(6,ny_ltc,2)        ! hollow config can be detached via ci
integer                    ::    n_cfgT_D(6)                ! number of top config can be detached via ci
integer                    ::    cfgT_D(6,ny_ltc,2)        ! top config can be detached via ci

! variables for the present observation window
integer                    ::    xPos_ltc                ! the start position of the present observation window (in x direction)
integer                    ::    pos_front                ! the largest x position of the present graphene configure
integer                    ::    pos_graphene            ! the largest x position of the fullfilled graphene column

! variables for KMC simulation
integer                    ::    nT = 0                    ! the present number of MC steps
integer                    ::    nT_stay = 0
real                    ::    dt = .0
real*8                    ::    time = .0                    ! the present time
real                    ::    CE_time = 1E3
real                    ::    c_ci(6)                    ! concentration of ci clusters in equilibrium region, unit: grid^-1
real                    ::    flux_ciH(6)                ! net flux of ci clusters per hollow site in growth region , unit: grid^-1
real                    ::    flux_ciT(6)                ! net flux of ci clusters per top site in growth region , unit: grid^-1
real                    ::    k_ciH(6)                ! the relative reaction probability of each process
real                    ::    k_ciT(6)                ! the relative reaction probability of each process
real                    ::    k_ciH_D(6)                ! the relative reaction probability of each process
real                    ::    k_ciT_D(6)                ! the relative reaction probability of each process
real                    ::    cnt_growth(6)
! parameters for KMC events
real                    ::    dEiH(6) = (/-0.86,-1.37,-2.57,-3.43,-4.06,-4.37/)
                                                    ! potential energy change for attachment of Ci on a hollow site
real                    ::    dEiH_B(6) = (/1.0,1.0,2.03,0.08,0.08,0.08/)    !(/1.0,1.0,2.03,1.51,1.51,1.51/)
real                    ::    dEiH_B2(6) = (/1.86,2.36,4.59,3.49,4.12,4.12/)    !(/1.0,1.0,2.03,1.51,1.51,1.51/)
real                    ::    dEiT(6) = (/2.13,1.44,0.12,-1.44,-4.06,-4.37/)    !(/2.13,1.25,0.12,-1.44,-1.85,-2.16/)
                                                    ! potential energy change for attachment of Ci on a top site
real                    ::    dEiT_B(6) = (/2.30,1.61,2.03,0.08,0.08,0.08/)    !(/2.30,1.61,2.03,0.08,0.08,0.08/)
real                    ::    dEiT_B2(6) = (/0.17,0.17,1.91,3.49,4.12,4.42/)    !(/0.17,0.17,1.91,1.52,4.12,4.42/)
real                    ::    D_Ci(6) = (/0.71,0.73,0.37,0.54,0.4,0.54/) !(/0.71,0.76,0.93,0.54,0.4,0.54/)        ! diffusion barrier for Ci
real,parameter            ::    T = 1170                    ! unit: K
real                    ::    beta = 1.6022/(1.3807E-4*T)! the inverse temperature, unit: eV^-1
real                    ::    timeScale(2) = (/2.E1,1.E-2/)

! parameters for drawing graphic
integer(2),parameter        ::    x_scr = 10                ! location of start point for drawing graphic, Unit: pixel
integer(2),parameter        ::    y_scr = 10
integer(2),parameter        ::    w_scr = 1200                ! width of the graphic, Unit: pixel
integer(2),parameter        ::    h_scr = 800                ! heigth of the graphic, Unit: pixel
integer                            ::    hh_grid                    ! half heigth of a grid on screen, Unit: pixel
integer                            ::    len_grid                ! length of a grid on screen, Unit: pixel
integer,parameter        ::    color_body = 1 !3                
integer,parameter        ::    color_surface = 1 !6                
integer,parameter        ::    color_active = 1 !2                    
integer,parameter        ::    color_inactive = 1 !2                
integer,parameter        ::    color_bond = 1 !7
integer,parameter        ::    color_background = 2 !1
integer,parameter        ::    color_spc = 4
integer,parameter        ::    showBody = 1
integer,parameter        ::    showSurface = 1
integer,parameter        ::    showActive = 1
integer,parameter        ::    showInactive = 1
integer,parameter        ::    showSPC = 1

!variables for drawing graphic
integer                    ::    f_redraw = 0

! variables for temporary use
integer                    ::    newInitial = 1
real                    ::    tmp
integer                    ::    P_count(6) = 0
character*50            ::    fName

contains
!******************************************************
subroutine initial()
integer :: i_grph = 2                                    ! the initial graphene width
hh_grid = int2(h_scr/ny_ltc)-1            
len_grid = int2(2.*hh_grid/3.**(0.5))+8
call getReactionProb()
!do ii = 1,6
!    prob_ci(ii) = flux_ci(ii)/sum(flux_ci)
!enddo !ii
! the initial graphene configure
status = 0
growthStatus = 0
if(newInitial = =1) then
    status(1:i_grph,:) = 1
    pos_graphene = i_grph
    pos_front = i_grph
!    do i = i_grph+1,nx_ltc
!        do j = 1,ny_ltc
!            call random_number(tmp)
!            status(i,j) = floor(tmp/0.8)
!        enddo !j
!    enddo !i
else
    open(11,file = 'initial.dat')
    do ii = 1,nx_ltc
        read(11,"(60I5)") status(ii,:)
    enddo !ii
    close(11)    
    xPos_ltc = 0
    pos_graphene = 0;pos_front=0
    do ii = 1,nx_ltc
        if(sum(status(ii,:)) = =ny_ltc.and.pos_graphene+1==ii) pos_graphene=ii
        do jj = 1,ny_ltc
            if(status(ii,jj) = =1.and.ii>pos_front) pos_front=ii
        enddo !jj
    enddo !ii 
endif
end subroutine initial
!******************************************************
subroutine getEqDis()
real,parameter            ::    pE_ci(6) = (/-7.43,-7.12,-7.07,-7.05,-7.14,-7.2/)    ! pE_ci(6)=(/-7.43,-7.12,-7.07,-7.05,-7.14,-7.2/)                
                                                    ! the potential energy of ci 
open(10,file = 'ci.dat')
write(10,"(G18.9)") c_ci(1)
do i = 2,6
    c_ci(i) = ((c_ci(1))**i)*exp(i*(pE_ci(1)-pE_ci(i))*beta)
    write(10,"(G18.9)") c_ci(i)
enddo !i
close(10)
end subroutine getEqDis
!******************************************************
subroutine getFlux()
call getEqDis()
open(10,file = 'flux_i.dat')
flux_ciH = (c_ci)*exp(-D_Ci*beta)/sum(exp(-D_Ci*beta)) !*timescale(1)    !*1.E13*CE_Diff    
flux_ciT = (c_ci)*exp(-D_Ci*beta)/sum(exp(-D_Ci*beta)) !*timescale(1)    !*1.E13*CE_Diff    
do i = 1,6
    write(10,*) i,flux_ciH(i),flux_ciT(i)
enddo !i
close(10)
end subroutine getFlux
!******************************************************
subroutine getReactionProb()
call getFlux()
open(10,file = 'k.dat')
k_ciH = flux_ciH*1.E13*exp(-dEiH_B*beta) !* timeScale(1)
k_ciH_D = exp(-dEiH_B2*beta)*1.E13    
k_ciT = flux_ciT*1.E13*exp(-dEiT_B*beta) !* timeScale(1)
k_ciT_D = exp(-dEiT_B2*beta)*1.E13
do i = 1,6
    write(10,"(I5,4G18.9)") i,k_ciH(i),k_ciH_D(i),k_ciT(i),k_ciT_D(i)
!                k_ciH(i)*timescale(2),k_ciH_D(i),k_ciT(i)*timescale(2),k_ciT_D(i)
enddo !i
close(10)
end subroutine getReactionProb
!******************************************************
end module common_module
