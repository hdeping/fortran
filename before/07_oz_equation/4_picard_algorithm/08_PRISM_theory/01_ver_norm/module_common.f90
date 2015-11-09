!   m:  matrice
!   f:  fluid 
module module_common

    implicit none 
!parameters{{{
    integer,parameter           ::   l    = 12
    integer,parameter           ::   n    = 2**l
    integer,parameter           ::  fre   =  int(1E2)
    real(8),parameter           ::   pi   = 3.141592653 
    real(8),parameter           :: deltar = 0.04 
    real(8),parameter           :: deltak = pi/dble(n)/deltar 
    real(8),parameter           :: error  = 1E-8      !  for the diwwerences
    ! the diameter
    real(8),parameter           :: dpp    = 1.0       !  
    real(8),parameter           :: dww    = 1.0      !  
    real(8),parameter           :: dwp    = (dpp + dww)/2.0  
    real(8),parameter           :: dnp    = 1.0       !  
    real(8),parameter           :: dnw    = 1.0       !  
    ! the number dendity
    real(8),parameter           :: rhop   = 0.4
    real(8),parameter           :: rhow   = 0.4
    integer,parameter           :: ntmp   = 3
    ! prism
    real(8),parameter           ::  b0    = 1.269*dpp    ! length of bond
    integer,parameter           ::  nnum  = 50 ! number of polymer monomer 
!}}}
!common variables{{{
    real(8)                     ::  k
    real(8)                     ::  r
    real(8)                     ::  dk(n)     ! k
    real(8)                     ::  dr(n)     ! r
    real(8)                     ::  t1
    real(8)                     ::  t2
    real(8)                     ::  eta
    real(8)                     ::  xtmp
    real(8)                     ::  a(ntmp,ntmp)
    real(8)                     ::  b(ntmp)
    real(8)                     ::  bnew(ntmp)
    real(8)                     ::  x(ntmp)
    real(8)                     ::  x1
    real                        ::  rate
    ! judge the convergence
    real(8)                     ::  lambda_test
    real(8)                     ::  lambda    
    real(8)                     ::  lambda_final
    real(8)                     ::  lambda1  
    real(8)                     ::  lambda2  
    real(8)                     ::  test(n)   
    real(8)                     ::  test1(n)   
    real(8)                     ::  test2(n)   
    ! mayer function  
    real(8)                     ::  maypp(n)  ! p-p
    real(8)                     ::  maywp(n) ! w-p
    real(8)                     ::  mayww(n)  ! w-w
    real(8)                     ::  maynp(n) ! w-p
    real(8)                     ::  maynw(n)  ! w-w
    integer                     ::  i        
    integer                     ::  j        
    integer                     ::  times        
    character(20)               ::  filename 
!}}}
    ! variables for OZ equation  
!OZ{{{
    ! for convenience, h for H, c for C
    ! H = r*h, C = r*c
    ! p for polymer
    ! w for water
    ! n for nano particle
    !  k space 
    real(8)                     ::  hkpp(n)      
    real(8)                     ::  hkwp(n)      
    real(8)                     ::  hkww(n)      
    real(8)                     ::  hknw(n)      
    real(8)                     ::  hknp(n)      

    real(8)                     ::  gamma_kpp(n)      
    real(8)                     ::  gamma_kwp(n)      
    real(8)                     ::  gamma_kww(n)      
    real(8)                     ::  gamma_knw(n)      
    real(8)                     ::  gamma_knp(n)      

    real(8)                     ::  ckpp(n)      
    real(8)                     ::  ckwp(n)      
    real(8)                     ::  ckww(n)      
    real(8)                     ::  cknw(n)      
    real(8)                     ::  cknp(n)      
    !  r space                       
    real(8)                     ::  hrpp(n)      
    real(8)                     ::  hrwp(n)      
    real(8)                     ::  hrww(n)      
    real(8)                     ::  hrnw(n)      
    real(8)                     ::  hrnp(n)      

    real(8)                     ::  gamma_rpp(n)      
    real(8)                     ::  gamma_rwp(n)      
    real(8)                     ::  gamma_rww(n)      
    real(8)                     ::  gamma_rnw(n)      
    real(8)                     ::  gamma_rnp(n)      

    real(8)                     ::  crpp(n)      
    real(8)                     ::  crwp(n)      
    real(8)                     ::  crww(n)      
    real(8)                     ::  crnw(n)      
    real(8)                     ::  crnp(n)      
    ! g(r)
    real(8)                     ::  grpp(n)      
    real(8)                     ::  grwp(n)      
    real(8)                     ::  grww(n)      
    real(8)                     ::  grnw(n)      
    real(8)                     ::  grnp(n)      
    ! structure factor
    real(8)                     ::  skpp(n)      
    real(8)                     ::  skwp(n)      
    real(8)                     ::  skww(n)      
    real(8)                     ::  sknw(n)      
    real(8)                     ::  sknp(n)      
!}}}
    ! variables of prism
    real(8)                     ::  wp(n) ! partition function



end module module_common
