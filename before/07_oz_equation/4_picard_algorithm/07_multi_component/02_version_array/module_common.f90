!   m:  matrice
!   f:  fluid 
module module_common

    implicit none 
!parameters{{{
    integer,parameter           ::   l    = 10
    integer,parameter           ::   n    = 2**l
    integer,parameter           ::   m    = 2
    integer,parameter           ::  fre   =  int(1E2)
    real(8),parameter           ::   pi   = 3.141592653 
    real(8),parameter           :: deltar = 0.01
    real(8),parameter           :: deltak = pi/dble(n)/deltar 
    real(8),parameter           :: error  = 1E-8      !  for the diwwerences
    ! the diameter
    real(8),parameter           :: d11    = 1.0       !  
    real(8),parameter           :: d22    = 1.0      !  
    real(8),parameter           :: d12    = (d11 + d22)/2.0  
    ! the number dendity
    real(8),parameter           :: rho1   = 0.01
    real(8),parameter           :: rho2   = 0.1
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
    real(8)                     ::  x1
    real(8)                     ::  rate
    ! judge the convergence
    real(8)                     ::  lambda
    real(8)                     ::  lambda1
    real(8)                     ::  test(m,m,n)   
    real(8)                     ::  test1(n)   
    real(8)                     ::  test2(n)   
    real(8)                     ::  test3(n)   
    ! mayer function  
    real(8)                     ::  mayfun(m,m,n)  ! p-p
    real(8)                     ::  may11(n)
    real(8)                     ::  may12(n)
    real(8)                     ::  may22(n)
    integer                     ::  i        
    integer                     ::  j        
    integer                     ::  times        
    character(20)               ::  filename 
!}}}
    ! variables for OZ equation  
!OZ{{{
    ! for convenience, h for H, c for C
    !  k space 
    real(8)                     ::  hk11(n)
    real(8)                     ::  hk12(n)
    real(8)                     ::  hk22(n)
    real(8)                     ::  ck11(n)
    real(8)                     ::  ck12(n)
    real(8)                     ::  ck22(n)
    real(8)                     ::  gamma_k11(n)
    real(8)                     ::  gamma_k12(n)
    real(8)                     ::  gamma_k22(n)
    !  r space                       
    real(8)                     ::  hr11(n)
    real(8)                     ::  hr12(n)
    real(8)                     ::  hr22(n)
    real(8)                     ::  cr11(n)
    real(8)                     ::  cr12(n)
    real(8)                     ::  cr22(n)
    real(8)                     ::  gamma_r11(n)
    real(8)                     ::  gamma_r12(n)
    real(8)                     ::  gamma_r22(n)
    ! g(r)
    real(8)                     ::  gr11(n)
    real(8)                     ::  gr12(n)
    real(8)                     ::  gr22(n)
    ! structure factor
    real(8)                     ::  sk11(n)
    real(8)                     ::  sk12(n)
    real(8)                     ::  sk22(n)
!}}}
    ! variables of prism



end module module_common
