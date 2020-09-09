!   m:  matrice
!   f:  fluid 
module module_common

    implicit none 
!parameters{{{
    integer,parameter           ::   l    = 10
    integer,parameter           ::   n    = 2**l
    integer,parameter           ::   m    = 2 !number components
    integer,parameter           ::  fre   =  int(1E3)
    real(8),parameter           ::   pi   = 3.141592653 
    real(8),parameter           :: deltar = 0.01
    real(8),parameter           :: deltak = pi/dble(n)/deltar 
    real(8),parameter           :: error  = 1E-8      !  for the diwwerences
    ! the diameter
    real(8),parameter           :: d11    = 1.0       !  
    real(8),parameter           :: d22    = 1.0      !  
    real(8),parameter           :: d12    = (d11 + d22)/2.0  
    real(8),parameter           :: d(m*m) = (/d11,d12,d22,0D0/)
    ! the number dendity
    real(8),parameter           :: rho1   = 0.5
    real(8),parameter           :: rho2   = 0.3
    real(8),parameter           :: rho(m) = (rho1,rho2)
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
    real(8)                     ::  deltafun(m,m) ! delta function
    real(8)                     ::  a(m,m)
    real(8)                     ::  b(m,m)
    real(8)                     ::  bnew(m)
    real(8)                     ::  x(m)
    real(8)                     ::  x1
    real(8)                     ::  rate
    ! judge the convergence
    real(8)                     ::  lambda
    real(8)                     ::  lambda1
    real(8)                     ::  test(m,m,n)   
    real(8)                     ::  test1(n)   
    real(8)                     ::  test2(n)   
    ! mayer function  
    real(8)                     ::  mayfun(m,m,n)  ! p-p
    integer                     ::  i        
    integer                     ::  j        
    integer                     ::  times        
    character(20)               ::  filename 
!}}}
    ! variables for OZ equation  
!OZ{{{
    ! for convenience, h for H, c for C
    !  k space 
    real(8)                     ::  hk(m,m,n)      
    real(8)                     ::  ck(m,m,n)      
    real(8)                     ::  gamma_k(m,m,n)      
    !  r space                       
    real(8)                     ::  hr(m,m,n)      
    real(8)                     ::  cr(m,m,n)      
    real(8)                     ::  gamma_r(m,m,n)      
    ! g(r)
    real(8)                     ::  gr(m,m,n)      
    ! structure factor
    real(8)                     ::  sk(m,m,n)      
!}}}
    ! variables of prism



end module module_common
