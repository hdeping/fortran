!   m:  matrice
!   f:  fluid 
module module_common

    implicit none 
!parameters{{{
    integer,parameter           ::   l    = 10
    integer,parameter           ::   n    = 2**l
    integer,parameter           ::  ncut  = 300
    ! cut-off of n
    integer,parameter           :: tmnum  = 200
    integer,parameter           ::   m    = 2 !number components
    integer,parameter           ::  fre   =  int(1E3)
    integer,parameter           ::  nnum  = 40
    real(8),parameter           ::   pi   = 3.141592653 
    real(8),parameter           :: deltar = 0.01
    real(8),parameter           ::    h   = pi/dble(n)/deltar ! deltak
    real(8),parameter           :: error  = 1E-10      !  for the diwwerences
    ! the diameter
    real(8),parameter           :: d11    = 1.0       !  
    real(8),parameter           :: d22    = 1.0      !  
    real(8),parameter           :: d12    = (d11 + d22)/2.0  
    real(8),parameter           :: d(m*m) = (/d11,d12,d22,0D0/)
    ! the number dendity
    real(8),parameter           :: rho    = 0.9
    !   matrix T and dimensionless
    real(8),parameter           :: atom   = 1.6605402E-27
    real(8),parameter           :: temper = 300.0           ! temperature 
    real(8),parameter           :: kb     = 1.3806488E-23
    real(8),parameter           :: tpture = kb*temper/atom
    real(8),parameter           :: mass1  = 18.0            ! water
    real(8),parameter           :: mass2  = 28.0*dble(nnum)  ! polymer
    real(8),parameter           :: v      = 1.0
    real(8),parameter           :: gold   = (sqrt(5.0) - 1.0)/2.0
!}}}
!common variables{{{
    real(8)                     ::  dt     ! step length
    real(8)                     ::  l1    ! in MCT
    real(8)                     ::  l2    ! in MCT
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
    real(8)                     ::  lambda2
    real(8)                     ::  test(m,m,n)   
    real(8)                     ::  test1(n)   
    real(8)                     ::  test2(n)   
    ! mayer function  
    real(8)                     ::  mayfun(n)  ! p-p
    integer                     ::  i
    integer                     ::  j
    integer                     ::  t
    integer                     ::  q  ! summation    
    integer                     ::  p  ! summation      
    integer                     ::  k  ! summation 
    integer                     ::  times        
    character(20)               ::  filename 
!}}}
    ! variables for OZ equation  
!OZ{{{
    ! for convenience, h for H, c for C
    !  k space 
    real(8)                     ::  hk(ncut)      
    real(8)                     ::  ck(ncut)      
    real(8)                     ::  gamma_k(ncut)      
    !  r space                       
    real(8)                     ::  hr(ncut)      
    real(8)                     ::  cr(ncut)      
    real(8)                     ::  gamma_r(ncut)      
    ! g(r)
    real(8)                     ::  gr(ncut)      
    ! structure factor
    real(8)                     ::  sk(ncut)      
!}}}
!  MCT equation
    ! self-scattering
    real(8)                     ::  f(ncut,tmnum) 
    real(8)                     ::  finalf(ncut) 
    ! self-scattering partial time
    real(8)                     ::  par_f(ncut,tmnum) 
    real(8)                     ::  memory(ncut,tmnum)! memory kernel
    real(8)                     ::  diffu      ! diffusion
    real(8)                     ::  mat_A(ncut,ncut,ncut)
    !real(8)                     ::  mat_B(ncut)
    !real(8)                     ::  mat_D(ncut)
    real(8)                     ::  mat_R(ncut)
    real(8)                     ::  mat_K(ncut)
    real(8)                     ::  mat_U(ncut)
    real(8)                     ::  mat_V(ncut)
    real(8)                     ::  inver_sk(ncut)



end module module_common
