!   m:  matrice
!   f:  fluid 
module module_common

    implicit none 
    integer,parameter           ::  lnum  = 10   ! cut-off of n
    integer,parameter           ::  n     = 2**lnum! cut-off of n
    integer,parameter           ::  ncut  = 200   ! cut-off of n
    integer,parameter           ::   m    = 2 !number components
    integer,parameter           ::  fre   =  int(1E3)
    real(8),parameter           ::   pi   = 3.141592653 
    real(8),parameter           :: deltar = 0.01
    real(8),parameter           ::    h   = pi/dble(n)/deltar ! deltak
    real(8),parameter           :: error  = 1E-8      !  for the diwwerences
    ! the number dendity
    real(8),parameter           :: rho    = 0.9
    ! judge the convergence
    real(8)                     ::  lambda
    integer                     ::  times
    integer                     ::  p
    integer                     ::  q
    integer                     ::  k
    integer                     ::  l
    character(20)               ::  filename 
    ! variables for OZ equation  
    real(8)                     ::  ck(ncut)      
    real(8)                     ::  sk(ncut)      
    real(8)                     ::  dk(n)      
    real(8)                     ::  deltafun(ncut,ncut)      
    !  MCT equation
    ! self-scattering
    real(8)                     ::  finalf(ncut) 
    real(8)                     ::  memory(ncut) 
    ! self-scattering partial time
    real(8)                     ::  diffu      ! diffusion
    real(8)                     ::  mat_A(ncut,ncut,ncut)
end module module_common
