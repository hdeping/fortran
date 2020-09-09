!   m:  matrice
!   f:  fluid 
module module_common

    implicit none 
!parameters{{{
    integer,parameter           ::   l    = 10
    integer,parameter           ::   n    = 2**l
    integer,parameter           ::  ncut  = 100   ! cut-off of n
    integer,parameter           :: tmnum  = 200
    integer,parameter           ::   m    = 2 !number components
    integer,parameter           ::  fre   =  int(1E3)
    integer,parameter           ::  nnum  = 40
    real(8),parameter           ::   pi   = 3.141592653 
    real(8),parameter           :: deltar = 0.01
    real(8),parameter           ::    h   = pi/dble(n)/deltar ! deltak
    real(8),parameter           :: error  = 1E-8      !  for the diwwerences
    ! the diameter
    real(8),parameter           :: d11    = 1.0       !  
    real(8),parameter           :: d22    = 1.0      !  
    real(8),parameter           :: d12    = (d11 + d22)/2.0  
    real(8),parameter           :: d(m*m) = (/d11,d12,d22,0D0/)
    ! the number dendity
    real(8),parameter           :: rho1   = 0.4
    real(8),parameter           :: rho2   = 0.4
    real(8),parameter           :: rhosum = rho1 + rho2
    real(8),parameter           :: xrate1 = rho1/rhosum
    real(8),parameter           :: xrate2 = 1 - xrate1
    real(8),parameter           ::xrate(m)= (/xrate1,xrate2/)
    real(8),parameter           :: rho(m) = (/rho1,rho2/)
    !   matrix T and dimensionless
    real(8),parameter           :: atom   = 1.6605402E-27
    real(8),parameter           :: temper = 300.0           ! temperature 
    real(8),parameter           :: kb     = 1.3806488E-23
    real(8),parameter           :: tpture = kb*temper/atom
    real(8),parameter           :: mass1  = 18.0            ! water
    real(8),parameter           :: mass2  = 28.0*dble(nnum)  ! polymer
    real(8),parameter           :: v1     = 1.0
    real(8),parameter           :: v2     = 1.0
    real(8),parameter           :: v(m)   = (/v1,v2/) 
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
    real(8)                     ::  mayfun(m,m,n)  ! p-p
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
    real(8)                     ::  hk(m,m,ncut)      
    real(8)                     ::  ck(m,m,ncut)      
    real(8)                     ::  gamma_k(m,m,ncut)      
    !  r space                       
    real(8)                     ::  hr(m,m,ncut)      
    real(8)                     ::  cr(m,m,ncut)      
    real(8)                     ::  gamma_r(m,m,ncut)      
    ! g(r)
    real(8)                     ::  gr(m,m,ncut)      
    ! structure factor
    real(8)                     ::  sk(m,m,ncut)      
!}}}
!  MCT equation
    ! self-scattering
    real(8)                     ::  f(m,m,ncut,tmnum) 
    ! self-scattering partial time
    real(8)                     ::  par_f(m,m,ncut,tmnum) 
    real(8)                     ::  memory(m,m,ncut,tmnum)! memory kernel
    real(8)                     ::  diffu(m,m)       ! diffusion
    real(8)                     ::  mat_A(m,m,ncut)
    real(8)                     ::  mat_B(m,m,ncut)
    real(8)                     ::  mat_D(m,m,ncut)
    real(8)                     ::  mat_R(m,m,ncut)
    real(8)                     ::  mat_K(m,m,ncut)
    real(8)                     ::  mat_U(m,m,ncut)
    real(8)                     ::  mat_V(m,m,ncut)
    real(8)                     ::  inver_sk(m,m,ncut)

    contains
!*********** judge ********************
!function judge{{{
! compare the difference between two arrays
! judge the convergence
function judge(a,b)
    real(8),intent(in)      :: a(n)
    real(8),intent(in)      :: b(n)
    real(8)                 :: judge
    integer                 :: ii
    
    judge = 0
    do ii = 1,n 
        judge = judge + abs(a(ii) - b(ii))
    end do
    !judge = judge/dble(n)
    
end function judge
!}}}
!function bisetion{{{
! compare the difference between two arrays
! judge the convergence
function bisetion(a,b)
    real(8),intent(in)      :: a(n)
    real(8),intent(inout)   :: b(n)
    real(8)                 :: bisetion
    integer                 :: ii
    
    bisetion = judge(a,b)
    do ii = 1,n 
        b(ii) = (a(ii) + b(ii))/2.0
    end do
    !judge = judge/dble(n)
    
end function bisetion
!}}}
!function setion_rate{{{
! compare the difference between two arrays
! judge the convergence
function setion_rate(a,b,rate)
    real(8),intent(in)      :: a(n)
    real(8),intent(inout)   :: b(n)
    real(8),intent(in)      :: rate
    real(8)                 :: setion_rate
    integer                 :: ii
    
    setion_rate = judge(a,b)
    do ii = 1,n 
        b(ii) = a(ii)*(1.0 - rate) + b(ii)*rate
    end do
    
end function setion_rate
!}}}
!function conver{{{
!  judge convergence with golden setion 
function conver(a,b)
    real(8),intent(in)      :: a(n)
    real(8),intent(inout)   :: b(n)
    real(8)                 :: conver
    integer                 :: ii

    do ii = 1,n
        test1(ii) = a(ii)*gold + b(ii)*(1.0 - gold)
        test2(ii) = a(ii)*(1.0 - gold) + b(ii)*gold 
    end do

    lambda1 = judge(a,test1)
    lambda2 = judge(a,test2)
    if(lambda1 < lambda2)then
        b      = test1
        conver = lambda1
    else
        b      = test2
        conver = lambda2
    endif
    
    
end function conver
!}}}


end module module_common
