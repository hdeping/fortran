!   m:  matrice
!   f:  fluid 
module module_common

    implicit none 
!*********   variables  ********************************
!variables{{{
!parameters{{{
    integer,parameter           ::   l    = 10
    integer,parameter           ::   n    = 2**l
    integer,parameter           ::  fre   =  int(1E3)
    real(8),parameter           ::   pi   = 3.141592653 
    !real(8),parameter           ::  top   = 10.24
    real(8),parameter           :: deltar = 0.01
    real(8),parameter           :: deltak = pi/deltar/dble(n)
    real(8),parameter           :: error  = 1E-8              !  for the differences
    ! diameters
    real(8),parameter           :: dmm    = 1.0                
    real(8),parameter           :: dff    = 1.0                
    real(8),parameter           :: dss    = 1.0
    !real(8),parameter           :: dfm    = (dmm + dff)/2.0   
    !real(8),parameter           :: dsm    = (dss + dmm)/2.0   
    !real(8),parameter           :: dsf    = (dss + dff)/2.0   
    real(8),parameter           :: dfm    = 1.0   
    real(8),parameter           :: dsm    = 1.0   
    real(8),parameter           :: dsf    = 1.0   

    real(8),parameter           :: rhom   = 0.5              !  the density of matrix
    real(8),parameter           :: rhof   = 0.5                !  the density of fluid 
    real(8),parameter           :: gold   = (sqrt(5.0) - 1.0)/2.0  ! golden rate
!}}}
!common variables{{{
    !  variables for fft
    integer                     :: status 
    !type(dfti_descriptor), pointer :: my_desc1_handle
    !type(dfti_descriptor), pointer :: my_desc2_handle
    ! variables for OZ equation  
    real(8)                     ::  k
    real(8)                     ::  r
    real(8)                     ::  t1
    real(8)                     ::  t2
    real(8)                     ::  eta
    real                        ::  rate
    ! mayer function
    real(8)                     ::  maymm(n)  
    real(8)                     ::  mayfm(n)  
    real(8)                     ::  mayff(n)  
    real(8)                     ::  maysm(n)  
    real(8)                     ::  maysf(n)  
    real(8)                     ::  mayss(n)  
    ! judge convergence
    real(8)                     ::  lambda_final
    real(8)                     ::  lambda    
    real(8)                     ::  lambda1  
    real(8)                     ::  lambda2  
    real(8)                     ::  test(n)   ! test for the convergence 
    real(8)                     ::  test1(n)  ! test for the convergence 
    real(8)                     ::  test2(n)  ! test for the convergence 
    real(8)                     ::  test_cr(n)  ! test for the convergence 
    real(8)                     ::  dk(n)     ! k
    real(8)                     ::  dr(n)     ! r
    real(8)                     ::  chik      ! chi
!}}}
    ! for convenience, h for H, c for C
    ! H = r*h, C = r*c
    ! h(r): the total  correlation function
    ! c(r): the direct correlation function
    ! g(r): the indirect correlation function
    ! g_(r): the pair correlation function
! normal part{{{
!********************* k space ****************************************
!k space{{{
    real(8)                     ::  hkmm(n)      !  matrix-matrix          
    real(8)                     ::  hkfm(n)      !  matrix-fluid 
    real(8)                     ::  hkffb(n)     !  fluid-fluid(connected)
    real(8)                     ::  hkffc(n)     !  fluid-fluid(block)
    real(8)                     ::  hkff(n)      !  total fluid-fluid

    real(8)                     ::  ckmm(n)      !  matrix-matrix          
    real(8)                     ::  ckfm(n)      !  matrix-fluid 
    real(8)                     ::  ckffb(n)     !  fluid-fluid(connected)
    real(8)                     ::  ckffc(n)     !  fluid-fluid(block)
    real(8)                     ::  ckff(n)      !  total fluid-fluid

    real(8)                     ::  gkmm(n)      !  matrix-matrix          
    real(8)                     ::  gkfm(n)      !  matrix-fluid 
    real(8)                     ::  gkffb(n)     !  fluid-fluid(connected)
    real(8)                     ::  gkffc(n)     !  fluid-fluid(block)
    real(8)                     ::  gkff(n)      !  total fluid-fluid

!}}}
!********************* r space ****************************************
!r space{{{
    real(8)                     ::  hrmm(n)      !  matrix-matrix          
    real(8)                     ::  hrfm(n)      !  matrix-fluid 
    real(8)                     ::  hrffb(n)     !  fluid-fluid(connected)
    real(8)                     ::  hrffc(n)     !  fluid-fluid(block)
    real(8)                     ::  hrff(n)      !  total fluid-fluid

    real(8)                     ::  crmm(n)      !  matrix-matrix          
    real(8)                     ::  crfm(n)      !  matrix-fluid 
    real(8)                     ::  crffb(n)     !  fluid-fluid(connected)
    real(8)                     ::  crffc(n)     !  fluid-fluid(block)
    real(8)                     ::  crff(n)      !  total fluid-fluid

    real(8)                     ::  grmm(n)      !  matrix-matrix          
    real(8)                     ::  grfm(n)      !  matrix-fluid 
    real(8)                     ::  grffb(n)     !  fluid-fluid(connected)
    real(8)                     ::  grffc(n)     !  fluid-fluid(block)
    real(8)                     ::  grff(n)      !  total fluid-fluid

!}}}
    real(8)                     ::  g_rmm(n)      !  matrix-matrix          
    real(8)                     ::  g_rfm(n)      !  matrix-fluid 
    real(8)                     ::  g_rffb(n)     !  fluid-fluid(connected)
    real(8)                     ::  g_rffc(n)     !  fluid-fluid(block)
    real(8)                     ::  g_rff(n)      !  total fluid-fluid

!}}}
!single particle part{{{
!************  s for single particle ********************
!********************** k space ******************************
!k space{{{
!   h(k)
    real(8)                     ::  hksm(n)      !  single-matrix          
    real(8)                     ::  hksf(n)      !  total single-fluid
    real(8)                     ::  hksfc(n)     !  single-fluid(connected)
    real(8)                     ::  hksfb(n)     !  single-single(block)
    real(8)                     ::  hkss(n)      !  total single-single 
    real(8)                     ::  hkssc(n)     !  single-single(connected)
    real(8)                     ::  hkssb(n)     !  single-single(block)
!   c(k)
    real(8)                     ::  cksm(n)      !  single-matrix          
    real(8)                     ::  cksf(n)      !  total single-fluid
    real(8)                     ::  cksfc(n)     !  single-fluid(connected)
    real(8)                     ::  cksfb(n)     !  single-single(block)
    real(8)                     ::  ckss(n)      !  total single-single 
    real(8)                     ::  ckssc(n)     !  single-single(connected)
    real(8)                     ::  ckssb(n)     !  single-single(block)
!   g(k)
    real(8)                     ::  gksm(n)      !  single-matrix          
    real(8)                     ::  gksf(n)      !  total single-fluid
    real(8)                     ::  gksfc(n)     !  single-fluid(connected)
    real(8)                     ::  gksfb(n)     !  single-single(block)
    real(8)                     ::  gkss(n)      !  total single-single 
    real(8)                     ::  gkssc(n)     !  single-single(connected)
    real(8)                     ::  gkssb(n)     !  single-single(block)
!}}}
!********************** r space ******************************
!r space{{{
!   h(r)
    real(8)                     ::  hrsm(n)      !  single-matrix          
    real(8)                     ::  hrsf(n)      !  total single-fluid
    real(8)                     ::  hrsfc(n)     !  single-fluid(connected)
    real(8)                     ::  hrsfb(n)     !  single-single(block)
    real(8)                     ::  hrss(n)      !  total single-single 
    real(8)                     ::  hrssc(n)     !  single-single(connected)
    real(8)                     ::  hrssb(n)     !  single-single(block)
!   c(r)
    real(8)                     ::  crsm(n)      !  single-matrix          
    real(8)                     ::  crsf(n)      !  total single-fluid
    real(8)                     ::  crsfc(n)     !  single-fluid(connected)
    real(8)                     ::  crsfb(n)     !  single-single(block)
    real(8)                     ::  crss(n)      !  total single-single 
    real(8)                     ::  crssc(n)     !  single-single(connected)
    real(8)                     ::  crssb(n)     !  single-single(block)
!   g(r)
    real(8)                     ::  grsm(n)      !  single-matrix          
    real(8)                     ::  grsf(n)      !  total single-fluid
    real(8)                     ::  grsfc(n)     !  single-fluid(connected)
    real(8)                     ::  grsfb(n)     !  single-single(block)
    real(8)                     ::  grss(n)      !  total single-single 
    real(8)                     ::  grssc(n)     !  single-single(connected)
    real(8)                     ::  grssb(n)     !  single-single(block)
!}}}

! pair correlation function
    real(8)                     ::  g_rsm(n)      !  single-matrix          
    real(8)                     ::  g_rsf(n)      !  total single-fluid
    real(8)                     ::  g_rsfc(n)     !  single-fluid(connected)
    real(8)                     ::  g_rsfb(n)     !  single-single(block)
    real(8)                     ::  g_rss(n)      !  total single-single 
    real(8)                     ::  g_rssc(n)     !  single-single(connected)
    real(8)                     ::  g_rssb(n)     !  single-single(block)
!}}}
! structure factors
    real(8)                     ::  skmm(n)      
    real(8)                     ::  skfm(n)      
    real(8)                     ::  skff(n)      
    real(8)                     ::  skss(n)      
    real(8)                     ::  sksf(n)      
    real(8)                     ::  sksm(n)      
    !  single particle part

    integer                     ::  serror
    integer                     ::  ierror
    integer                     ::  i
    integer                     ::  j
    integer                     ::  times
    character(20)               ::  filename 
    real(8)                     ::  xtmp
    real(8)                     ::  ytmp
    real(8)                     ::  tmp
    real(8)                     ::  ctmp
    
!}}}

end module module_common
