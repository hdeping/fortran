!   m:  matrice
!   f:  fluid 
module module_common

    implicit none 
!*********   variables  ********************************
!variables{{{
    integer,parameter           ::   l    = 10
    integer,parameter           ::   n    = 2**l
    integer,parameter           ::  fre   =  int(1E6)
    real(8),parameter           ::   pi   = 3.141592653 
    !real(8),parameter           ::  top   = 10.24
    real(8),parameter           :: deltar = 0.01
    real(8),parameter           :: deltak = pi/deltar/dble(n)
    real(8),parameter           :: error  = 1E-10              !  for the differences
    real(8),parameter           :: dmm    = 1.0                !  m-m 
    real(8),parameter           :: dff    = 1.0                !  f-f
    real(8),parameter           :: dfm    = (dmm + dff)/2.0    !  f-m
    real(8),parameter           :: rhom   = 1.1              !  the density of matrix
    real(8),parameter           :: rhof   = 1.0                !  the density of fluid 
    real(8),parameter           :: gold   = (sqrt(5.0) - 1.0)/2.0  ! golden rate
    !  variables for fft
    integer                     :: status 
    !type(dfti_descriptor), pointer :: my_desc1_handle
    !type(dfti_descriptor), pointer :: my_desc2_handle
    ! variables for OZ equation  
    real(8)                     ::  rho  ! new density
    real(8)                     ::  k
    real(8)                     ::  r
    real(8)                     ::  t1
    real(8)                     ::  t2
    real(8)                     ::  eta
    real(8)                     ::  lambda1  ! hardsphere c(x)
    real(8)                     ::  lambda2  ! hardsphere c(x)
    real(8)                     ::  lambda    ! judge the convergence
    real(8)                     ::  maymm(n)  ! mayer function for m-m 
    real(8)                     ::  mayfm(n)  ! mayer function for f-m 
    real(8)                     ::  mayff(n)  ! mayer function for f-f 
    real(8)                     ::  test(n)   ! test for the convergence 
    real(8)                     ::  test1(n)  ! test for the convergence 
    real(8)                     ::  test2(n)  ! test for the convergence 
    real(8)                     ::  test_cr(n)  ! test for the convergence 
    real(8)                     ::  dk(n)     ! k
    real(8)                     ::  dr(n)     ! r
    real(8)                     ::  chik      ! chi

    ! for convenience, h for H, c for C
    ! H = r*h, C = r*c
! the total correlation function
    !  k space  ,h(k)
    real(8)                     ::  hkmm(n)      !  matrix-matrix          
    real(8)                     ::  hkfm(n)      !  matrix-fluid 
    real(8)                     ::  hkffb(n)     !  fluid-fluid(connected)
    real(8)                     ::  hkffc(n)     !  fluid-fluid(block)
    !  r space   h(r)
    real(8)                     ::  hrmm(n)      !  matrix-matrix          
    real(8)                     ::  hrfm(n)      !  matrix-fluid 
    real(8)                     ::  hrffb(n)     !  fluid-fluid(connected)
    real(8)                     ::  hrffc(n)     !  fluid-fluid(block)
! the indirect correlation function
    !  k space  gamma(k)
    real(8)                     ::  gkmm(n)      !  matrix-matrix
    real(8)                     ::  gkfm(n)      !  matrix-fluid 
    real(8)                     ::  gkffb(n)     !  fluid-fluid(connected)
    real(8)                     ::  gkffc(n)     !  fluid-fluid(block)
    real(8)                     ::  gkff(n)      !  total fluid-fluid
    !  r space  gamma(r)
    real(8)                     ::  grmm(n)      !  matrix-matrix
    real(8)                     ::  grfm(n)      !  matrix-fluid 
    real(8)                     ::  grffb(n)     !  fluid-fluid(connected)
    real(8)                     ::  grffc(n)     !  fluid-fluid(block)
    real(8)                     ::  grff(n)      !  total fluid-fluid

! the direct correlation function
    !  k space   c(k)
    real(8)                     ::  ckmm(n)      !  matrix-matrix
    real(8)                     ::  ckfm(n)      !  matrix-fluid 
    real(8)                     ::  ckff(n)      !  fluid-fluid
    real(8)                     ::  ckffb(n)     !  fluid-fluid(connected)
    real(8)                     ::  ckffc(n)     !  fluid-fluid(block)
    !  r space   c(r)
    real(8)                     ::  crmm(n)      !  matrix-matrix
    real(8)                     ::  crfm(n)      !  matrix-fluid 
    real(8)                     ::  crff(n)      !  fluid-fluid
    real(8)                     ::  crffb(n)     !  fluid-fluid(connected)
    real(8)                     ::  crffc(n)     !  fluid-fluid(block)
! pair correlation function
    real(8)                     ::  g_rmm(n)      !  matrix-matrix
    real(8)                     ::  g_rfm(n)      !  matrix-fluid 
    real(8)                     ::  g_rff(n)      !  fluid-fluid
    real(8)                     ::  g_rffb(n)     !  fluid-fluid(connected)
    real(8)                     ::  g_rffc(n)     !  fluid-fluid(block)

    integer                     ::  jj
    integer                     ::  i
    integer                     ::  j        
    integer                     ::  times        
    integer                     ::  ierror
    real(8)                     ::  xtmp
    real(8)                     ::  ytmp
    real(8)                     ::  tmp
    real(8)                     ::  ctmp
    character(20)               ::  filename 
!}}}

end module module_common
