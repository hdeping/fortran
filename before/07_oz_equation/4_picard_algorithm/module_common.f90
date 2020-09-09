!   m:  matrice
!   f:  fluid 
module module_common

    implicit none 
    integer,parameter           ::   l    = 5
    integer,parameter           ::   n    = 2**l
    integer,parameter           ::  fre   =  int(1E3)
    real(8),parameter           ::   pi   = 3.141592653 
    real(8),parameter           :: deltar = 0.04
    real(8),parameter           :: deltak = pi/dble(n)/deltar 
    real(8),parameter           :: error  = 1E-8      !  for the differences
    real(8),parameter           :: dmm    = 1.0       !  m-m 
    real(8),parameter           :: dfm    = 1.0       !  f-m
    real(8),parameter           :: dff    = 1.0       !  f-f
    real(8),parameter           :: rhom   = 0.6       !  the density of matrix
    real(8),parameter           :: rhof   = 0.1       !  the density of fluid 
    !  variables for fft
    integer                     :: status 
    !type(dfti_descriptor), pointer :: my_desc1_handle
    !type(dfti_descriptor), pointer :: my_desc2_handle
    ! variables for OZ equation  
    real(8)                     ::  k
    real(8)                     ::  r
    real(8)                     ::  t1
    real(8)                     ::  t2
    real(8)                     ::  t
    real(8)                     ::  lambda    ! judge the convergence
    real(8),dimension(0:(n-1))  ::  maymm(n)  ! test for the convergence
    real(8),dimension(0:(n-1))  ::  mayfm(n)  ! test for the convergence
    real(8),dimension(0:(n-1))  ::  mayff(n)  ! test for the convergence
    real(8),dimension(0:(n-1))  ::  test(n)   ! mayer function 
    real(8),dimension(0:(n-1))  ::  dk(n)     ! k
    real(8),dimension(0:(n-1))  ::  dr(n)     ! r
    real(8),dimension(0:(n-1))  ::  chik      ! chi

    ! for convenience, h for H, c for C
    ! H = r*h, C = r*c
! the total correlation function
    !  k space  ,h(k)
    real(8),dimension(0:(n-1))  ::  hkmm(n)      !  matrix-matrix          
    real(8),dimension(0:(n-1))  ::  hkfm(n)      !  matrix-fluid 
    real(8),dimension(0:(n-1))  ::  hkffb(n)     !  fluid-fluid(connected)
    real(8),dimension(0:(n-1))  ::  hkffc(n)     !  fluid-fluid(block)
    !  r space   h(r)
    real(8),dimension(0:(n-1))  ::  hrmm(n)      !  matrix-matrix          
    real(8),dimension(0:(n-1))  ::  hrfm(n)      !  matrix-fluid 
    real(8),dimension(0:(n-1))  ::  hrffb(n)     !  fluid-fluid(connected)
    real(8),dimension(0:(n-1))  ::  hrffc(n)     !  fluid-fluid(block)
! the indirect correlation function
    !  k space  gamma(k)
    real(8),dimension(0:(n-1))  ::  gkmm(n)      !  matrix-matrix
    real(8),dimension(0:(n-1))  ::  gkfm(n)      !  matrix-fluid 
    real(8),dimension(0:(n-1))  ::  gkffb(n)     !  fluid-fluid(connected)
    real(8),dimension(0:(n-1))  ::  gkffc(n)     !  fluid-fluid(block)
    real(8),dimension(0:(n-1))  ::  gkff(n)      !  total fluid-fluid
    !  r space  gamma(r)
    real(8),dimension(0:(n-1))  ::  grmm(n)      !  matrix-matrix
    real(8),dimension(0:(n-1))  ::  grfm(n)      !  matrix-fluid 
    real(8),dimension(0:(n-1))  ::  grffb(n)     !  fluid-fluid(connected)
    real(8),dimension(0:(n-1))  ::  grffc(n)     !  fluid-fluid(block)
    real(8),dimension(0:(n-1))  ::  grff(n)      !  total fluid-fluid

! the direct correlation function
    !  k space   c(k)
    real(8),dimension(0:(n-1))  ::  ckmm(n)      !  matrix-matrix
    real(8),dimension(0:(n-1))  ::  ckfm(n)      !  matrix-fluid 
    real(8),dimension(0:(n-1))  ::  ckff(n)      !  fluid-fluid
    real(8),dimension(0:(n-1))  ::  ckffb(n)     !  fluid-fluid(connected)
    real(8),dimension(0:(n-1))  ::  ckffc(n)     !  fluid-fluid(block)
    !  r space   c(r)
    real(8),dimension(0:(n-1))  ::  crmm(n)      !  matrix-matrix
    real(8),dimension(0:(n-1))  ::  crfm(n)      !  matrix-fluid 
    real(8),dimension(0:(n-1))  ::  crff(n)      !  fluid-fluid
    real(8),dimension(0:(n-1))  ::  crffb(n)     !  fluid-fluid(connected)
    real(8),dimension(0:(n-1))  ::  crffc(n)     !  fluid-fluid(block)
! pair correlation function
    real(8),dimension(0:(n-1))  ::  g_rmm(n)      !  matrix-matrix
    real(8),dimension(0:(n-1))  ::  g_rfm(n)      !  matrix-fluid 
    real(8),dimension(0:(n-1))  ::  g_rff(n)      !  fluid-fluid
    real(8),dimension(0:(n-1))  ::  g_rffb(n)     !  fluid-fluid(connected)
    real(8),dimension(0:(n-1))  ::  g_rffc(n)     !  fluid-fluid(block)

    integer                     ::  i        
    integer                     ::  j        
    integer                     ::  times        
    character(20)               ::  filename 

end module module_common
