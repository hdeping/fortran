!   m:  matrice
!   f:  fluid 
module module_common

    implicit none 
!variables{{{
    integer,parameter           ::   l    = 6
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
    complex(8)                  :: x(n)
    complex(8)                  :: y(n)
    ! variables for OZ equation  
    real(8)                     ::  k
    real(8)                     ::  r
    real(8)                     ::  t1
    real(8)                     ::  t2
    real(8)                     ::  lambda    ! judge the convergence 
    real(8)                     ::  maymm(n)  ! test for the convergence
    real(8)                     ::  mayfm(n)  ! test for the convergence
    real(8)                     ::  mayff(n)  ! test for the convergence
    real(8)                     ::  test(n)   ! mayer function 
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
    real(8)                     ::  grff(n)      !  total 
! the direct correlation function
    !  k space   c(k)
    real(8)                     ::  ckmm(n)      !  matrix-matrix
    real(8)                     ::  ckfm(n)      !  matrix-fluid 
    real(8)                     ::  ckffb(n)     !  fluid-fluid(connected)
    real(8)                     ::  ckffc(n)     !  fluid-fluid(block)
    real(8)                     ::  ckff(n)      !  total  
    !  r space   c(r)
    real(8)                     ::  crmm(n)      !  matrix-matrix
    real(8)                     ::  crfm(n)      !  matrix-fluid 
    real(8)                     ::  crff(n)      !  fluid-fluid
    real(8)                     ::  crffb(n)     !  fluid-fluid(connected)
    real(8)                     ::  crffc(n)     !  fluid-fluid(block)
! pair correlation function
    real(8)                     ::  g_rmm(n)      !  matrix-matrix
    real(8)                     ::  g_rfm(n)      !  matrix-fluid 
    real(8)                     ::  g_rffc(n)      !  fluid-fluid
    real(8)                     ::  g_rffb(n)     !  fluid-fluid(connected)
    real(8)                     ::  g_rff(n)      !  total 

    integer                     ::  i        
    integer                     ::  j        
    integer                     ::  times        
    character(20)               ::  filename 

    integer                     ::  itmp
    real(8)                     ::  a(n)
    real(8)                     ::  b(n)
    real(8)                     ::  c(n)
    real(8)                     ::  a1(n)
    real(8)                     ::  a2(n)
    real(8)                     ::  a3(n)
    real(8)                     ::  a4(n)
    complex(8)                  ::  b1(n)
    complex(8)                  ::  b2(n)
!}}}


end module module_common
