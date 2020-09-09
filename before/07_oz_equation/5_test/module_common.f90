!   m:  matrice
!   f:  fluid 
module module_common
    !use mkl_dfti
    implicit none 
    integer,parameter              ::   l   = 20
    integer,parameter              ::   n   = 2**7 
    integer,parameter              ::   m   = n + 2 
    real(8),parameter              :: delta = 0.05
    real(8),parameter              :: error = 1E-8      !  for the differences
    real(8),parameter              :: dmm   = 0.8       !  m-m 
    real(8),parameter              :: dfm   = 0.8       !  f-m
    real(8),parameter              :: dff   = 0.8       !  f-f
    real(8)                        :: rhom  = 0.4       !  the density of matrix
    real(8)                        :: rhof  = 0.4       !  the density of fluid 
    integer,parameter              ::  fre  = int(n/10000)
    !  variables for fft
    integer                        :: status 
    !type(dfti_descriptor), pointer :: my_desc1_handle
    !type(dfti_descriptor), pointer :: my_desc2_handle
    ! variables for OZ equation  
    real(8)                        ::  q
    real(8)                        ::  r
    real(8)                        ::  lambda    ! judge the convergence
    real(8)                     ::  test(n)      ! test for the convergence
    real(8)                     ::  mayer(n)     ! mayer function 
    ! for convenience, h for H, c for C
    ! H = r*h, C = r*c
    ! the total correlation function
    !  q space  ,h(q)
    real(8)                     ::  hqmm(n)      !  matrix-matrix
    real(8)                     ::  hqfm(n)      !  matrix-fluid 
    real(8)                     ::  hqffb(n)     !  fluid-fluid(connected)
    real(8)                     ::  hqffc(n)     !  fluid-fluid(block)
    !  r space   h(r)
    real(8)                     ::  hrmm(n)      !  matrix-matrix
    real(8)                     ::  hrfm(n)      !  matrix-fluid 
    real(8)                     ::  hrffb(n)     !  fluid-fluid(connected)
    real(8)                     ::  hrffc(n)     !  fluid-fluid(block)

    ! the indirect correlation function
    !  q space  gamma(q)
    real(8)                     ::  gqmm(n)      !  matrix-matrix
    real(8)                     ::  gqfm(n)      !  matrix-fluid 
    real(8)                     ::  gqffb(n)     !  fluid-fluid(connected)
    real(8)                     ::  gqffc(n)     !  fluid-fluid(block)
    !  r space  gamma(r)
    real(8)                     ::  grmm(n)      !  matrix-matrix
    real(8)                     ::  grfm(n)      !  matrix-fluid 
    real(8)                     ::  grffb(n)     !  fluid-fluid(connected)
    real(8)                     ::  grffc(n)     !  fluid-fluid(block)
    ! the direct correlation function
    !  q space   c(q)
    real(8)                     ::  cqmm(n)      !  matrix-matrix
    real(8)                     ::  cqfm(n)      !  matrix-fluid 
    real(8)                     ::  cqffb(n)     !  fluid-fluid(connected)
    real(8)                     ::  cqffc(n)     !  fluid-fluid(block)
    !  r space   c(r)
    real(8)                     ::  crmm(n)      !  matrix-matrix
    real(8)                     ::  crfm(n)      !  matrix-fluid 
    real(8)                     ::  crff(n)      !  fluid-fluid
    real(8)                     ::  crffb(n)     !  fluid-fluid(connected)
    real(8)                     ::  crffc(n)     !  fluid-fluid(block)

    integer                        ::  i        
    integer                        ::  j        
    integer                        ::  times        
    character(20)                  ::  filename 

    contains
!!subroutine fft{{{
!      subroutine fft(y,x,n)
!          integer,intent(in)           :: n
!          real(8),intent(in)        :: x(n)
!          real(8),intent(out)       :: y(n)
!      y = x
!      status = dfticreatedescriptor( my_desc1_handle, dfti_single, &
!          dfti_real, 1, n) 
!      status = dfticommitdescriptor( my_desc1_handle ) 
!      status = dfticomputeforward( my_desc1_handle, y ) 
!      status = dftifreedescriptor(my_desc1_handle) 
!      end subroutine fft
!!}}}
!!subroutine ifft{{{
!      subroutine ifft(y,x,n)
!          integer,intent(in)           :: n
!          real(8),intent(in)        :: x(n)
!          real(8),intent(out)       :: y(n)
!      y = x
!      status = dfticreatedescriptor( my_desc1_handle, dfti_single, &
!          dfti_real, 1, n) 
!      status = dfticommitdescriptor( my_desc1_handle ) 
!      status = dfticomputebackward( my_desc1_handle, y ) 
!      status = dftifreedescriptor(my_desc1_handle) 
!      y = y/n
!      end subroutine ifft
!!}}}
!subroutine may{{{
    subroutine may(mayer,d)
        real(8)         ::  mayer(n)
        real(8)         ::  d
        real(8)         ::  x
        integer         ::  itmp
        do itmp = 1,n
            x = i*delta
            if(x <= d)then
                mayer(i) = -1
            else
                mayer(i) = 0
            endif
        end do
    end subroutine may
!}}}
!subroutine converg{{{
    !  convergence judgement
    subroutine converg(lambda,a,b)
        real(8),intent(out)           :: lambda
        real(8),intent(in)            :: a(n) 
        real(8),intent(in)            :: b(n) 
        integer                       :: itmp
        real(8)                       :: tmp 
        
        tmp = 0
        do i = 1,n
           tmp = tmp + abs(a(i) - b(i))
        end do
        tmp = tmp/n
        lambda = tmp

    end subroutine converg
!}}}
    
     
end module module_common
