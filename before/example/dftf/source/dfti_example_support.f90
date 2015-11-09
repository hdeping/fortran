!*****************************************************************************
!                              INTEL CONFIDENTIAL
! Copyright(C) 2003-2009 Intel Corporation. All Rights Reserved.
! The source code contained  or  described herein and all documents related to
! the source code ("Material") are owned by Intel Corporation or its suppliers
! or licensors.  Title to the  Material remains with  Intel Corporation or its
! suppliers and licensors. The Material contains trade secrets and proprietary
! and  confidential  information of  Intel or its suppliers and licensors. The
! Material  is  protected  by  worldwide  copyright  and trade secret laws and
! treaty  provisions. No part of the Material may be used, copied, reproduced,
! modified, published, uploaded, posted, transmitted, distributed or disclosed
! in any way without Intel's prior express written permission.
! No license  under any  patent, copyright, trade secret or other intellectual
! property right is granted to or conferred upon you by disclosure or delivery
! of the Materials, either expressly, by implication, inducement, estoppel or
! otherwise.  Any  license  under  such  intellectual property  rights must be
! express and approved by Intel in writing.
!
!*****************************************************************************
! Content:
!       MKL DFTI example support functions (Fortran-interface)
!
!*****************************************************************************

      subroutine ZERO_INIT_COMPLEX_C(X, n)
!     input parameters
      COMPLEX(4) :: X(*)
      integer,intent(in) :: n

!     Local parameters
      integer   i

      do i=1,n
          X(i)=cmplx(0.0,0.0)
      end do
      end subroutine ZERO_INIT_COMPLEX_C

!---------------
      subroutine ZERO_INIT_REAL_S(X, n)
!     input parameters
      REAL(4):: X(*)
      integer,intent(in) :: n

!     Local parameters
      integer   i

      do i=1,n
          X(i)=0.0
      end do
      end subroutine ZERO_INIT_REAL_S

!---------------
      subroutine ZERO_INIT_COMPLEX_Z(X, n)
!     input parameters
      COMPLEX(8) :: X(*)
      integer,intent(in) :: n

!     Local parameters
      integer   i

      do i=1,n
          X(i)=dcmplx(0.0,0.0)
      end do
      end subroutine ZERO_INIT_COMPLEX_Z

!---------------
      subroutine ZERO_INIT_REAL_D(X, n)
!     input parameters
      REAL(8):: X(*)
      integer,intent(in) :: n

!     Local parameters
      integer   i

      do i=1,n
          X(i)=0.0
      end do
      end subroutine ZERO_INIT_REAL_D

!---------------
      subroutine INIT_FORW_TONE_AND_EXP_RES_Z(X_IN, X_EXP, n)
!     input parameters
      COMPLEX(8) :: X_IN(*)
      COMPLEX(8) :: X_EXP(*)
      integer,intent(in) :: n

!     Local parameters
      REAL(KIND=8)    MATH_PI
      integer         i
      REAL(KIND=8)    f_step

      if (n .EQ. 1) RETURN

      MATH_PI = real(1, KIND=8)
      MATH_PI = real(4, KIND=8) * atan(MATH_PI)

      do i=1,n
          f_step = real(2, KIND=8) * MATH_PI * real(i-1, KIND=8) / real(n, KIND=8)
          X_IN(i)=dcmplx( dcos(f_step), dsin(f_step) )
      end do

      do i=1,n
          X_EXP(i) = dcmplx( 0.0, 0.0 )
      end do
      X_EXP(2) = dcmplx( 1.0, 0.0 )
      end subroutine INIT_FORW_TONE_AND_EXP_RES_Z

!---------------
      subroutine INIT_BACKW_TONE_AND_EXP_RES_Z(X_IN, X_EXP, n)
!     input parameters
      COMPLEX(8) :: X_IN(*)
      COMPLEX(8) :: X_EXP(*)
      integer,intent(in) :: n

!     Local parameters
      REAL(KIND=8)    MATH_PI
      integer        i

      REAL(KIND=8)   f_step

      if (n .EQ. 1) RETURN

      MATH_PI = real(1, KIND=8)
      MATH_PI = real(4, KIND=8) * atan(MATH_PI)

      do i=1,n
          f_step = (real(2, KIND=8) * MATH_PI * real(i-1, KIND=8)) / real(n, KIND=8)
          X_IN(i)=dcmplx( dcos(f_step), -dsin(f_step) )
      end do

      do i=1,n
          X_EXP(i) = dcmplx( 0.0, 0.0 )
      end do
      X_EXP(2) = dcmplx( 1.0, 0.0 )
      end subroutine INIT_BACKW_TONE_AND_EXP_RES_Z

!---------------
      subroutine INIT_COMPLEX_VECTOR_Z(X, n)
!     input parameters
      COMPLEX(8) :: X(*)
      integer,intent(in) :: n

!     Local parameters
      integer   i
      real      f_step, pi
      parameter (pi=3.1416)

      do i=1,n
          f_step = float(i)
          X(i)=dcmplx((sin(f_step)*sqrt(3.0))/float(2),sin(f_step) &
               /sqrt(3.0))
      end do
      end subroutine INIT_COMPLEX_VECTOR_Z

!---------------
      subroutine INIT_REAL_VECTOR_D(X, n)
!     input parameters
      REAL(8):: X(*)
      integer,intent(in) :: n

!     Local parameters
      integer   i
      real      f_step, pi
      parameter (pi=3.1416)

      do i=1,n
          f_step = float(i)
          X(i)=sin(f_step)*sqrt(3.0)/float(2)
      end do
      end subroutine INIT_REAL_VECTOR_D

!---------------
      subroutine INIT_MULTIPLE_VECTOR_Z(X, n, multiple, dist, strides)
!     input parameters
      COMPLEX(8) :: X(*)
      integer,intent(in) ::  n
      integer,intent(in) ::  multiple
      integer,intent(in) ::  dist
      integer,intent(in) ::  strides(*)

!     Local parameters
      integer   i, j
      integer   first_index
      integer   stride_x

      real      f_step, pi
      parameter (pi=3.1416)

      first_index = strides(1)
      stride_x    = strides(2)

      do j=1, multiple*dist, dist
          do i=j+first_index, j+first_index+n*stride_x-1, stride_x
              f_step = float(i)* float(j)
              X(i)=dcmplx((sin(f_step)*sqrt(3.0))/float(2),sin(f_step) &
                   /sqrt(3.0))
          end do
      end do
      end subroutine INIT_MULTIPLE_VECTOR_Z

!---------------
      subroutine INIT_MULTIPLE_REAL_VECTOR_D(X, n, multiple, dist, strides)
!     input parameters
      REAL(8):: X(*)
      integer,intent(in) :: n
      integer,intent(in) :: multiple
      integer,intent(in) :: dist
      integer,intent(in) :: strides(*)

!     Local parameters
      integer   i, j
      integer   first_index
      integer   stride_x

      integer   k
      real      f_step, pi
      parameter (pi=3.1416)

      first_index = strides(1)
      stride_x    = strides(2)

      k=1
      do j=1, multiple*dist, dist
          do i=j+first_index, j+first_index+n*stride_x-1, stride_x
                k=k+1
               f_step = float(k)
!              f_step = float(i)* float(j)
              X(i)=sin(f_step)*sqrt(3.0)/float(2)
          end do
      end do
      end subroutine INIT_MULTIPLE_REAL_VECTOR_D

!---------------
      subroutine INIT_MULTIPLE_2D_COLUMNS_Z(X, m, n, strides)
!     input parameters
      COMPLEX(8) :: X(*)
      integer,intent(in) :: m
      integer,intent(in) :: n
      integer,intent(in) :: strides(*)

!     Local parameters
      integer   i, j
      integer   first_index
      integer   stride_m
      integer   stride_n

      real      f_step, pi
      parameter (pi=3.1416)

      first_index = strides(1)
      stride_m    = strides(2)
      stride_n    = strides(3)

      do j=1+first_index, n*stride_n+first_index, stride_n
          do i=j, j+m*stride_m-1, stride_m
              f_step = float(i)* float(j)
              X(i)=dcmplx((sin(f_step)*sqrt(3.0))/float(2),sin(f_step) &
                   /sqrt(3.0))
          end do
      end do
      end subroutine INIT_MULTIPLE_2D_COLUMNS_Z

!---------------
      subroutine INIT_3D_COLUMNS_Z(X, m, n, k, strides)
!     input parameters
      COMPLEX(8) :: X(*)
      integer,intent(in) :: m
      integer,intent(in) :: n
      integer,intent(in) :: k
      integer,intent(in) :: strides(*)

!     Local parameters
      integer   i, j, l
      integer   first_index
      integer   stride_m
      integer   stride_n
      integer   stride_k

      real      f_step, pi
      parameter (pi=3.1416)

      first_index = strides(1)
      stride_m    = strides(2)
      stride_n    = strides(3)
      stride_k    = strides(4)

      do l=1+first_index, k*stride_k+first_index, stride_k
          do j=l, l+n*stride_n-1, stride_n
              do i=j, j+m*stride_m-1, stride_m
                  f_step = float(i)* float(j)* float(l)
                  X(i)=dcmplx((sin(f_step)*sqrt(3.0))/float(2),sin(f_step) &
                       /sqrt(3.0))
              end do
          end do
      end do
      end subroutine INIT_3D_COLUMNS_Z

!---------------
      subroutine INIT_3D_COLUMNS_D(X, m, n, k, strides)
!     input parameters
      REAL(8) :: X(*)
      integer,intent(in) :: m
      integer,intent(in) :: n
      integer,intent(in) :: k
      integer,intent(in) :: strides(*)

!     Local parameters
      integer   i, j, l
      integer   first_index
      integer   stride_m
      integer   stride_n
      integer   stride_k

      real      f_step, pi
      parameter (pi=3.1416)

      first_index = strides(1)
      stride_m    = strides(2)
      stride_n    = strides(3)
      stride_k    = strides(4)

      do l=1+first_index, k*stride_k+first_index, stride_k
          do j=l, l+n*stride_n-1, stride_n
              do i=j, j+m*stride_m-1, stride_m
                  f_step = float(i)* float(j)* float(l)
                  X(i) = sin(f_step)*sqrt(3.0)/float(2)
              end do
          end do
      end do
      end subroutine INIT_3D_COLUMNS_D

!---------------
      subroutine INIT_FORW_TONE_AND_EXP_RES_C(X_IN, X_EXP, n)
!     input parameters
      COMPLEX(4) :: X_IN(*)
      COMPLEX(4) :: X_EXP(*)
      integer,intent(in) :: n

!     Local parameters
      REAL(KIND=8)    MATH_PI
      PARAMETER       (MATH_PI = 3.14159265358979323846)
      integer   i
      real      f_step

      if (n .EQ. 1) RETURN

      do i=1,n
          f_step = 2 * MATH_PI * real((i-1), KIND=4)/real(n, KIND=4)
          X_IN(i)=cmplx( cos(f_step), sin(f_step) )
      end do

      do i=1,n
          X_EXP(i) = cmplx( 0.0, 0.0 )
      end do
      X_EXP(2) = cmplx( 1.0, 0.0 )
      end subroutine INIT_FORW_TONE_AND_EXP_RES_C

!---------------
      subroutine INIT_BACKW_TONE_AND_EXP_RES_C(X_IN, X_EXP, n)
!     input parameters
      COMPLEX(4) :: X_IN(*)
      COMPLEX(4) :: X_EXP(*)
      integer,intent(in) :: n

!     Local parameters
      REAL(KIND=8)    MATH_PI
      PARAMETER       (MATH_PI = 3.14159265358979323846)
      integer   i
      real      f_step

      if (n .EQ. 1) RETURN

      do i=1,n
          f_step = 2 * MATH_PI * float(i-1) / float(n)
          X_IN(i)=cmplx( cos(f_step), -sin(f_step) )
      end do

      do i=1,n
          X_EXP(i) = cmplx( 0.0, 0.0 )
      end do
      X_EXP(2) = cmplx( 1.0, 0.0 )
      end subroutine INIT_BACKW_TONE_AND_EXP_RES_C

!---------------
      subroutine INIT_COMPLEX_VECTOR_C(X, n)
!     input parameters
      COMPLEX(4) :: X(*)
      integer,intent(in) :: n

!     Local parameters
      integer   i
      real      f_step, pi
      parameter (pi=3.1416)

      do i=1,n
          f_step = float(i)
          X(i)=cmplx((sin(f_step)*sqrt(3.0))/float(2),sin(f_step) &
               /sqrt(3.0))
      end do
      end subroutine INIT_COMPLEX_VECTOR_C

!---------------
      subroutine INIT_REAL_VECTOR_S(X, n)
!     input parameters
      real(4):: X(*)
      integer,intent(in) :: n

!     Local parameters
      integer   i
      real      f_step, pi
      parameter (pi=3.1416)

      do i=1,n
          f_step = float(i)
          X(i)=sin(f_step)*sqrt(3.0)/float(2)
      end do
      end subroutine INIT_REAL_VECTOR_S

!---------------
      subroutine INIT_MULTIPLE_VECTOR_C(X, n, multiple, dist, strides)
!     input parameters
      COMPLEX(4) :: X(*)
      integer,intent(in) :: n
      integer,intent(in) :: multiple
      integer,intent(in) :: dist
      integer,intent(in) :: strides(*)

!     Local parameters
      integer   i, j
      integer   first_index
      integer   stride_x

      real      f_step, pi
      parameter (pi=3.1416)

      first_index = strides(1)
      stride_x    = strides(2)

      do j=1, multiple*dist, dist
          do i=j+first_index, j+first_index+n*stride_x-1, stride_x
              f_step = float(i)* float(j)
              X(i)=cmplx((sin(f_step)*sqrt(3.0))/float(2),sin(f_step) &
                   /sqrt(3.0))
          end do
      end do
      end subroutine INIT_MULTIPLE_VECTOR_C

!---------------
      subroutine INIT_MULTIPLE_REAL_VECTOR_S(X, n, multiple, dist, strides)
!     input parameters
      REAL(4):: X(*)
      integer,intent(in) :: n
      integer,intent(in) :: multiple
      integer,intent(in) :: dist
      integer,intent(in) :: strides(*)

!     Local parameters
      integer   i, j
      integer   first_index
      integer   stride_x
      integer   k

      real      f_step

      first_index = strides(1)
      stride_x    = strides(2)

      k=1
      do j=1, multiple*dist, dist
          do i=j+first_index, j+first_index+n*stride_x-1, stride_x
!              f_step = float(i)* float(j)
              f_step = float(k)
              X(i)=sin(f_step)*sqrt(3.0)/float(2)
           k=k+1
         end do
      end do
      end subroutine INIT_MULTIPLE_REAL_VECTOR_S

!---------------
      subroutine INIT_MULTIPLE_2D_COLUMNS_S(X, m, n, strides)
!     input parameters
      REAL(4) :: X(*)
      integer,intent(in) :: m
      integer,intent(in) :: n
      integer,intent(in) :: strides(*)

!     Local parameters
      integer   i, j
      integer   first_index
      integer   stride_m
      integer   stride_n

      real      f_step

      first_index = strides(1)
      stride_m    = strides(2)
      stride_n    = strides(3)

      do j=1+first_index, n*stride_n+first_index, stride_n
          do i=j, j+m*stride_m-1, stride_m
              f_step = float(i)* float(j)
              X(i)=sin(f_step)*sqrt(3.0)/float(2)
          end do
      end do
      end subroutine INIT_MULTIPLE_2D_COLUMNS_S

!---------------
      subroutine INIT_MULTIPLE_2D_COLUMNS_D(X, m, n, strides)
!     input parameters
      REAL(8) :: X(*)
      integer,intent(in) :: m
      integer,intent(in) :: n
      integer,intent(in) :: strides(*)

!     Local parameters
      integer   i, j
      integer   first_index
      integer   stride_m
      integer   stride_n

      REAL(8)   f_step

      first_index = strides(1)
      stride_m    = strides(2)
      stride_n    = strides(3)

      do j=1+first_index, n*stride_n+first_index, stride_n
          do i=j, j+m*stride_m-1, stride_m
              f_step = float(i)* float(j)
              X(i)=sin(f_step)*sqrt(3.0)/float(2)
          end do
      end do
      end subroutine INIT_MULTIPLE_2D_COLUMNS_D

!---------------
      subroutine INIT_MULTIPLE_2D_COLUMNS_C(X, m, n, strides)
!     input parameters
      COMPLEX(4) :: X(*)
      integer,intent(in) :: m
      integer,intent(in) :: n
      integer,intent(in) :: strides(*)

!     Local parameters
      integer   i, j
      integer   first_index
      integer   stride_m
      integer   stride_n

      real      f_step, pi
      parameter (pi=3.1416)

      first_index = strides(1)
      stride_m    = strides(2)
      stride_n    = strides(3)

      do j=1+first_index, n*stride_n+first_index, stride_n
          do i=j, j+m*stride_m-1, stride_m
              f_step = float(i)* float(j)
              X(i)=cmplx((sin(f_step)*sqrt(3.0))/float(2),sin(f_step) &
                   /sqrt(3.0))
          end do
      end do
      end subroutine INIT_MULTIPLE_2D_COLUMNS_C

!---------------
      subroutine INIT_3D_COLUMNS_C(X, m, n, k, strides)
!     input parameters
      COMPLEX(4) :: X(*)
      integer,intent(in) :: m
      integer,intent(in) :: n
      integer,intent(in) :: k
      integer,intent(in) :: strides(*)

!     Local parameters
      integer   i, j, l
      integer   first_index
      integer   stride_m
      integer   stride_n
      integer   stride_k

      real      f_step, pi
      parameter (pi=3.1416)

      first_index = strides(1)
      stride_m    = strides(2)
      stride_n    = strides(3)
      stride_k    = strides(4)

      do l=1+first_index, k*stride_k+first_index, stride_k
          do j=l, l+n*stride_n-1, stride_n
              do i=j, j+m*stride_m-1, stride_m
                  f_step = float(i)* float(j)* float(l)
                  X(i)=cmplx((sin(f_step)*sqrt(3.0))/float(2),sin(f_step) &
                       /sqrt(3.0))
              end do
          end do
      end do
      end subroutine INIT_3D_COLUMNS_C

!---------------
      subroutine INIT_3D_COLUMNS_S(X, m, n, k, strides)
!     input parameters
      REAL(4) :: X(*)
      integer,intent(in) :: m
      integer,intent(in) :: n
      integer,intent(in) :: k
      integer,intent(in) :: strides(*)

!     Local parameters
      integer   i, j, l
      integer   first_index
      integer   stride_m
      integer   stride_n
      integer   stride_k

      real      f_step, pi
      parameter (pi=3.1416)

      first_index = strides(1)
      stride_m    = strides(2)
      stride_n    = strides(3)
      stride_k    = strides(4)

      do l=1+first_index, k*stride_k+first_index, stride_k
          do j=l, l+n*stride_n-1, stride_n
              do i=j, j+m*stride_m-1, stride_m
                  f_step = float(i)* float(j)* float(l)
                  X(i)=sin(f_step)*sqrt(3.0)/float(2)
              end do
          end do
      end do
      end subroutine INIT_3D_COLUMNS_S

!---------------
      subroutine PRINT_VECTOR_C( X, n)
!     input parameters
      COMPLEX(4),INTENT(IN) :: X(*)
      integer,intent(in) :: n

!     Local parameters
      integer   i

      do i=1,n
         print 902, i, X(i)
      end do

902   format(' X(', I4, ')= (', F8.3, ',', F8.3, ')')
      end subroutine PRINT_VECTOR_C

!---------------
      subroutine PRINT_VECTOR_S( X, n)
!     input parameters
      REAL(4),INTENT(IN) :: X(*)
      integer,intent(in) :: n

!     Local parameters
      integer   i

      do i=1,n
         print 902, i, X(i)
      end do

902   format(' X(', I4, ')= ', F8.3)
      end subroutine PRINT_VECTOR_S

!---------------
      subroutine PRINT_THREE_VECTORS_C(X, n, dist, strides)
!     input parameters
      COMPLEX(4),INTENT(IN) :: X(*)
      integer,intent(in) :: n
      integer,intent(in) :: dist
      integer,intent(in) :: strides(*)

!     Local parameters
      integer   i, j
      integer   first_index
      integer   stride_x

      first_index = strides(1)
      stride_x    = strides(2)

      do i=1, n
         print 903, ( X(1+first_index+(j-1)*dist+(i-1)*stride_x), j=1, 3)
      end do

903   format(3('   (', F8.3, ',', F8.3, ')'))
      end subroutine PRINT_THREE_VECTORS_C

!---------------
      subroutine PRINT_THREE_2D_COLUMNS_C(X, m, strides)
!     input parameters
      COMPLEX(4),INTENT(IN) :: X(*)
      integer,intent(in) :: m
      integer,intent(in) :: strides(*)

!     Local parameters
      integer   i, j
      integer   first_index
      integer   stride_m
      integer   stride_n

      first_index = strides(1)
      stride_m    = strides(2)
      stride_n    = strides(3)

      do i=1, m
         print 903, ( X(1+first_index+(j-1)*stride_n+(i-1)*stride_m), j=1, 3)
      end do

903   format(3('   (', F8.3, ',', F8.3, ')'))
      end subroutine PRINT_THREE_2D_COLUMNS_C

!---------------
      subroutine PRINT_THREE_2D_COLUMNS_S(X, m, strides)
!     input parameters
      REAL(4),INTENT(IN) :: X(*)
      integer,intent(in) :: m
      integer,intent(in) :: strides(*)

!     Local parameters
      integer   i, j
      integer   first_index
      integer   stride_m
      integer   stride_n

      first_index = strides(1)
      stride_m    = strides(2)
      stride_n    = strides(3)

      do i=1, m
         print 903, ( X(1+first_index+(j-1)*stride_n+(i-1)*stride_m), j=1, 3)
      end do

903   format(3('   (', F8.3, ')'))
      end subroutine PRINT_THREE_2D_COLUMNS_S

!---------------
      subroutine PRINT_THREE_VECTORS_S(X, n, dist, strides)
!     input parameters
      REAL(4),INTENT(IN) :: X(*)
      integer,intent(in) :: n
      integer,intent(in) :: dist
      integer,intent(in) :: strides(*)

!     Local parameters
      integer   i, j
      integer   first_index
      integer   stride_x

      first_index = strides(1)
      stride_x    = strides(2)

      do i=1, n
         print 903, ( X(1+first_index+(j-1)*dist+(i-1)*stride_x), j=1, 3)
      end do

903   format(3(   F8.3  ))
      end subroutine PRINT_THREE_VECTORS_S

!---------------
      subroutine PRINT_THREE_FIRST_AND_THREE_LAST_VECTORS_2D_S(X, lengths, strides)
!     input parameters
      REAL(4),INTENT(IN) :: X(*)
      integer,intent(in) :: lengths(*)
      integer,intent(in) :: strides(*)

!     Local parameters
      real      maxerr
      integer   m, n, i, j
      integer   i0
      integer   step_1d, step_2d

      m = lengths(1)+2
      n = lengths(2)

      i0      = strides(1)
      step_1d = strides(2)
      step_2d = strides(3)

      maxerr = 0.0
      do i=1,m
         print 903, ( X(1+i0+(j-1)*step_2d+(i-1)*step_1d), j=1,3), ( X(1+i0+(n-2+j)*step_2d+(i-1)*step_1d), j=1,3)
      end do

903   format(6(   F8.3  ))
      end subroutine PRINT_THREE_FIRST_AND_THREE_LAST_VECTORS_2D_S

!---------------
      subroutine PRINT_THREE_FIRST_AND_THREE_LAST_VECTORS_2D_D(X, lengths, strides)
!     input parameters
      REAL(8),INTENT(IN) :: X(*)
      integer,intent(in) :: lengths(*)
      integer,intent(in) :: strides(*)

!     Local parameters
      real      maxerr
      integer   m, n, i, j
      integer   i0, index
      integer   step_1d, step_2d

      m = lengths(1)+2
      n = lengths(2)

      i0      = strides(1)
      step_1d = strides(2)
      step_2d = strides(3)

      maxerr = 0.0
      do i=1,m
         index = i0+i*step_1d + (j-1)*step_2d
         print 903, ( X(1+i0+(j-1)*step_2d+(i-1)*step_1d), j=1,3), ( X(1+i0+(n-2+j)*step_2d+(i-1)*step_1d), j=1,3)
      end do

903   format(6(   F8.3  ))
      end subroutine PRINT_THREE_FIRST_AND_THREE_LAST_VECTORS_2D_D

!---------------
      subroutine PRINT_THREE_FIRST_VECTORS_2D_S(X, lengths, strides)
!     input parameters
      REAL(4),INTENT(IN) :: X(*)
      integer,intent(in) :: lengths(*)
      integer,intent(in) :: strides(*)

!     Local parameters
      integer   m, n, i, j
      integer   i0
      integer   step_1d, step_2d

      m = lengths(1)
      n = lengths(2)

      i0      = strides(1)
      step_1d = strides(2)
      step_2d = strides(3)

      do i=1,m
         print 903, ( X(1+i0+(j-1)*step_2d+(i-1)*step_1d), j=1,3)
      end do

903   format(6(   F8.3  ))
      end subroutine PRINT_THREE_FIRST_VECTORS_2D_S

!---------------
      subroutine PRINT_THREE_FIRST_VECTORS_2D_D(X, lengths, strides)
!     input parameters
      REAL(8),INTENT(IN) :: X(*)
      integer,intent(in) :: lengths(*)
      integer,intent(in) :: strides(*)

!     Local parameters
      integer   m, n, i, j
      integer   i0, index
      integer   step_1d, step_2d

      m = lengths(1)
      n = lengths(2)

      i0      = strides(1)
      step_1d = strides(2)
      step_2d = strides(3)

      do i=1,m
         index = i0+i*step_1d + (j-1)*step_2d
         print 903, ( X(1+i0+(j-1)*step_2d+(i-1)*step_1d), j=1,3)
      end do

903   format(6(   F8.3  ))
      end subroutine PRINT_THREE_FIRST_VECTORS_2D_D

!---------------
      subroutine PRINT_VECTOR_Z( X, n)
!     input parameters
      COMPLEX(8),INTENT(IN) :: X(*)
      integer,intent(in) :: n

!     Local parameters
      integer   i

      do i=1,n
         print 902, i, X(i)
      end do

902   format(' X(', I4, ')= (', F8.3, ',', F8.3, ')')
      end subroutine PRINT_VECTOR_Z

!---------------
      subroutine PRINT_VECTOR_D( X, n)
!     input parameters
      REAL(8),INTENT(IN) :: X(*)
      integer,intent(in) :: n

!     Local parameters
      integer   i

      do i=1,n
         print 902, i, X(i)
      end do

902   format(' X(', I4, ')= ', F8.3)
      end subroutine PRINT_VECTOR_D

!---------------
      subroutine PRINT_THREE_VECTORS_Z(X, n, dist, strides)
!     input parameters
      COMPLEX(8),INTENT(IN) :: X(*)
      integer,intent(in) :: n
      integer,intent(in) :: dist
      integer,intent(in) :: strides(*)

!     Local parameters
      integer   i, j
      integer   first_index
      integer   stride_x

      first_index = strides(1)
      stride_x    = strides(2)

      do i=1, n
         print 903, ( X(1+first_index+(j-1)*dist+(i-1)*stride_x), j=1, 3)
      end do

903   format(3('   (', F8.3, ',', F8.3, ')'))
      end subroutine PRINT_THREE_VECTORS_Z

!---------------
      subroutine PRINT_THREE_VECTORS_D(X, n, dist, strides)
!     input parameters
      REAL(8),INTENT(IN) :: X(*)
      integer,intent(in) :: n
      integer,intent(in) :: dist
      integer,intent(in) :: strides(*)

!     Local parameters
      integer   i, j
      integer   first_index
      integer   stride_x

      first_index = strides(1)
      stride_x    = strides(2)

      do i=1, n
         print 903, ( X(1+first_index+(j-1)*dist+(i-1)*stride_x), j=1, 3)
      end do

903   format(3(   F8.3  ))
      end subroutine PRINT_THREE_VECTORS_D

!---------------
      subroutine PRINT_THREE_2D_COLUMNS_Z(X, m, strides)
!     input parameters
      COMPLEX(8),INTENT(IN) :: X(*)
      integer,intent(in) :: m
      integer,intent(in) :: strides(*)

!     Local parameters
      integer   i, j
      integer   first_index
      integer   stride_m
      integer   stride_n

      first_index = strides(1)
      stride_m    = strides(2)
      stride_n    = strides(3)

      do i=1, m
         print 903, ( X(1+first_index+(j-1)*stride_n+(i-1)*stride_m), j=1, 3)
      end do

903   format(3('   (', F8.3, ',', F8.3, ')'))
      end subroutine PRINT_THREE_2D_COLUMNS_Z

!---------------
      subroutine PRINT_THREE_2D_COLUMNS_D(X, m, strides)
!     input parameters
      REAL(8),INTENT(IN) :: X(*)
      integer,intent(in) :: m
      integer,intent(in) :: strides(*)

!     Local parameters
      integer   i, j
      integer   first_index
      integer   stride_m
      integer   stride_n

      first_index = strides(1)
      stride_m    = strides(2)
      stride_n    = strides(3)

      do i=1, m
         print 903, ( X(1+first_index+(j-1)*stride_n+(i-1)*stride_m), j=1, 3)
      end do

903   format(3('   (', F8.3, ')'))
      end subroutine PRINT_THREE_2D_COLUMNS_D

!---------------
      function CHECK_RESULT_C(X_IN, X_EXP, n)
      REAL(4) CHECK_RESULT_C
!     Input parameters
      COMPLEX(4),INTENT(IN) :: X_IN(*)
      COMPLEX(4),INTENT(IN) :: X_EXP(*)
      integer,intent(in) :: n

!     Local parameters
      complex   d
      real(4)   maxerr
      integer   i

      maxerr = 0.0
       do i=1,n
        d = X_EXP(i) - X_IN(i)
        if (ABS(d) .GT. maxerr) then
            maxerr = ABS(d)
        end if
      end do

      CHECK_RESULT_C = maxerr
      end function CHECK_RESULT_C

!---------------
      function CHECK_RESULT_S(X_IN, X_EXP, n)
      REAL(4)  CHECK_RESULT_S
!     Input parameters
      REAL(4),INTENT(IN) :: X_IN(*)
      REAL(4),INTENT(IN) :: X_EXP(*)
      integer,intent(in) :: n

!     Local parameters
      real(4)   d
      real(4)   maxerr
      integer   i

      maxerr = 0.0
       do i=1,n
        d = X_EXP(i) - X_IN(i)
        if (ABS(d) .GT. maxerr) then
            maxerr = ABS(d)
        end if
      end do

      CHECK_RESULT_S = maxerr
      end function CHECK_RESULT_S

!---------------
      function CHECK_RESULT_2D_S(X_IN, X_EXP, lengths, strides_in )
      REAL(4)  CHECK_RESULT_2D_S
!     Input parameters
      REAL(4),INTENT(IN) :: X_IN (*)
      REAL(4),INTENT(IN) :: X_EXP(*)
      integer,intent(in) :: lengths(*)
      integer,intent(in) :: strides_in(*)

!     Local parameters
      real(4)   d
      real(4)   maxerr
      integer   m, n, i, j
      integer   i0, index
      integer   step_1d, step_2d

      m = lengths(1)
      n = lengths(2)
      i0 = strides_in(1)
      step_1d = strides_in(2)
      step_2d = strides_in(3)

      maxerr = 0.0
      do j=1,n
          do i=1,m
            index = 1+i0+(i-1)*step_1d + (j-1)*step_2d
            d = X_EXP(index) - X_IN(index)
            if (ABS(d) .GT. maxerr) then
                maxerr = ABS(d)
            end if
          end do
      end do

      CHECK_RESULT_2D_S = maxerr
      end function CHECK_RESULT_2D_S

!---------------
      function CHECK_RESULT_Z(X_IN, X_EXP, n)
      REAL(8)  CHECK_RESULT_Z
!     Input parameters
      COMPLEX(8),INTENT(IN) :: X_IN(*)
      COMPLEX(8),INTENT(IN) :: X_EXP(*)
      integer,intent(in) :: n

!     Local parameters
      complex(8) d
      real(8)    maxerr
      integer    i

      maxerr = 0.0
       do i=1,n
        d = X_EXP(i) - X_IN(i)
        if (ABS(d) .GT. maxerr) then
            maxerr = ABS(d)
        end if
      end do

      CHECK_RESULT_Z = maxerr
      end function CHECK_RESULT_Z

!---------------
      function CHECK_RESULT_D(X_IN, X_EXP, n)
      REAL(8) CHECK_RESULT_D
!     Input parameters
      REAL(8),INTENT(IN) :: X_IN(*)
      REAL(8),INTENT(IN) :: X_EXP(*)
      integer,intent(in) :: n

!     Local parameters
      real(8)   d
      real(8)   maxerr
      integer   i

      maxerr = 0.0
       do i=1,n
        d = X_EXP(i) - X_IN(i)
        if (dabs(d) .GT. maxerr) then
            maxerr = dabs(d)
        end if
      end do

      CHECK_RESULT_D = maxerr
      end function CHECK_RESULT_D

!---------------
      function CHECK_RESULT_2D_D(X_IN, X_EXP, lengths, strides_in )
      REAL(8)  CHECK_RESULT_2D_D
!     Input parameters
      REAL(8),INTENT(IN) :: X_IN (*)
      REAL(8),INTENT(IN) :: X_EXP(*)
      integer,intent(in) :: lengths(*)
      integer,intent(in) :: strides_in(*)

!     Local parameters
      real(8)   d
      real(8)   maxerr
      integer   m, n, i, j
      integer   i0, index
      integer   step_1d, step_2d

      m = lengths(1)
      n = lengths(2)
      i0 = strides_in(1)
      step_1d = strides_in(2)
      step_2d = strides_in(3)

      maxerr = 1.0E-15
      do j=1,n
          do i=1,m
            index = 1+i0+(i-1)*step_1d + (j-1)*step_2d
            d = X_EXP(index) - X_IN(index)
            if (d .LT. 0.0) then
                d = -d
            end if
            if (d .GT. maxerr) then
                maxerr = d
            end if
          end do
      end do

      CHECK_RESULT_2D_D = maxerr
      end function CHECK_RESULT_2D_D
