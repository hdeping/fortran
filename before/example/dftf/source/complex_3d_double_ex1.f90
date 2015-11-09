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
!       MKL DFTI interface example program (Fortran-interface)
!
!       Forward-Backward 3D complex transform for double precision data inplace.
!
!*****************************************************************************
! Configuration parameters:
!           DFTI_FORWARD_DOMAIN = DFTI_COMPLEX  (obligatory)
!           DFTI_PRECISION      = DFTI_DOUBLE   (obligatory)
!           DFTI_DIMENSION      = 3             (obligatory)
!           DFTI_LENGTHS        = {m,n,k}       (obligatory)
!           DFTI_PLACEMENT      = DFTI_INPLACE  (default)
!           DFTI_INPUT_STRIDES  = {0, 1, M_MAX, M_MAX*N_MAX}
!                                               (default={0,1,m,m*n})
!           DFTI_FORWARD_SCALE  = 1.0           (default)
!           DFTI_BACKWARD_SCALE = 1.0/(m*n*k)   (default=1.0)
!
! Other default configuration parameters are in the mkl_dfti.f90 interface file
!*****************************************************************************

      Program COMPLEX_3D_DOUBLE_EX1

      Use MKL_DFTI
      include 'mkl_dfti_examples.fi'

      integer   m, n, k

      complex(8) :: X_IN_3D(M_MAX,N_MAX,K_MAX)
      complex(8) :: X_IN   (M_MAX*N_MAX*K_MAX)
      complex(8) :: X_EXP  (M_MAX*N_MAX*K_MAX)

      equivalence (X_IN, X_IN_3D)

      type(DFTI_DESCRIPTOR), POINTER :: hand
      integer   status
      real(8)   Scale
      integer   lengths(3)
      integer   strides_in(4)

      real(8)   maxerr
      real(8)   eps
      parameter (eps=DOUBLE_EPS)
      integer   i
      logical   failure

      failure = .FALSE.
!
!     Read input parameters from input file
!     m - size of transform  along first dimension
!     n - size of transform  along second dimension
!     k - size of transform  along third dimension
!
      read *
      read *, m
      read *, n
      read *, k

!
!     Put transform parameters
!
      lengths(1) = m
      lengths(2) = n
      lengths(3) = k

      strides_in(1) = 0
      strides_in(2) = 1
      strides_in(3) = M_MAX
      strides_in(4) = M_MAX*N_MAX

      if (LEGEND_PRINT) then
          print *
          print *, 'COMPLEX_3D_DOUBLE_EX1'
          print *, 'Forward-Backward 3D complex transform for double precision data'
          print *
          print *, 'Configuration parameters:'
          print *
          print *, 'DFTI_FORWARD_DOMAIN       = DFTI_COMPLEX'
          print *, 'DFTI_PRECISION            = DFTI_DOUBLE'
          print *, 'DFTI_DIMENSION            = 3'
          print 903, m, n, k
          print *, 'DFTI_PLACEMENT            = DFTI_INPLACE'
          print 907, (strides_in(i), i=1,4)
          print *, 'DFTI_FORWARD_SCALE        = 1.0'
          print *, 'DFTI_BACKWARD_SCALE       = 1.0/(m*n*k)'
          print *
      end if

!
!     Check input parameters
!
      if ((m*n*k) .GT. (M_MAX*N_MAX*K_MAX)) then
          print *, 'Error input parameters: (m*n*k) > (M_MAX*N_MAX*K_MAX)'
          print *, 'Please see mkl_dfti_examples.fi file'
          print *, 'TEST FAILED'
          failure = .TRUE.
          goto 101
      end if

!
!     Initialize X_IN and copy to expected X_EXP
!
      call ZERO_INIT_COMPLEX_Z(X_IN, M_MAX*N_MAX*K_MAX)
      call INIT_3D_COLUMNS_Z(X_IN, m, n, k, strides_in)
      call ZCOPY(M_MAX*N_MAX*K_MAX, X_IN, 1, X_EXP, 1)

      if (ADVANCED_DATA_PRINT) then
        print *
        do i=1,k
!           INPUT X(1:m, 1:3, i) (3D columns)
            print 911, i
            call PRINT_THREE_2D_COLUMNS_Z(X_IN_3D(1,1,i), m, strides_in)
        end do
      end if

!
!     Create DFTI descriptor
!
      Status = DftiCreateDescriptor(hand, DFTI_DOUBLE, &
                                    DFTI_COMPLEX, 3, lengths)
      if (.NOT. DftiErrorClass(Status, DFTI_NO_ERROR)) then
          call dfti_example_status_print(Status)
          print *, 'TEST FAILED : DftiCreateDescriptor(hand, DFTI_DOUBLE, ...'
          failure = .TRUE.
          goto 101
      end if

!
!     In case of data allocation in the 3D array(M_MAX*N_MAX)(Fortran interface):
!     srides(3) = M_MAX is not default parameter if m is not equal to M_MAX
!
      Status = DftiSetValue(hand, DFTI_INPUT_STRIDES, strides_in)
      if (.NOT. DftiErrorClass(Status, DFTI_NO_ERROR)) then
          call dfti_example_status_print(Status)
          print *, 'TEST FAILED : DftiSetValue(hand, DFTI_INPUT_STRIDES, strides_in)'
          failure = .TRUE.
          goto 100
      end if

!
!     Commit DFTI descriptor
!
      Status = DftiCommitDescriptor(hand)
      if (.NOT. DftiErrorClass(Status, DFTI_NO_ERROR)) then
          call dfti_example_status_print(Status)
          print *, 'TEST FAILED : DftiCommitDescriptor(hand)'
          failure = .TRUE.
          goto 100
      end if

!
!     Compute Forward transform
!
      print *
      print *, 'Compute DftiComputeForward'
      Status = DftiComputeForward(hand, X_IN)
      if (.NOT. DftiErrorClass(Status, DFTI_NO_ERROR) ) then
          call dfti_example_status_print(Status)
          print *, 'TEST FAILED : DftiComputeForward(hand, X_IN)'
          failure = .TRUE.
          goto 100
      end if

      if (ADVANCED_DATA_PRINT) then
        print *
        do i=1,k
!           Forward OUTPUT X(1:m, 1:3, i) (3D columns)
            print 912, i
            call PRINT_THREE_2D_COLUMNS_Z(X_IN_3D(1,1,i), m, strides_in)
        end do
      end if

!
!     Set Scale number for Backward transform
!
      Scale = 1.0/real(m*n*k, KIND=8)
      Status = DftiSetValue(hand, DFTI_BACKWARD_SCALE, Scale)
      if (.NOT. DftiErrorClass(Status, DFTI_NO_ERROR)) then
          call dfti_example_status_print(Status)
          print *, 'TEST FAILED : DftiSetValue(hand, DFTI_BACKWARD_SCALE, Scale)'
          failure = .TRUE.
          goto 100
      end if

!
!     Commit DFTI descriptor
!
      Status = DftiCommitDescriptor(hand)
      if (.NOT. DftiErrorClass(Status, DFTI_NO_ERROR)) then
          call dfti_example_status_print(Status)
          print *, 'TEST FAILED : DftiCommitDescriptor(hand)'
          failure = .TRUE.
          goto 100
      end if

!
!     Compute Backward transform
!
      print *
      print *, 'Compute DftiComputeBackward'
      Status = DftiComputeBackward(hand, X_IN)
      if (.NOT. DftiErrorClass(Status, DFTI_NO_ERROR) ) then
          call dfti_example_status_print(Status)
          print *, 'TEST FAILED : DftiComputeBackward(hand, X_IN)'
          failure = .TRUE.
          goto 100
      end if

      if (ADVANCED_DATA_PRINT) then
        print *
        do i=1,k
!           Backward OUTPUT X(1:m, 1:3, i) (3D columns)
            print 913, i
            call PRINT_THREE_2D_COLUMNS_Z(X_IN_3D(1,1,i), m, strides_in)
        end do
      end if

!
!     Check result
!
      maxerr = CHECK_RESULT_Z(X_IN, X_EXP, M_MAX*N_MAX*K_MAX)
      if (ACCURACY_PRINT) then
        print *
        print 904, maxerr
      end if

      if (maxerr .LT. eps) then
        print *, 'TEST PASSED'
      else
        print *, 'TEST FAILED'
        failure = .TRUE.
      end if

 100  continue

!
!     Free DFTI descriptor
!
      Status = DftiFreeDescriptor(hand)
      if (.NOT. DftiErrorClass(Status, DFTI_NO_ERROR) ) then
          call dfti_example_status_print(Status)
          print *, 'TEST FAILED : DftiFreeDescriptor(hand)'
          failure = .TRUE.
      end if

 903  format(' DFTI_LENGTHS              = {', I4,',',I4,',',I4,'}' )
 904  format(' Accuracy                  = ', G15.6)
 907  format(' DFTI_INPUT_STRIDES        = {',I4,',',I4,',',I4,',',I4,'}')
 911  format(' INPUT X(1:m, 1:3,', I2,')')
 912  format(' Forward OUTPUT X(1:m, 1:3,', I2,')')
 913  format(' Backward OUTPUT X(1:m, 1:3,', I2,')')

 101  continue

      if (failure) then
          STOP 1
      end if

      print *
      print *, 'END OF TEST'

      end
