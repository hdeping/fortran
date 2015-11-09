!*****************************************************************************
!                              INTEL CONFIDENTIAL
! Copyright(C) 2006-2009 Intel Corporation. All Rights Reserved.
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
!       Forward-Backward 2D real transform for double precision data inplace.
!       Multiple 2D transform along 2-nd dimension of 3D array.
!
!*****************************************************************************
! Configuration parameters:
!           DFTI_FORWARD_DOMAIN = DFTI_REAL        (obligatory)
!           DFTI_PRECISION      = DFTI_DOUBLE      (obligatory)
!           DFTI_DIMENSION      = 2                (obligatory)
!           DFTI_LENGTHS        = {m,n}            (obligatory)
!           DFTI_PACKED_FORMAT  = DFTI_PACK_FORMAT (default=DFTI_CCS_FORMAT)
!           DFTI_PLACEMENT      = DFTI_INPLACE     (default)
!           DFTI_INPUT_STRIDES  = {first_index, stride_in_m, stride_in_n}
!                                                  (default={0,1,m})
!           DFTI_FORWARD_SCALE  = 1.0              (default)
!           DFTI_BACKWARD_SCALE = 1.0/(m*n)        (default=1.0)
!           DFTI_NUMBER_OF_TRANSFORMS = multiple   (default=1)
!           DFTI_INPUT_DISTANCE       = dist_in    (obligatory,
!                                              if NUMBER_OF_TRANSFORMS >1)
!
! Other default configuration parameters are in the mkl_dfti.f90 interface file
!*****************************************************************************

      Program REAL_2D_PACK_DOUBLE_EX5

      Use MKL_DFTI
      include 'mkl_dfti_examples.fi'

      integer   m, n

      real(8) :: X_IN_3D(M_MAX,N_MAX,K_MAX)
      real(8) :: X_IN   (M_MAX*N_MAX*K_MAX)
      real(8) :: X_EXP  (M_MAX*N_MAX*K_MAX)
      real(8) :: X_EXP_3D(M_MAX,N_MAX,K_MAX)

      equivalence (X_EXP, X_EXP_3D)
      equivalence (X_IN, X_IN_3D)

      type(DFTI_DESCRIPTOR), POINTER :: hand
      integer   status
      real(8)   Scale
      integer   lengths(2)
      integer   strides_in(3)

      integer   first_index
      integer   multiple
      integer   dist_in

      real(8)   maxerr, err
      real(8)   eps
      parameter (eps=DOUBLE_EPS)
      integer   i
      logical   failure

      failure = .FALSE.

!
!     Read input parameters from input file
!     m - size of transform  along first dimension
!     n - size of transform  along second dimension
!     first_index - displacement from the first element of data array
!     multiple - number of multiple transform
!
      read *
      read *, m
      read *, n
      read *, first_index
      read *, multiple

!
!     Put transform parameters
!
      lengths(1) = m
      lengths(2) = n

      strides_in(1) = first_index
      strides_in(2) = 1
      strides_in(3) = M_MAX*N_MAX

      dist_in = M_MAX

      if (LEGEND_PRINT) then
          print *
          print *, 'REAL_2D_PACK_DOUBLE_EX5'
          print *, 'Forward-Backward 2D real transform for double precision data'
          print *, 'Multiple 2D transform along 2-nd dimension of 3D array'
          print *
          print *, 'Configuration parameters:'
          print *
          print *, 'DFTI_FORWARD_DOMAIN       = DFTI_REAL'
          print *, 'DFTI_PRECISION            = DFTI_DOUBLE'
          print *, 'DFTI_DIMENSION            = 2'
          print 903, m, n
          print 905, multiple
          print 906, dist_in
          print *, 'DFTI_PACKED_FORMAT        = DFTI_PACK_FORMAT'
          print *, 'DFTI_PLACEMENT            = DFTI_INPLACE'
          print 907, (strides_in(i), i=1,3)
          print *, 'DFTI_FORWARD_SCALE        = 1.0'
          print *, 'DFTI_BACKWARD_SCALE       = 1.0/(m*n)'
          print *
      end if

!
!     Check input parameters
!
      if ((m .GT. M_MAX) .OR. (n .GT. K_MAX)) then
          print *, 'Error input parameters: (m .GT. M_MAX) .OR. (n .GT. K_MAX)'
          print *, 'Please see mkl_dfti_examples.fi file'
          print *, 'TEST FAILED'
          failure = .TRUE.
          goto 101
      end if

      if ((first_index+m*multiple) .GT. (M_MAX*N_MAX)) then
          print *, 'Error input parameter first_index'
          print *, 'Please see mkl_dfti_examples.fi file'
          print *, 'TEST FAILED'
          failure = .TRUE.
          goto 101
      end if

!
!     Initialize X_IN and copy to expected X_EXP
!
      call ZERO_INIT_REAL_D(X_IN, M_MAX*N_MAX*K_MAX)

      do i=1, multiple
        call INIT_MULTIPLE_2D_COLUMNS_D(X_IN_3D(1,i,1), m, n, strides_in)
      end do

      call DCOPY(M_MAX*N_MAX*K_MAX, X_IN, 1, X_EXP, 1)

      if (ADVANCED_DATA_PRINT) then
        print *
        print *, 'INPUT vector X (2D columns)'
          do i=1, multiple
              print *
              print 908, i
              call PRINT_THREE_FIRST_VECTORS_2D_D(X_IN_3D(1,i,1), lengths, strides_in)
          end do
      end if

!
!     Create DFTI descriptor
!
      Status = DftiCreateDescriptor(hand, DFTI_DOUBLE, &
                                    DFTI_REAL, 2, lengths)
      if (.NOT. DftiErrorClass(Status, DFTI_NO_ERROR)) then
          call dfti_example_status_print(Status)
          print *, 'TEST FAILED : DftiCreateDescriptor(hand, DFTI_DOUBLE, ...'
          failure = .TRUE.
          goto 101
      end if

!
!     Set packed format for output complex conjugate-symmetric data
!
      Status = DftiSetValue(hand, DFTI_PACKED_FORMAT, DFTI_PACK_FORMAT)
      if (.NOT. DftiErrorClass(Status, DFTI_NO_ERROR)) then
          call dfti_example_status_print(Status)
          print *, 'TEST FAILED : DftiSetValue(hand, DFTI_PACKED_FORMAT, DFTI_PACK_FORMAT)'
          failure = .TRUE.
          goto 100
      end if

!
!     Set parameters for multiple transform mode
!
      if ( multiple .GT. 1 ) then
          Status = DftiSetValue(hand, DFTI_NUMBER_OF_TRANSFORMS, multiple)
          if (.NOT. DftiErrorClass(Status, DFTI_NO_ERROR)) then
              call dfti_example_status_print(Status)
              print *, 'TEST FAILED : DftiSetValue(hand, DFTI_NUMBER_OF_TRANSFORMS, multiple)'
              failure = .TRUE.
              goto 100
          end if

          Status = DftiSetValue(hand, DFTI_INPUT_DISTANCE, dist_in)
          if (.NOT. DftiErrorClass(Status, DFTI_NO_ERROR)) then
              call dfti_example_status_print(Status)
              print *, 'TEST FAILED : DftiSetValue(hand, DFTI_INPUT_DISTANCE, dist_in)'
              failure = .TRUE.
              goto 100
          end if
      end if

!
!     Set strides
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
        print *, 'Forward OUTPUT vector X (2D columns)'
          do i=1, multiple
              print *
              print 908, i
              call PRINT_THREE_FIRST_VECTORS_2D_D(X_IN_3D(1,i,1), lengths, strides_in)
          end do
      end if

!
!     Set Scale number for Backward transform
!
      Scale = 1.0/real(m*n, KIND=8)
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
        print *, 'Backward OUTPUT vector X (2D columns)'
          do i=1, multiple
              print *
              print 908, i
              call PRINT_THREE_FIRST_VECTORS_2D_D(X_IN_3D(1,i,1), lengths, strides_in)
          end do
      end if

!
!     Check result
!
      maxerr = 0.0
      do i=1, multiple
          err = CHECK_RESULT_2D_D(X_IN_3D(1,i,1), X_EXP_3D(1,i,1), lengths, strides_in)
          maxerr = max(err, maxerr)
      end do

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

 903  format(' DFTI_LENGTHS              = {', I4,',',I4,'}' )
 904  format(' Accuracy                  = ', G15.6)
 905  format(' DFTI_NUMBER_OF_TRANSFORMS = ', I4)
 906  format(' DFTI_INPUT_DISTANCE       = ', I4)
 907  format(' DFTI_INPUT_STRIDES        = {',I4,',',I4,',',I4,'}')
 908  format(' Transform number          = ', I2)

 101  continue

      if (failure) then
          STOP 1
      end if

      print *
      print *, 'END OF TEST'

      end
