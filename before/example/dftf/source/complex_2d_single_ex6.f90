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
!       Forward-Backward 2D complex transform for single precision data not inplace.
!       Multiple 2D transform along 2-nd dimension of 3D array.
!
!*****************************************************************************
! Configuration parameters:
!
!           DFTI_FORWARD_DOMAIN = DFTI_COMPLEX     (obligatory)
!           DFTI_PRECISION      = DFTI_SINGLE      (obligatory)
!           DFTI_DIMENSION      = 2                (obligatory)
!           DFTI_LENGTHS        = {m,n}            (obligatory)
!           DFTI_PLACEMENT      = DFTI_NOT_INPLACE (default=DFTI_INPLACE)
!           DFTI_INPUT_STRIDES  = {first_index, stride_in_m, stride_in_n}
!           DFTI_OUTPUT_STRIDES = {first_index, stride_out_m, stride_out_n}
!                                                  (default={0,1,m})
!           DFTI_FORWARD_SCALE  = 1.0              (default)
!           DFTI_BACKWARD_SCALE = 1.0/(m*n)        (default=1.0)
!           DFTI_NUMBER_OF_TRANSFORMS = multiple   (default=1)
!           DFTI_INPUT_DISTANCE       = dist_in    (obligatory,
!                                              if NUMBER_OF_TRANSFORMS >1)
!           DFTI_OUTPUT_DISTANCE      = dist_out   (obligatory,
!                                              if NUMBER_OF_TRANSFORMS >1)
!
! Other default configuration parameters are in the mkl_dfti.f90 interface file
!*****************************************************************************

      Program COMPLEX_2D_SINGLE_EX6

      Use MKL_DFTI
      include 'mkl_dfti_examples.fi'

      integer    m, n
      integer    first_index
      integer    multiple

      complex :: X_IN_3D (M_MAX,N_MAX,K_MAX)
      complex :: X_IN    (M_MAX*N_MAX*K_MAX)
      complex :: X_OUT_3D(M_MAX,N_MAX,K_MAX)
      complex :: X_OUT   (M_MAX*N_MAX*K_MAX)
      complex :: X_EXP   (M_MAX*N_MAX*K_MAX)

      equivalence (X_IN, X_IN_3D)
      equivalence (X_OUT, X_OUT_3D)

      type(DFTI_DESCRIPTOR), POINTER :: hand
      integer    status
      real       Scale
      integer    lengths(2)
      integer    strides_in(3)
      integer    dist_in
      integer    strides_out(3)
      integer    dist_out

      real       maxerr
      real       eps
      parameter  (eps=SINGLE_EPS)
      integer    i
      logical    failure

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

      strides_out(1) = first_index
      strides_out(2) = 1
      strides_out(3) = M_MAX*N_MAX

      dist_out = M_MAX

      if (LEGEND_PRINT) then
          print *
          print *, 'COMPLEX_2D_SINGLE_EX6'
          print *, 'Forward-Backward 2D complex transform for single precision data'
          print *, 'Multiple 2D transform along 2-nd dimension of 3D array'
          print *
          print *, 'Configuration parameters:'
          print *
          print *, 'DFTI_FORWARD_DOMAIN       = DFTI_COMPLEX'
          print *, 'DFTI_PRECISION            = DFTI_SINGLE'
          print *, 'DFTI_DIMENSION            = 2'
          print 903, m, n
          print *, 'DFTI_PLACEMENT            = DFTI_NOT_INPLACE'
          print 905, multiple
          print 906, dist_in
          print 909, dist_out
          print 907, strides_in(1), strides_in(2), strides_in(3)
          print 910, strides_out(1), strides_out(2), strides_out(3)
          print *, 'DFTI_FORWARD_SCALE        = 1.0'
          print *, 'DFTI_BACKWARD_SCALE       = 1.0/(m*n)'
          print *
      end if

!
!     Check input parameters
!     for multiple 2D transform along the 2-nd dimension of 3D array
!
      if ((first_index + multiple*m*n) .GT. (M_MAX* N_MAX*K_MAX)) then
          print *, 'Error input parameters:'
          print *, '(first_index + multiple*m*n) > (M_MAX*N_MAX*K_MAX)'
          print *, 'Please see mkl_dfti_examples.fi file'
          print *, 'TEST FAILED'
          failure = .TRUE.
          goto 101
      end if

!
!     Initialize X_IN and copy to expected X_EXP
!
      call ZERO_INIT_COMPLEX_C(X_IN, M_MAX*N_MAX*K_MAX)
      call ZERO_INIT_COMPLEX_C(X_OUT, M_MAX*N_MAX*K_MAX)

      do i=1, multiple
          call INIT_MULTIPLE_2D_COLUMNS_C(X_IN_3D(1,i,1), m, n, strides_in)
      end do

      call CCOPY(M_MAX*N_MAX*K_MAX, X_IN, 1, X_EXP, 1)

      if (ADVANCED_DATA_PRINT) then
        print *
        print *, 'INPUT X three 2D columns'
        print *
          do i=1, multiple
             print 908, i
             call PRINT_THREE_2D_COLUMNS_C(X_IN_3D(1,i,1), m, strides_in)
          end do
      end if

!
!     Create DFTI descriptor
!
      Status = DftiCreateDescriptor(hand, DFTI_SINGLE, &
                                    DFTI_COMPLEX, 2, lengths)
      if (.NOT. DftiErrorClass(Status, DFTI_NO_ERROR)) then
          call dfti_example_status_print(Status)
          print *, 'TEST FAILED : DftiCreateDescriptor(hand, DFTI_SINGLE, ...'
          failure = .TRUE.
          goto 101
      end if

!
!     Set placement of result DFTI_NOT_INPLACE
!
      Status = DftiSetValue(hand, DFTI_PLACEMENT, &
                            DFTI_NOT_INPLACE)
      if (.NOT. DftiErrorClass(Status, DFTI_NO_ERROR)) then
          call dfti_example_status_print(Status)
          print *, 'TEST FAILED : DftiSetValue(hand, DFTI_PLACEMENT, ...'
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

          Status = DftiSetValue(hand, DFTI_OUTPUT_DISTANCE, dist_out)
          if (.NOT. DftiErrorClass(Status, DFTI_NO_ERROR)) then
              call dfti_example_status_print(Status)
              print *, 'TEST FAILED : DftiSetValue(hand, DFTI_OUTPUT_DISTANCE, dist_out)'
              failure = .TRUE.
              goto 100
          end if
      end if

!
!     Set strides parameters
!
      Status = DftiSetValue(hand, DFTI_INPUT_STRIDES, strides_in)
      if (.NOT. DftiErrorClass(Status, DFTI_NO_ERROR)) then
         call dfti_example_status_print(Status)
         print *, 'TEST FAILED : DftiSetValue(hand, DFTI_INPUT_STRIDES, strides_in)'
         failure = .TRUE.
         goto 100
      end if

      Status = DftiSetValue(hand, DFTI_OUTPUT_STRIDES, strides_out)
      if (.NOT. DftiErrorClass(Status, DFTI_NO_ERROR)) then
         call dfti_example_status_print(Status)
         print *, 'TEST FAILED : DftiSetValue(hand, DFTI_OUTPUT_STRIDES, strides_out)'
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
      Status = DftiComputeForward(hand, X_IN, X_OUT)
      if (.NOT. DftiErrorClass(Status, DFTI_NO_ERROR) ) then
          call dfti_example_status_print(Status)
          print *, 'TEST FAILED : DftiComputeForward(hand, X_IN, X_OUT)'
          failure = .TRUE.
          goto 100
      end if

      if (ADVANCED_DATA_PRINT) then
          print *
          print *, 'Forward OUTPUT X three 2D columns'
          print *
          do i=1, multiple
             print 908, i
             call PRINT_THREE_2D_COLUMNS_C(X_OUT_3D(1,i,1), m, strides_out)
          end do
      end if

!
!     Set Scale number for Backward transform
!
      Scale = 1.0/real(m*n, KIND=4)
      Status = DftiSetValue(hand, DFTI_BACKWARD_SCALE, Scale)
      if (.NOT. DftiErrorClass(Status, DFTI_NO_ERROR)) then
          call dfti_example_status_print(Status)
          print *, 'TEST FAILED : DftiSetValue(hand, DFTI_BACKWARD_SCALE, Scale)'
          failure = .TRUE.
          goto 100
      end if

!
!     Change DFTI_INPUT_DISTANCE and DFTI_OUTPUT_DISTANCE values
!     if dist_in is not equal to dist_out
!
      if (dist_in .NE. dist_out) then
          Status = DftiSetValue(hand, DFTI_INPUT_DISTANCE, dist_out)
          if (.NOT. DftiErrorClass(Status, DFTI_NO_ERROR)) then
              call dfti_example_status_print(Status)
              print *, 'TEST FAILED : DftiSetValue(hand, DFTI_INPUT_DISTANCE, dist_out)'
              failure = .TRUE.
              goto 100
          end if

          Status = DftiSetValue(hand, DFTI_OUTPUT_DISTANCE, dist_in)
          if (.NOT. DftiErrorClass(Status, DFTI_NO_ERROR)) then
              call dfti_example_status_print(Status)
              print *, 'TEST FAILED : DftiSetValue(hand, DFTI_OUTPUT_DISTANCE, dist_in)'
              failure = .TRUE.
              goto 100
          end if
      end if

!
!     Change DFTI_INPUT_STRIDES and DFTI_OUTPUT_STRIDES values
!     if strides_in is not equal to strides_out
!
      do i=1,3
        if (strides_in(i) .NE. strides_out(i)) then
            Status = DftiSetValue(hand, DFTI_OUTPUT_STRIDES, strides_in)
            if (.NOT. DftiErrorClass(Status, DFTI_NO_ERROR)) then
                call dfti_example_status_print(Status)
                print *, 'TEST FAILED : DftiSetValue(hand, DFTI_OUTPUT_STRIDES, strides_in)'
                failure = .TRUE.
                goto 100
            end if

            Status = DftiSetValue(hand, DFTI_INPUT_STRIDES, strides_out)
            if (.NOT. DftiErrorClass(Status, DFTI_NO_ERROR)) then
                call dfti_example_status_print(Status)
                print *, 'TEST FAILED : DftiSetValue(hand, DFTI_INPUT_STRIDES, strides_out)'
                failure = .TRUE.
                goto 100
            end if

            EXIT
        end if
      end do

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
      Status = DftiComputeBackward(hand, X_OUT, X_IN)
      if (.NOT. DftiErrorClass(Status, DFTI_NO_ERROR) ) then
          call dfti_example_status_print(Status)
          print *, 'TEST FAILED : DftiComputeBackward(hand, X_OUT, X_IN)'
          failure = .TRUE.
          goto 100
      end if

      if (ADVANCED_DATA_PRINT) then
          print *
          print *, 'Backward OUTPUT X three 2D columns'
          print *
          do i=1, multiple
             print 908, i
             call PRINT_THREE_2D_COLUMNS_C(X_IN_3D(1,i,1), m, strides_in)
          end do
      end if

!
!     Check result
!
      maxerr = CHECK_RESULT_C(X_IN, X_EXP, M_MAX*N_MAX*K_MAX)
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
 909  format(' DFTI_OUTPUT_DISTANCE      = ', I4)
 910  format(' DFTI_OUTPUT_STRIDES       = {',I4,',',I4,',',I4,'}')

 101  continue

      if (failure) then
          STOP 1
      end if

      print *
      print *, 'END OF TEST'

      end
