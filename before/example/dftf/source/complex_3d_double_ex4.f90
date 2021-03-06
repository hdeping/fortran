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
!       Forward-Backward 3D complex transform for double precision data not inplace.
!
!*****************************************************************************
! Configuration parameters:
!           DFTI_FORWARD_DOMAIN = DFTI_COMPLEX     (obligatory)
!           DFTI_PRECISION      = DFTI_DOUBLE      (obligatory)
!           DFTI_DIMENSION      = 3                (obligatory)
!           DFTI_LENGTHS        = {m,n,k}          (obligatory)
!           DFTI_PLACEMENT      = DFTI_NOT_INPLACE (default=DFTI_INPLACE)
!           DFTI_INPUT_STRIDES  = {first_index, step_in_m, step_in_n, step_in_k}
!           DFTI_OUTPUT_STRIDES = {first_index, step_out_m, step_out_n, step_out_k}
!                                                  (default={0,1,m,m*n})
!           DFTI_FORWARD_SCALE  = 1.0              (default)
!           DFTI_BACKWARD_SCALE = 1.0/(m*n*k)      (default=1.0)
!
! Other default configuration parameters are in the mkl_dfti.f90 interface file
!*****************************************************************************

      Program COMPLEX_3D_DOUBLE_EX4

      Use MKL_DFTI
      include 'mkl_dfti_examples.fi'

      integer   m, n, k
      integer   index_strides_in(4)
      integer   index_strides_out(4)

      complex(8) :: X_IN_3D (M_MAX,N_MAX,K_MAX)
      complex(8) :: X_IN    (M_MAX*N_MAX*K_MAX)
      complex(8) :: X_OUT_3D(M_MAX,N_MAX,K_MAX)
      complex(8) :: X_OUT   (M_MAX*N_MAX*K_MAX)
      complex(8) :: X_EXP   (M_MAX*N_MAX*K_MAX)

      equivalence (X_IN, X_IN_3D)
      equivalence (X_OUT, X_OUT_3D)

      type(DFTI_DESCRIPTOR), POINTER :: hand
      integer   status
      real(8)   Scale
      integer   lengths(3)
      integer   strides_in(4)
      integer   strides_out(4)

      integer   first_index_in
      integer   step_in_m
      integer   step_in_n
      integer   step_in_k
      integer   first_index_out
      integer   step_out_m
      integer   step_out_n
      integer   step_out_k

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
!     index_strides_in[0] - displacement from the first element of input array
!     index_strides_in[1] - index step along first dimension of input array
!     index_strides_in[2] - index step along second dimension of input array
!     index_strides_in[3] - index step along third dimension of input array
!     index_strides_out[0] - displacement from the first element of output array
!     index_strides_out[1] - index step along first dimension of output array
!     index_strides_out[2] - index step along second dimension of output array
!     index_strides_out[3] - index step along third dimension of output array
!
      read *
      read *, m
      read *, n
      read *, k
      read *, index_strides_in(1)
      read *, index_strides_in(2)
      read *, index_strides_in(3)
      read *, index_strides_in(4)
      read *, index_strides_out(1)
      read *, index_strides_out(2)
      read *, index_strides_out(3)
      read *, index_strides_out(4)

!
!     Put transform parameters
!
      lengths(1) = m
      lengths(2) = n
      lengths(3) = k

      first_index_in  = index_strides_in(1);
      step_in_m       = index_strides_in(2);
      step_in_n       = index_strides_in(3);
      step_in_k       = index_strides_in(4);

      strides_in(1) = first_index_in
      strides_in(2) = step_in_m
      strides_in(3) = step_in_n * M_MAX
      strides_in(4) = step_in_k * M_MAX*N_MAX

      first_index_out  = index_strides_out(1);
      step_out_m       = index_strides_out(2);
      step_out_n       = index_strides_out(3);
      step_out_k       = index_strides_out(4);

      strides_out(1) = first_index_out
      strides_out(2) = step_out_m
      strides_out(3) = step_out_n * M_MAX
      strides_out(4) = step_out_k * M_MAX*N_MAX

      if (LEGEND_PRINT) then
          print *
          print *, 'COMPLEX_3D_DOUBLE_EX4'
          print *, 'Forward-Backward 3D complex transform for double precision data'
          print *
          print *, 'Configuration parameters:'
          print *
          print *, 'DFTI_FORWARD_DOMAIN       = DFTI_COMPLEX'
          print *, 'DFTI_PRECISION            = DFTI_DOUBLE'
          print *, 'DFTI_DIMENSION            = 3'
          print 903, m, n, k
          print *, 'DFTI_PLACEMENT            = DFTI_NOT_INPLACE'
          print 907, (strides_in(i), i=1,4)
          print 908, (strides_out(i), i=1,4)
          print *, 'DFTI_FORWARD_SCALE        = 1.0'
          print *, 'DFTI_BACKWARD_SCALE       = 1.0/(m*n*k)'
          print *
      end if

!
!     Check test input parameters
!
      if ((first_index_in+m*step_in_m .GT. M_MAX) .OR. &
          (first_index_in+n*step_in_n .GT. N_MAX) .OR. &
          (first_index_in+k*step_in_k .GT. K_MAX) .OR. &
          (first_index_out+m*step_out_m .GT. M_MAX) .OR. &
          (first_index_out+n*step_out_n .GT. N_MAX) .OR. &
          (first_index_out+k*step_out_k .GT. K_MAX)) then
          print *, 'Error input parameters'
          print *, 'TEST FAILED'
          failure = .TRUE.
          goto 101
      end if

!
!     Initialize X_IN and copy to expected X_EXP
!
      call ZERO_INIT_COMPLEX_Z(X_IN, M_MAX*N_MAX*K_MAX)
      call ZERO_INIT_COMPLEX_Z(X_OUT, M_MAX*N_MAX*K_MAX)
      call INIT_3D_COLUMNS_Z(X_IN, m, n, k, strides_in)
      call ZCOPY(M_MAX*N_MAX*K_MAX, X_IN, 1, X_EXP, 1)

      if (ADVANCED_DATA_PRINT) then
        print *
        do i=1,k*step_in_k, step_in_k
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
!     Set placement of result DFTI_NOT_INPLACE
!
      Status = DftiSetValue(hand, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
      if (.NOT. DftiErrorClass(Status, DFTI_NO_ERROR)) then
          call dfti_example_status_print(Status)
          print *, 'TEST FAILED : DftiSetValue(hand, DFTI_PLACEMENT, DFTI_NOT_INPLACE)'
          failure = .TRUE.
          goto 100
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
        do i=1,k*step_in_k, step_in_k
!           Forward OUTPUT X(1:m, 1:3, i) (3D columns)
            print 912, i
            call PRINT_THREE_2D_COLUMNS_Z(X_OUT_3D(1,1,i), m, strides_out)
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
!     Change DFTI_INPUT_STRIDES and DFTI_OUTPUT_STRIDES values
!     if strides_in is not equal to strides_out
!
      do i=1,4
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
        do i=1,k*step_in_k, step_in_k
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
 908  format(' DFTI_OUTPUT_STRIDES       = {',I4,',',I4,',',I4,',',I4,'}')
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
