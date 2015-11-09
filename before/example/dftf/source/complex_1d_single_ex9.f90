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
!       Backward 1D complex transform for single precision data inplace
!       Accuracy test with sinusoid tone input data
!
!*****************************************************************************
! Configuration parameters:
!           DFTI_FORWARD_DOMAIN = DFTI_COMPLEX (obligatory)
!           DFTI_PRECISION      = DFTI_SINGLE  (obligatory)
!           DFTI_DIMENSION      = 1            (obligatory)
!           DFTI_LENGTHS        = n            (obligatory)
!           DFTI_PLACEMENT      = DFTI_INPLACE (default)
!           DFTI_BACKWARD_SCALE = 1.0/n        (default=1.0)
!
! Other default configuration parameters are in the mkl_dfti.f90 interface file
!*****************************************************************************

      Program COMPLEX_1D_SINGLE_EX9

      Use MKL_DFTI
      include 'mkl_dfti_examples.fi'

      integer    n

      complex :: X_IN (M_MAX)
      complex :: X_EXP(M_MAX)

      type(DFTI_DESCRIPTOR), POINTER :: hand
      integer   status
      real      Scale

      real      maxerr
      real      eps
      parameter (eps=SINGLE_EPS)
      logical    failure

      failure = .FALSE.

!
!     Read input n - size of transform from input file
!
      read *
      read *, n

      if (LEGEND_PRINT) then
          print *
          print *, 'COMPLEX_1D_SINGLE_EX9'
          print *, 'Forward 1D complex transform for single precision data'
          print *, 'Accuracy test with sinusoid tone input data'
          print *
          print *, 'Configuration parameters:'
          print *
          print *, 'DFTI_FORWARD_DOMAIN = DFTI_COMPLEX'
          print *, 'DFTI_PRECISION      = DFTI_SINGLE'
          print *, 'DFTI_DIMENSION      = 1'
          print 903, n
          print *, 'DFTI_PLACEMENT      = DFTI_INPLACE'
          print *, 'DFTI_BACKWARD_SCALE = 1.0/n'
          print *
      end if

!
!     Check input parameters
!
      if (n .GT. M_MAX) then
          print *, 'Error input parameters n > M_MAX'
          print *, 'Please see mkl_dfti_examples.fi file'
          print *, 'TEST FAILED'
          failure = .TRUE.
          goto 101
      end if

!
!     Initialize X_IN and copy to expected X_EXP
!
      call INIT_BACKW_TONE_AND_EXP_RES_C(X_IN, X_EXP, n)

      if (ADVANCED_DATA_PRINT) then
        print *
        print *, 'INPUT vector X'
        call PRINT_VECTOR_C(X_IN, n)
      end if

!
!     Create DFTI descriptor for 1D single precision transform
!
      Status = DftiCreateDescriptor(hand, DFTI_SINGLE, &
                                    DFTI_COMPLEX, 1, n)
      if (.NOT. DftiErrorClass(Status, DFTI_NO_ERROR)) then
          call dfti_example_status_print(Status)
          print *, 'TEST FAILED : DftiCreateDescriptor(hand, DFTI_SINGLE, ...'
          failure = .TRUE.
          goto 101
      end if

!
!     Set Scale number
!
      Scale = 1.0/real(n, KIND=4)
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
!     Compute Forward transform
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
        print *, 'Backward OUTPUT vector X'
        call PRINT_VECTOR_C( X_IN, n)
      end if

!
!     Check result
!
      maxerr = CHECK_RESULT_C(X_IN, X_EXP, n)
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

 903  format(' DFTI_LENGTHS        = ', I4)
 904  format(' Accuracy            = ', G15.6)

 101  continue

      if (failure) then
          STOP 1
      end if

      print *
      print *, 'END OF TEST'

      end
