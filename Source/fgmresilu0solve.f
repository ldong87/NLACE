!*******************************************************************************
!                              INTEL CONFIDENTIAL
!   Copyright(C) 2006-2010 Intel Corporation. All Rights Reserved.
!   The source code contained  or  described herein and all documents related to
!   the source code ("Material") are owned by Intel Corporation or its suppliers
!   or licensors.  Title to the  Material remains with  Intel Corporation or its
!   suppliers and licensors. The Material contains trade secrets and proprietary
!   and  confidential  information of  Intel or its suppliers and licensors. The
!   Material  is  protected  by  worldwide  copyright  and trade secret laws and
!   treaty  provisions. No part of the Material may be used, copied, reproduced,
!   modified, published, uploaded, posted, transmitted, distributed or disclosed
!   in any way without Intel's prior express written permission.
!   No license  under any  patent, copyright, trade secret or other intellectual
!   property right is granted to or conferred upon you by disclosure or delivery
!   of the Materials,  either expressly, by implication, inducement, estoppel or
!   otherwise.  Any  license  under  such  intellectual property  rights must be
!   express and approved by Intel in writing.
!*******************************************************************************
!  Content:
!  Intel MKL example of RCI Flexible Generalized Minimal RESidual method with
!  ILU0 Preconditioner
!*******************************************************************************

!---------------------------------------------------------------------------
!  Example program for solving non-degenerate system of equations.
!  Full functionality of RCI FGMRES solver is exploited. Example shows how
!  ILU0 preconditioner accelerates the solver by reducing the number of
!  iterations.
!---------------------------------------------------------------------------

c      PROGRAM FGMRES_FULL_FUNCT_F (from example in the MKL library)
!     modified by JFD - 2010-12-21
      subroutine fgmresilu0solve 
      USE IOUNIT!access the IO integer for screen-printing
      USE MAINMEM! access neqns, atang, jdiag, jcsr, aload, tvect, nnz, tol

      IMPLICIT NONE

      INCLUDE "mkl_rci.fi"

c      INTEGER neqns!N
c      PARAMETER(neqns=4)
      INTEGER SIZE! size of the IPAR and DPAR vectors for commanding the GRMES
      PARAMETER (SIZE=128)! fixed by Intel
!---------------------------------------------------------------------------
! Define arrays for the upper triangle of the coefficient matrix
! Compressed sparse row storage is used for sparse representation
!---------------------------------------------------------------------------
c      INTEGER jdiag(5)!IA
c      DATA jdiag /1,4,7,10,13/
c      INTEGER jcsr(12)!JA
c      DATA jcsr /1,2,3,1,2,4,1,3,4,2,3,4/
c      DOUBLE PRECISION atang(12)!A
c      DATA atang / 4.,-1.,-1.,-1.,4.,-1.,-1.,4.,-1.,-1.,-1.,4./
      DOUBLE PRECISION BILU0(nnz),TRVEC(neqns)
!---------------------------------------------------------------------------
! Allocate storage for the ?par parameters and the solution/rhs/residual vectors
!---------------------------------------------------------------------------
      INTEGER IPAR(SIZE),IERR
      DOUBLE PRECISION DPAR(SIZE)
      DOUBLE PRECISION TMP(neqns*(2*neqns+1)+(neqns*(neqns+9))/2+1)
c      DOUBLE PRECISION EXPECTED_SOLUTION(neqns)
c      DATA EXPECTED_SOLUTION /1.0,1.0,1.0,1.0/
c      DOUBLE PRECISION aload(neqns)!RHS
      DOUBLE PRECISION B(neqns)
c      DOUBLE PRECISION tvect(neqns)! COMPUTED_SOLUTION
      DOUBLE PRECISION RESIDUAL(neqns)

      INTEGER INCX
c      INTEGER REF_NIT
c      DOUBLE PRECISION REF_NORM2
      DOUBLE PRECISION NRM2
      PARAMETER (INCX=1)
c      PARAMETER (REF_NIT=2,REF_NORM2=7.7723d0)
!---------------------------------------------------------------------------
! Some additional variables to use with the RCI (P)FGMRES solver
!---------------------------------------------------------------------------
      INTEGER ITERCOUNT
      INTEGER RCI_REQUEST, I
      DOUBLE PRECISION DVAR
!---------------------------------------------------------------------------
! An external BLAS function is taken from MKL BLAS to use
! with the RCI (P)FGMRES solver
!---------------------------------------------------------------------------
      DOUBLE PRECISION DNRM2
      EXTERNAL DNRM2

c      WRITE( *,'(A,A)') '---------------------------------------------',
c     1 '----------------------'
c      WRITE(*,'(A,A)') 'The FULLY ADVANCED example RCI FGMRES with',
c     1 ' ILU0 preconditioner'
c      WRITE(*,'(A,A)') 'to solve the non-degenerate algebraic system',
c     1 ' of linear equations'
c      WRITE( *,'(A,A)') '---------------------------------------------',
c     1 '----------------------'
!---------------------------------------------------------------------------
! Initialize variables and the right hand side through matrix-vector product
!---------------------------------------------------------------------------
c      CALL MKL_DCSRGEMV('N',neqns,atang,jdiag,jcsr,
c     &                  EXPECTED_SOLUTION,aload)
!---------------------------------------------------------------------------
! Save the right-hand side in vector B for future use
!---------------------------------------------------------------------------
      CALL DCOPY(neqns, aload, 1, B, 1)
!---------------------------------------------------------------------------
! Initialize the initial guess
!---------------------------------------------------------------------------
      tvect(1)=0.00001d0
      DO I=2,neqns
        tvect(I)=0.0d0
      ENDDO

!---------------------------------------------------------------------------
! Initialize the solver
!---------------------------------------------------------------------------
      CALL DFGMRES_INIT(neqns, tvect, aload, 
     &                  RCI_REQUEST, IPAR, DPAR, TMP)
      IF (RCI_REQUEST.NE.0) GOTO 999

!---------------------------------------------------------------------------
! Calculate ILU0 preconditioner.
!                      !ATTENTION!
! DCSRILU0 routine uses some IPAR, DPAR set by DFGMRES_INIT routine.
! Important for DCSRILU0 default entries set by DFGMRES_INIT are
! ipar(2) = 6 - output of error messages to the screen,
! ipar(6) = 1 - allow output of error messages,
! ipar(31)= 0 - abort DCSRILU0 calculations if routine meets zero diagonal element.
!
! If ILU0 is going to be used out of MKL FGMRES context, than the values
! of ipar(2), ipar(6), ipar(31), dpar(31), and dpar(32) should be user
! provided before the DCSRILU0 routine call.
!
! In this example, specific for DCSRILU0 entries are set in turn:
! ipar(31)= 1 - change small diagonal value to that given by dpar(32),
! dpar(31)= 1.D-20 instead of the default value set by DFGMRES_INIT.
!                  It is a small value to compare a diagonal entry with it.
! dpar(32)= 1.D-16 instead of the default value set by DFGMRES_INIT.
!                  It is the target value of the diagonal value if it is
!                  small as compared to dpar(31) and the routine should change
!                  it rather than abort DCSRILU0 calculations.
!---------------------------------------------------------------------------

        IPAR(31)=1
        DPAR(31)=1.D-20
        DPAR(32)=1.D-16
        CALL DCSRILU0(neqns,atang,jdiag,jcsr,BILU0,IPAR,DPAR,IERR)
        NRM2=DNRM2(nnz, BILU0, INCX)

        IF(IERR.ne.0) THEN
          WRITE(*,'(A,A,I1)') ' Error after calculation of the',
     1    ' preconditioner DCSRILU0',IERR
          GOTO 998
        ENDIF

!---------------------------------------------------------------------------
! Set the desired parameters: (see the Intel MKL website, section solvers)
!---------------------------------------------------------------------------
c      IPAR(1)=nnz! size of the problem; set by DFGMRES_INIT
c      IPAR(2)=6! all error and warning messages are printed on the screen
c      IPAR(3)=1! current stage of RCI FGRMES algorithm
c      IPAR(4)=0! current iteration number
      IPAR(5)=1200! maximum number of iterations
c      IPAR(6)=1! Rule for error printing; set by DFGMRES_INIT
c      IPAR(7)=1! Rule for warning printing; set by DFGMRES_INIT
      IPAR(8)=1! >0: stop after IPAR(5) iterations [namely IPAR(4)>IPAR(5)]; =0: do not stop
c      IPAR(9)=0! >0: perform the residual stopping test; =0: do not perform
c      IPAR(10)=1! >0: perform the user-defined stopping test [RCI_REQUEST=2]; =0: do not perform
      IPAR(11)=1! >0: run the preconditioned GMRES [RCI_REQUEST=3]; =0: do not use preconditioning
c      IPAR(12)=0! >0: test the vector norm==0; =0: do not test
c      IPAR(13)=0! =0: DGMRES_GET updates the solution; >0: DGMRES_GET copies sol to rhs; <0: DGMRES returns the iteration number (no update)
c      IPAR(14)=0! number of iterations before the restart
      IPAR(15)=250! restart after IPAR(15) iterations; IPAR(15) is also the size of the Krylov space
c      IPAR(16:SIZE)! service or unused variables in the current version of MKL
      DPAR(1)=1.0d-6! relative tolerance
c      DPAR(2)=0.0d-0! absolute tolerance
c      DPAR(3)=0.0d-0! euclidian norm of initial residual
c      DPAR(4)=0.0d-0! service variable
c      DPAR(5)=0.0d-0! euclidian norm of the current residual
c      DPAR(6)=0.0d-0! euclidian norm of the previous iteration residual
c      DPAR(7)=0.0d-0! norm of the generated vector
c      DPAR(8)=1.0d-12! tolerance for the zero norm vector
c      DPAR(9:SIZE)! service of unused variables in the current version of MKL

!---------------------------------------------------------------------------
! Check the correctness and consistency of the newly set parameters
!---------------------------------------------------------------------------
      CALL DFGMRES_CHECK(neqns, tvect, aload, RCI_REQUEST,
     1 IPAR, DPAR, TMP)
      IF (RCI_REQUEST.NE.0) GOTO 999
!---------------------------------------------------------------------------
! Print the info about the RCI FGMRES method
!---------------------------------------------------------------------------
c      WRITE( *,'(A)') ' '
c      WRITE( *,'(A,A)') 'Some info about the current run of RCI FGMRES',
c     1 ' method:'
      if (countFGcalls.eq.1) then! to avoid printing all the time
      WRITE( *,'(A,A)') '        Options of RCI FGMRES with ILU0'
c      WRITE( *,'(A)') ' '
      IF (IPAR(8).NE.0) THEN
        WRITE(*,'(A,I1,A,A)') '        IPAR(8)=',IPAR(8),', performing',
     &  ' the automatic test for the max number of iter.'
      ELSE
        WRITE(*,'(A,I1,A,A)') '        IPAR(8)=',IPAR(8),', skipping',
     &  ' the automatic test for the max number of iter.'
      ENDIF
c      WRITE( *,'(A)') '+++'
      IF (IPAR(9).NE.0) THEN
        WRITE(*,'(A,I1,A,A)') '        IPAR(9)=',IPAR(9),', performing',
     &  ' the automatic residual test'
      ELSE
        WRITE(*,'(A,I1,A,A)') '        IPAR(9)=',IPAR(9),', skipping',
     &  ' the automatic residual test'
      ENDIF
c      WRITE( *,'(A)') '+++'
      IF (IPAR(10).NE.0) THEN
        WRITE(*,'(A,I1,A,A)') '        IPAR(10)=',IPAR(10),', the',
     1 ' user-defined stopping test will be requested via RCI_REQUEST=2'
      ELSE
        WRITE(*,'(A,I1,A,A,A)') '        IPAR(10)=',IPAR(10),', the',
     1  ' user-defined stopping test will not be requested, thus,',
     1  ' RCI_REQUEST will not take the value 2'
      ENDIF
c      WRITE( *,'(A)') '+++'
      IF (IPAR(11).NE.0) THEN
        WRITE(*,'(A,I1,A,A)') '        IPAR(11)=',IPAR(11),', the',
     1  ' Preconditioned FGMRES iterations will be performed, thus,'
        WRITE(*,'(A,A)') '          the preconditioner action will',
     &  ' be requested via RCI_REQUEST=3'
      ELSE
        WRITE(*,'(A,I1,A,A)') '        IPAR(11)=',IPAR(11),', the',
     1  ' Preconditioned FGMRES iterations will not be performed,'
        WRITE(*,'(A,A)') '          thus, RCI_REQUEST will not take',
     &  ' the value 3'
      ENDIF
c      WRITE( *,'(A)') '+++'
      IF (IPAR(12).NE.0) THEN
        WRITE(*,'(A,I1,A,A)')'        IPAR(12)=',IPAR(12),', the ',
     &  'automatic test for the norm of the next generated vector'
        WRITE(*,'(A,A)') '          is not equal to zero up to ',
     &  'rounding and computational errors will be performed,'
        WRITE(*,'(A,A)') '          thus, RCI_REQUEST will not',
     &  ' take the value 4'
      ELSE
        WRITE(*,'(A,I1,A,A)')'        IPAR(12)=',IPAR(12),', the',
     &  ' automatic test for the norm of the next generated vector'
        WRITE(*,'(A,A)') '          is not equal to zero up to ',
     &  'rounding and computational errors will be skipped,'
        WRITE(*,'(A,A)') '          thus, the user-defined test',
     &  ' will be requested via RCI_REQUEST=4'
      ENDIF
      endif
c      WRITE( *,'(A)') '+++'
!---------------------------------------------------------------------------
! Compute the solution by RCI (P)FGMRES solver with preconditioning
! Reverse Communication starts here
!---------------------------------------------------------------------------
1     CALL DFGMRES(neqns, tvect, aload, RCI_REQUEST, IPAR,
     1 DPAR, TMP)
!---------------------------------------------------------------------------
! If RCI_REQUEST=0, then the solution was found with the required precision
!---------------------------------------------------------------------------
      IF (RCI_REQUEST.EQ.0) GOTO 3
!---------------------------------------------------------------------------
! If RCI_REQUEST=1, then compute the vector A*TMP(IPAR(22))
! and put the result in vector TMP(IPAR(23))
!---------------------------------------------------------------------------
      IF (RCI_REQUEST.EQ.1) THEN
        CALL MKL_DCSRGEMV('N',neqns,atang,jdiag,jcsr, 
     &                    TMP(IPAR(22)), TMP(IPAR(23)))
        GOTO 1
      ENDIF
!---------------------------------------------------------------------------
! If RCI_request=2, then do the user-defined stopping test
! The residual stopping test for the computed solution is performed here
!---------------------------------------------------------------------------
! NOTE: from this point vector B(N) is no longer containing the right-hand
! side of the problem! It contains the current FGMRES approximation to the
! solution. If you need to keep the right-hand side, save it in some other
! vector before the call to DFGMRES routine. Here we saved it in vector
! aload(N). The vector B is used instead of aload to preserve the original
! right-hand side of the problem and guarantee the proper restart of FGMRES
! method. Vector B will be altered when computing the residual stopping
! criterion!
!---------------------------------------------------------------------------
      IF (RCI_REQUEST.EQ.2) THEN
! Request to the DFGMRES_GET routine to put the solution into B(N) via IPAR(13)
        IPAR(13)=1
! Get the current FGMRES solution in the vector B(N)
        CALL DFGMRES_GET(neqns, tvect,B,RCI_REQUEST, IPAR,
     1 DPAR, TMP, ITERCOUNT)
! Compute the current true residual via MKL (Sparse) BLAS routines
        CALL MKL_DCSRGEMV('N',neqns,atang,jdiag,jcsr,B,RESIDUAL)
        CALL DAXPY(neqns, -1.0D0, aload, 1, RESIDUAL, 1)
        DVAR=DNRM2(neqns, RESIDUAL, 1)
c        write(*,*)'the computed residual is: ',DVAR! added JFD
        if (IPAR(4).eq.1) then
          write(*,'(A,A,1p,E10.4)') '        ILU0+FGMRES: at iteration',
     &       '  1,    the (initial) residual: ',DVAR
        endif
        IF (DVAR.LT.(5.0d-2*tol)) THEN
          GOTO 3
        ELSE
          GOTO 1
        ENDIF
      ENDIF
!---------------------------------------------------------------------------
! If RCI_REQUEST=3, then apply the preconditioner on the vector
! TMP(IPAR(22)) and put the result in vector TMP(IPAR(23))
! Here is the recommended usage of the result produced by ILU0 routine
! via standard MKL Sparse Blas solver routine mkl_dcsrtrsv.
!---------------------------------------------------------------------------
      IF (RCI_REQUEST.EQ.3) THEN
        CALL MKL_DCSRTRSV('L','N','U',neqns,BILU0,jdiag,
     &                    jcsr,TMP(IPAR(22)),TRVEC)
        CALL MKL_DCSRTRSV('U','N','N',neqns,BILU0,jdiag,
     &                    jcsr,TRVEC,TMP(IPAR(23)))
       GOTO 1
      ENDIF
!---------------------------------------------------------------------------
! If RCI_REQUEST=4, then check if the norm of the next generated vector is
! not zero up to rounding and computational errors. The norm is contained
! in DPAR(7) parameter
!---------------------------------------------------------------------------
      IF (RCI_REQUEST.EQ.4) THEN
        IF (DPAR(7).LT.(5.0d-2*tol)) THEN
          GOTO 3
        ELSE
          GOTO 1
        ENDIF
!---------------------------------------------------------------------------
! If RCI_REQUEST=anything else, then DFGMRES subroutine failed
! to compute the solution vector: tvect(N)
!---------------------------------------------------------------------------
      ELSE
        GOTO 999
      ENDIF
!---------------------------------------------------------------------------
! Reverse Communication ends here
! Get the current iteration number and the FGMRES solution. (DO NOT FORGET to
! call DFGMRES_GET routine as computed_solution is still containing
! the initial guess!). Request to DFGMRES_GET to put the solution into
! vector tvect(N) via IPAR(13)
!---------------------------------------------------------------------------
3     IPAR(13)=0
      CALL DFGMRES_GET(neqns, tvect, aload, RCI_REQUEST,
     &        IPAR, DPAR, TMP, ITERCOUNT)
      DO i = 1, neqns
         aload(i) = tvect(i)
c         WRITE(16,*) ' aload(',i,') = ', aload(i)
      END DO
!---------------------------------------------------------------------------
! Print solution vector: tvect(N) and
! the number of iterations: ITERCOUNT
!---------------------------------------------------------------------------
c      WRITE( *,'(A)') ' '
c      WRITE( *,'(A)') 'The system has been solved'
c      WRITE( *,'(A)') ' '
c      WRITE( *,'(A)') 'The following solution has been obtained:'
c      DO I=1,neqns
c         WRITE(*,'(A18,I1,A2,E10.3)') 'COMPUTED_SOLUTION(',I,')=',
c     1 tvect(I)
c      ENDDO
c      WRITE( *,'(A)') ' '
c      WRITE( *,'(A)') 'The expected solution is:'
c      DO I=1,neqns
c         WRITE(*,'(A18,I1,A2,E10.3)') 'EXPECTED_SOLUTION(',I,')=',
c     1 EXPECTED_SOLUTION(I)
c      ENDDO
c      WRITE( *,'(A)') ' '
      WRITE(*,'(A,A,I6,A,1p,E10.4)')'        ILU0+FGMRES: # of ',
     &   'iterations: ',ITERCOUNT,' and final residual: ',DVAR
c      WRITE( *,'(A)') ' '

!---------------------------------------------------------------------------
! Release internal MKL memory that might be used for computations
! NOTE: It is important to call the routine below to avoid memory leaks
! unless you disable MKL Memory Manager
!---------------------------------------------------------------------------
      CALL MKL_FREEBUFFERS
      return

c      IF(ITERCOUNT.EQ.REF_NIT.AND.DABS(REF_NORM2-NRM2).LT.1.D-6) THEN
c         WRITE( *,'(A)') ' '
cc         WRITE( *,'(A)') '---------------------------------------------'
c         WRITE( *,'(A,A)') 'Fortran example of FGMRES with ILU0',
c     1 ' preconditioner '
c         WRITE( *,'(A,A)') 'has successfully PASSED all stages of',
c     1 ' computations'
cc         WRITE( *,'(A)') '---------------------------------------------'
cc         WRITE( *,'(A)') ' '
cc         STOP 0
c      ELSE
c         WRITE( *,'(A,A)') 'Probably, the preconditioner was computed',
c     1 ' incorrectly:'
c         WRITE( *,'(A,F9.6,A,F9.6)')
c     1 'Either the preconditioner norm',NRM2,
c     2 ' differs from the expected norm',REF_NORM2
c         WRITE( *,'(A,I2,A,I2)'),
c     1 'and/or the number of iterations ', ITERCOUNT,
c     2 ' differs from the expected number ', REF_NIT
c         WRITE( *,'(A)') ' '
c         WRITE( *,'(A,A)') '------------------------------------------',
c     1 '----------------------'
c         WRITE( *,'(A,A)') 'Unfortunately, FGMRES+ILU0 Fortran example',
c     1 ' has FAILED'
c         WRITE( *,'(A,A)') '------------------------------------------',
c     1 '----------------------'
c         WRITE( *,'(A)') ' '
c         STOP 1
c      END IF
999   WRITE( *,'(A,I2)') 'The solver has returned the ERROR code ',
     1 RCI_REQUEST
!---------------------------------------------------------------------------
! Release internal MKL memory that might be used for computations
! NOTE: It is important to call the routine below to avoid memory leaks
! unless you disable MKL Memory Manager
!---------------------------------------------------------------------------
998   WRITE( *,'(A)') ' '
      WRITE( *,'(A,A)') '---------------------------------------------',
     1 '----------------------'
      WRITE( *,'(A,A)') 'Unfortunately, FGMRES+ILU0 Fortran example',
     1 ' has FAILED'
      WRITE( *,'(A,A)') '---------------------------------------------',
     1 '----------------------'
      WRITE( *,'(A)') ' '
      CALL MKL_FREEBUFFERS
      STOP 1

      END
