!    ERGO: Efficient Robust Global Optimization 
!    Copyright (C) 2007  Klaus Keller and Brian Tuttle 
!
!    This program is free software; you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation; either version 2 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program; if not, write to the Free Software
!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!    Klaus Keller, kkeller@geosc.psu.edu
!    Brian Tuttle, btuttle@psu.edu
!    David McInerny, dmcinern@geosc.psu.edu
!    Alex Robinson, ajr225@psu.edu
!

!---------------------------------------------------------------------------
MODULE ERGO

    USE global

    implicit none

#ifdef HAVE_MPI
   include "mpif.h"

#endif
    real(DP) :: pre_time, ST1_time, ST2_time, ST3_time

CONTAINS
!---------------------------------------------------------------------------
SUBROUTINE ergo_init()
!  ========================================================================
! |  This routine initializes ERGO by reading the input file, calculating  |
! |  sizes of global variables, and initializing them.                     |
!  ========================================================================

    USE fileIO
    USE input_vars
    USE timer
    USE rndseed
    USE model

    implicit none

! Orient for multiple processors.
#ifdef HAVE_MPI
! Initialize MPI
    call MPI_Init(ierr)
    print*, ierr
    call MPI_Comm_rank(MPI_COMM_WORLD, iam, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, n_procs, ierr)
#else
    iam = MASTER
    n_procs = 1
#endif

    print*, n_procs
! Start the timer.
  if (iam == MASTER) then
    call init_timer()
  end if

! Read the input file.
    call read_input()

! Initialize random number seed.
    seeded = .false.
  if (iam == MASTER) then
    call get_seed(rand_seed)
  end if

! Convert time strings to hours.
  if (iam == MASTER) then
    call parse_time(rtimePRE_str, pre_time)
    call parse_time(rtimeST1_str, ST1_time)
    call parse_time(rtimeST2_str, ST2_time)
    call parse_time(rtimeST3_str, ST3_time)
  end if

! Find out sizes of variables.
    call model_size()
!  Find number of parameters (nop_1) for BAU without learning.
!    nop_1 = int(real(mstop - mstart) / real(mstep)) + 1

!  Find nop with learning (nop_2) and multiple states of world.
!    call init_learn()

! Allocate model variables.
!    call alloc_ergo()

! Initialize model variables.

! Report the state of the model.
    if (iam == MASTER) then
        call print_input()
        call init_outfiles()
    end if

    RETURN

END SUBROUTINE ergo_init
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
SUBROUTINE optimize_pre(Prereq_out)
!  =========================================================================
! |  This routine prepares the Dice model for a Business-As-Usual scenario  |
! |  (no damages due to MOC or WAIS collapse and abatement policy of zero)  |
! |  and optimizes the investment policy for use later.                     |
!  =========================================================================

    USE timer
    USE input_vars 
    USE converge
    USE ceo
    USE de
    USE model

    implicit none

    real(DP), dimension(n_param_pre), intent(OUT) :: Prereq_out
    real(DP), dimension(n_param_pre,2) :: coarse_soln
    real(DP), dimension(n_param_pre,2) :: sres_sigma1, sres_sigma2
    real(DP), dimension(n_param_pre,1) :: Solution_out
!    integer(i4b), dimension(n_conv_pre) :: conv_indx
    integer(i4b) :: maxiter_null = 0    ! dummy placeholder maximum iterations.

!    print *, "Starting optimize_pre."

#ifdef ERGO
! Set DE parameters to their original values since this part runs fine, and
! we are interested in comparing results from the next optimization.
    inDE_pop = 10
    inDE_rndscale = 0
    inDE_strat = 6
    inDE_F_CR = 0.8d0
#endif
!    writing_Cout = .false.
    nfeval = 0

! Setup model for prerequisite calculation.
    nop = n_param_pre
    n_g = n_pen_pre
    cx_len = n_conv_pre
    call alloc_conv_indx()
    call model_init(conv_indx)
!    print *, "PRE- conv_indx:", conv_indx

! Reset total number of model iterations.
    model_iters = 0

! Choose optimization strategy and run it.
    select case (opt_strategy_pre)
      case (1)      ! DE(2)
        call run_DE(Solution_out, maxiter_PRE, pre_time, RMStgt_pre, 2)

      case (2)      ! SRES(2)
        call run_sres(Solution_out, sres_sigma1, maxiter_PRE, pre_time, &
                        RMStgt_pre, 2)

      case (3)      ! DE(2) => SRES(2)
        print *, "Please note that this optimization stategy is not necessarily"
        print *, "repeatable for calculation of prerequisite solution."
        call run_DE(coarse_soln, maxiter_null, pre_time, RMStgt_de, 2)

        call update_timer(pre_time)
        call run_sres(Solution_out, sres_sigma1, maxiter_null, pre_time, &
                        RMStgt_sres, 3, INIT_POP = coarse_soln)

      case (4)      ! SRES(2) => ISRES(2)
        print *, "Please note that this optimization stategy is not necessarily"
        print *, "repeatable for calculation of prerequisite solution."
        call run_sres(coarse_soln, sres_sigma1, maxiter_null, pre_time, &
                        RMStgt_sres, 2)

        call update_timer(pre_time)
        call run_sres(Solution_out, sres_sigma2, maxiter_null, pre_time, &
                        RMStgt_pre, 2, INIT_POP = coarse_soln, &
                        INIT_SIGMA = sres_sigma1, &
                        IN_ALPHA = isres_alpha)

      case (5)
        print *, "Please note that this optimization stategy is not necessarily"
        print *, "repeatable for calculation of prerequisite solution."
        call run_sres(coarse_soln, sres_sigma1, maxiter_null, pre_time, &
                        RMStgt_sres, 3)

        call update_timer(pre_time)
        call run_sres(Solution_out, sres_sigma2, maxiter_null, pre_time, &
                        RMStgt_pre, 2, INIT_POP = coarse_soln, &
                        INIT_SIGMA = sres_sigma1, &
                        IN_ALPHA = isres_alpha)

      case (6)
        print *, "Please note that this optimization stategy is not necessarily"
        print *, "repeatable for calculation of prerequisite solution."
        call run_DE(coarse_soln, maxiter_null, pre_time, RMStgt_de, 3)

        call update_timer(pre_time)
        call run_sres(coarse_soln, sres_sigma1, maxiter_null, pre_time, &
                        RMStgt_sres, 3, INIT_POP = coarse_soln)

        call update_timer(pre_time)
        call run_sres(Solution_out, sres_sigma2, maxiter_null, pre_time, &
                        RMStgt_pre, 2, INIT_POP = coarse_soln, &
                        INIT_SIGMA = sres_sigma1, &
                        IN_ALPHA = isres_alpha)

    end select

! Clean up Dice and Convergence variables.
    call dealloc_conv_indx()
    call model_finalize()

! Output the solution.
    Prereq_out = Solution_out(:,1)

!    print *, "proc", iam, "Finished optimize_pre."

    RETURN

END SUBROUTINE optimize_pre
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
SUBROUTINE optimize(Prereq_in, Soln_out)
!  =========================================================================
! |  This subroutine sets up the Dice model with a prescribed investment    |
! |  policy and optimizes with respect to the abatement policy.             |
!  =========================================================================

    USE timer
    USE input_vars 
    USE converge
    USE ceo
    USE de
    USE model

    implicit none

    real(DP), dimension(n_param_pre), intent(IN) :: Prereq_in
    real(DP), dimension(n_param), intent(OUT) :: Soln_out
    real(DP), dimension(n_param,2) :: coarse_soln
    real(DP), dimension(n_param,2) :: sres_sigma1, sres_sigma2 
    real(DP), dimension(n_param,1) :: Solution_out

#ifdef ERGO
    inDE_pop = DE_pop
    inDE_rndscale = DE_rndscale 
    inDE_strat = DE_strat
    inDE_F_CR = DE_F_CR
#endif

! Reset the timer.
    call reset_timer()

!    writing_Cout = .true.
    nfeval = 0

! Setup model for optimization.
    nop = n_param
    n_g = n_pen 
    cx_len = n_conv
    call alloc_conv_indx()
    call model_init(conv_indx, Prereq_in)
!    print *, "conv_indx:", conv_indx

! Set up parameters for endogenous detection of convergence.
!    call init_conv(YOLx)

! Reset total number of model iterations.
    model_iters = 0

! Do the optimization.
    !print*, 'testo1'
    select case (opt_strategy)
      case (1)      ! DE(2)
        call run_DE(Solution_out, maxiter_ST1, ST1_time, RMStgt_de, 2)

      case (2)      ! SRES(2)
        call run_sres(Solution_out, sres_sigma1, maxiter_ST1, ST1_time, &
                        RMStgt_sres, 2)

      case (3)      ! DE(3)
        call run_DE(Solution_out, maxiter_ST1, ST1_time, RMStgt_de, 3)

      case (4)      ! SRES(3)
        call run_sres(Solution_out, sres_sigma1, maxiter_ST1, ST1_time, &
                        RMStgt_sres, 3)

      case (5)      ! DE(2) => SRES(2)
        call run_DE(coarse_soln, maxiter_ST1, ST1_time, RMStgt_de, 2)
        call reset_timer()
        call run_sres(Solution_out, sres_sigma1, maxiter_ST2, ST2_time, &
                        RMStgt_sres, 2, INIT_POP = coarse_soln)

      case (6)      ! SRES(2) => ISRES(2)
        call run_sres(coarse_soln, sres_sigma1, maxiter_ST1, ST1_time, &
                        RMStgt_sres, 2)
        call reset_timer()
        call run_sres(Solution_out, sres_sigma2, maxiter_ST2, ST2_time, &
                        RMStgt_isres, 2, INIT_POP = coarse_soln, &
                        INIT_SIGMA = sres_sigma1, &
                        IN_ALPHA = isres_alpha)

      case (7)      ! SRES(3) => ISRES(2)
        call run_sres(coarse_soln, sres_sigma1, maxiter_ST1, ST1_time, &
                        RMStgt_sres, 3)
        
        call reset_timer()
        call run_sres(Solution_out, sres_sigma2, maxiter_ST2, ST2_time, &
                        RMStgt_isres, 2, INIT_POP = coarse_soln, &
                        INIT_SIGMA = sres_sigma1, &
                        IN_ALPHA = isres_alpha)

      case (8)      ! DE(3) => SRES(3) => ISRES(2)
        call run_DE(coarse_soln, maxiter_ST1, ST1_time, RMStgt_de, 3)

        call reset_timer()
        call run_sres(coarse_soln, sres_sigma1, maxiter_ST2, ST2_time, &
                        RMStgt_sres, 3, INIT_POP = coarse_soln)

        call reset_timer()
        call run_sres(Solution_out, sres_sigma2, maxiter_ST3, ST3_time, &
                        RMStgt_isres, 2, INIT_POP = coarse_soln, &
                        INIT_SIGMA = sres_sigma1, &
                        IN_ALPHA = isres_alpha)

      case (9)      ! ISRES(3)
        call run_sres(Solution_out, sres_sigma1, maxiter_ST1, ST1_time, &
                        RMStgt_isres, 3, IN_ALPHA = isres_alpha)

    end select
    prinT*, 'testo2'

! Clean up Dice and Convergence variables.
    call dealloc_conv_indx()
    call model_finalize()

! Output the solution.
    Soln_out = Solution_out(:,1)

    RETURN

END SUBROUTINE optimize
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
SUBROUTINE test_model()
!  ========================================================================
! |  Run the objective function with a known solution and print out the    |
! |  objective function value.  This is called from main.f90 in place of   | 
! |  report(Prerequisite, Solution) if all optimization is skipped.        |
!  ========================================================================

    USE converge
    USE model

    implicit none

    real(DP), dimension(n_param_pre) :: Prereq_in
    real(DP), dimension(n_param) :: X_in
    real(DP), dimension(n_pen) :: g
    real(DP) :: f

    call setup_test_objf(Prereq_in, X_in)

    call report(Prereq_in, X_in, f)

    if (iam == MASTER) then
        print *, "f = ", f
    end if

!! Setup model for optimization.
!    nop = n_param
!    n_g = n_pen 
!    cx_len = n_conv
!    print *, "Calling alloc_conv_indx..."
!    call alloc_conv_indx()
!    print *, "Done alloc_conv_indx."
!    print *, "Calling model_init..."
!    call model_init(conv_indx, Prereq_in)
!    print *, "Done model_init."
!
!
!    print *, "Calling model_objf..."
!    call model_objf(X_in, f, g)
!    print *, "Done model_objf."
!
!    print *, "f = ", f
!
!    call dealloc_conv_indx()
!    call model_finalize()

    RETURN

END SUBROUTINE test_model
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
SUBROUTINE report(Prereq, Soln, Objf)
!  =========================================================================
! |  Write out the final solution in model variables and the create the     |
! |  input file that will reproduce the model run.                          |
!  =========================================================================

    USE fileIO
    USE input_vars
    USE converge
    USE model

    implicit none

    real(DP), dimension(n_param_pre), intent(IN) :: Prereq
    real(DP), dimension(n_param), intent(IN) :: Soln
    real(DP), intent(OUT), OPTIONAL :: Objf
    
!    writing_Cout = .false.

! Set up the model.
    nop = n_param
    n_g = n_pen 
    cx_len = n_conv
    call alloc_conv_indx()

    call model_init(conv_indx, Prereq)


! Write out the solution.
    if (present(Objf)) then
      Objf = 0.0d0          ! Initialize for non-Masters
      if (iam == MASTER) then
        call model_out(Soln, out_path, Objf)
      end if
    else
      if (iam == MASTER) then
        call model_out(Soln, out_path)
      end if
    end if

! Clean up model variables.
    call model_finalize()

! Write archive of input parameters.
    if (iam == MASTER) then
        call write_in_arch()
    end if

    RETURN

END SUBROUTINE report
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
SUBROUTINE finalize()
!  ========================================================================
! |  Clean up ERGO variables and finalize.                                 |
!  ========================================================================

    implicit none

! Deallocate model variables.
    call dealloc_ergo()

#ifdef HAVE_MPI
    Call Mpi_Finalize(ierr)
#endif

    RETURN

END SUBROUTINE finalize
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
SUBROUTINE alloc_ergo()
!  ========================================================================
! |  Allocate ERGO and global module variables.                           |
!  ========================================================================

    implicit none

    RETURN

END SUBROUTINE alloc_ergo
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
SUBROUTINE dealloc_ergo()
!  ========================================================================
! |  Deallocate ERGO and global module variables.                         |
!  ========================================================================

    implicit none

    RETURN

END SUBROUTINE dealloc_ergo
!---------------------------------------------------------------------------
END MODULE ERGO
