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
!==========================================================================
!  This module provides an interface between the ERGO optimization package
!  and the numerical model that serves as an objective function.  Its
!  purpose is to centralize calls to the necessary model subroutines and
!  facilitate the incorporation of different models with the ERGO package.
!  The model itself may be defined here, or may be included as an established
!  module via a USE statement at the top with calls to its particular 
!  subroutines in the appropriate places.
!==========================================================================
MODULE model

    USE global

    implicit none

    integer(i4b) :: nop
    real(DP), dimension(:), allocatable :: x_lower, x_upper

CONTAINS
!---------------------------------------------------------------------------
SUBROUTINE read_model_input(fid)
!  =========================================================================
! |  This subroutine is called from read_input() in module fileIO which     |
! |  opens the input file, reads the optimization section, and positions    |
! |  the pointer at the beginning of the model_specific section.  This      | 
! |  routine continues reading the formatted data and assigns values to     |
! |  model variables.                                                       |
!  =========================================================================

    implicit none

! ID of a currently open file.
    integer(i4b), intent(IN) :: fid


    RETURN

END SUBROUTINE read_model_input
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
SUBROUTINE print_model_input()
!  ========================================================================
! |  This routine prints out the variables read from the input file and    |
! |  those computed in the model_init routine.  It is called from          |
! |  print_input() in ergo_fileIO.f90 and performs the model-specific      |
! |  part of that function.                                                |
!  ========================================================================
    
    implicit none


    RETURN

END SUBROUTINE print_model_input
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
SUBROUTINE model_write_arch(FID)
!  ========================================================================
! |  This routine writes the input variables to a file archive that can be |
! |  used as an input file to reproduce the calculation.  This is called   |
! |  from write_in_arch() in ergo_fileIO.f90.                              |
!  ========================================================================
    
    implicit none

! ID of a currently open file.
    integer(i4b), intent(IN) :: FID


    RETURN

END SUBROUTINE model_write_arch
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
SUBROUTINE model_size()
!  =========================================================================
! |  This subroutine conveys to the optimization package the sizes of the   |
! |  one-dimensional arrays for model parameters, penalty, and convergence. |
!  =========================================================================

    implicit none

! Default values for model sizes.
    n_param = 1
    n_param_pre = 1
    n_pen = 1
    n_pen_pre = 1
    n_conv = n_param
    n_conv_pre = n_param_pre


    RETURN

END SUBROUTINE model_size
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
SUBROUTINE model_init(eval_index, prereq_in)
!  ========================================================================
! |  This routine does any initialization that the model may require.  Any |
! |  prerequisite solution array may be included in the second argument.   |
! |  The output array is a list of parameter indices to be included by the |
! |  optimization program in evaluating the convergence of the solution.   |
! |  The default is to use all parameters for convergence evaluation.      |
!  ========================================================================

    implicit none

    integer(i4b), dimension(:), intent(OUT) :: eval_index
    real(DP), dimension(:), intent(IN), optional :: prereq_in
    integer :: i, astat

! Initialize evaluation index with default values.
    eval_index = (/ (i, i=1, size(eval_index)) /)

! Initialize the model for the PRE-CALCULATION.
    if (.not. present(prereq_in)) then
        nop = n_param_pre
        n_g = n_pen_pre
    ! Allocate bounds.
        call model_allocate(nop)
    ! Set lower and upper limits.
        x_lower = 0.0d0
        x_upper = 1.0d0

! Initialize the model for the MAIN OPTIMIZATION.
    else
        nop = n_param
        n_g = n_pen
    ! Allocate bounds.
        call model_allocate(nop)
    ! Set lower and upper limits.
        x_lower = 0.0d0
        x_upper = 1.0d0
    end if

    RETURN

END SUBROUTINE model_init
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
SUBROUTINE model_objf(x,f,g)
!  ========================================================================
! |  This is the model objective function.  The input is an array of       |
! |  decision parameters x.  The output is the single fitness value that   | 
! |  rates the choice of parameters.  The constraint vector g provides a   |
! |  mechanism to penalize infeasible parameter choices.                   |
!  ========================================================================

    implicit none

    real(DP), dimension(nop), intent(IN) :: x
    real(DP), intent(OUT) :: f
    real(DP), dimension(n_g) intent(OUT) :: g
    integer :: nfe

! Initialize output.
    f = 0d0
    g = 0d0

! Default number of function evaluations is 1 for each call to model_objf.
    nfe = 1

    nfeval = nfeval + nfe

    RETURN

END SUBROUTINE model_objf
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
SUBROUTINE model_init_outfiles()
!  ========================================================================
! |  This routine provides an opportunity to initialize output files that  |
! |  will be updated by the model throughout the optimization process.     |
! |  File IDs for open files may be established here as well as header     |
! |  information for the output files.  It is called from init_outfiles()  |
! |  in ergo_fileIO.f90.                                                   |
!  ========================================================================

    implicit none

    RETURN

END SUBROUTINE model_init_outfiles
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
SUBROUTINE model_out(Solution)
!  ========================================================================
! |  This output routine is called by the optimization program as it       |
! |  progresses toward convergence on a solution.  It allows the model to  |
! |  report the state of its various parameters given the latest best      |
! |  solution.  It is also called at the end of the optimization process.  |
!  ========================================================================

    implicit none

    real(DP), dimension(:), intent(IN) :: Solution

    RETURN

END SUBROUTINE model_out
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
SUBROUTINE model_allocate(n_bounds)
!  ========================================================================
! |  Allocate upper and lower constraint arrays.                           |
!  ========================================================================

    implicit none

    integer, intent(IN) :: n_bounds
    integer :: astat

    astat = 0
    allocate(x_upper(n_bounds), STAT=astat)
    if (astat > 0) call error_msg(1,"x_upper")
    astat = 0
    allocate(x_lower(n_bounds), STAT=astat)
    if (astat > 0) call error_msg(1,"x_lower")

    RETURN

END SUBROUTINE model_allocate
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
SUBROUTINE model_finalize()
!  ========================================================================
! |  This routine cleans up memory allocated for model variables.          |
!  ========================================================================

    implicit none

    integer :: astat


! Deallocate upper and lower constraint arrays.
    astat = 0
    deallocate(x_upper, STAT=astat)
    if (astat > 0) call error_msg(1,"x_upper")
    astat = 0
    deallocate(x_lower, STAT=astat)
    if (astat > 0) call error_msg(1,"x_lower")

    RETURN

END SUBROUTINE model_finalize
!---------------------------------------------------------------------------
END MODULE model
