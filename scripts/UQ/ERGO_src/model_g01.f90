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
    real(DP) :: DE_pen_val

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
    n_param = 13
    n_param_pre = 1
    n_pen = 9
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

! Default value for pen_val in DE.
    DE_pen_val = 10000d0

! Initialize evaluation index with default values.
    eval_index = (/ (i, i=1, size(eval_index)) /)

    if (present(prereq_in)) then
        nop = n_param
        n_g = n_pen
    ! Allocate bounds.
        call model_allocate(nop)
    ! Initialize the model for the main optimization.
        x_lower = (/ 0,0,0,0,0,0,0,0,0,0,0,0,0 /)
        x_upper = (/ 1,1,1,1,1,1,1,1,1,100,100,100,1 /)

    else
        nop = n_param_pre
        n_g = n_pen_pre
    ! Allocate bounds.
        call model_allocate(nop)
    ! Initialize the model for the pre-calculation.

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
    real(DP), dimension(n_g), optional, intent(OUT) :: g
    integer :: nfe, i

! Initialize output.
    f = 0d0
    g = 0d0

! Default number of function evaluations is 1 for each call to model_objf.
    nfe = 1

! objective function
    f = 5.0*(x(1)+x(2)+x(3)+x(4)) - 5.0*(x(1)**2+x(2)**2+x(3)**2+x(4)**2)
    Do i = 5, 13
      f = f - x(i)
    End Do
! constraints g<=0
    g(1) =  2.0*x(1)+2.0*x(2)+x(10)+x(11)-10.
    g(2) =  2.0*x(1)+2.0*x(3)+x(10)+x(12)-10.
    g(3) =  2.0*x(2)+2.0*x(3)+x(11)+x(12)-10.
    g(4) = -8.0*x(1)+x(10)
    g(5) = -8.0*x(2)+x(11)
    g(6) = -8.0*x(3)+x(12)
    g(7) = -2.0*x(4)-x(5)+x(10)
    g(8) = -2.0*x(6)-x(7)+x(11)
    g(9) = -2.0*x(8)-x(9)+x(12)

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
SUBROUTINE setup_test_objf(test_prereq_soln, test_solution)
!  ========================================================================
! |  Provide a known set of parameters to test the model objective         |
! |  function.                                                             |
!  ========================================================================

    implicit none

    real(DP), dimension(n_param_pre) :: test_prereq_soln
    real(DP), dimension(n_param) :: test_solution


    RETURN

END SUBROUTINE setup_test_objf
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
SUBROUTINE model_out(Solution, obf_out)
!  ========================================================================
! |  This output routine is called by the optimization program as it       |
! |  progresses toward convergence on a solution.  It allows the model to  |
! |  report the state of its various parameters given the latest best      |
! |  solution.  It is also called at the end of the optimization process.  |
!  ========================================================================

    implicit none

    real(DP), dimension(:), intent(IN) :: Solution
    real(DP), intent(OUT), OPTIONAL :: obf_out
    character(len=12) :: fname = "solution_out"
    real(DP), dimension(n_g) :: g
    real(DP) :: fout
    integer :: FID = 10
    integer :: istat, i

    call model_objf(Solution,fout,g)

    open(UNIT=FID, FILE=fname, STATUS='REPLACE', ACTION='WRITE', &
        ACCESS='SEQUENTIAL', IOSTAT=istat)

    write(FID,"(a,f6.2)") '# g01:  f(x) = ', fout
    write(FID,"(a)") '# i    x(i)'

    do i = 1,size(Solution)
        write(FID, "(i3,f6.1)") i, Solution(i)
    end do

    close(FID)

    if (present(obf_out)) obf_out = fout

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
    if (astat > 0) call error_msg(2,"x_upper")
    astat = 0
    deallocate(x_lower, STAT=astat)
    if (astat > 0) call error_msg(2,"x_lower")

    RETURN

END SUBROUTINE model_finalize
!---------------------------------------------------------------------------
END MODULE model
