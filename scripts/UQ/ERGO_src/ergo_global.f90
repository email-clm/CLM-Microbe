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

MODULE global

    implicit none

    integer, parameter :: DP = Kind(1.0d0)
    integer, parameter :: SP = kind(1.0)
    integer, parameter :: i4b = Selected_Int_Kind(9)
    integer, parameter :: i8b = Selected_Int_Kind(18)
    integer, parameter :: MASTER = 0
    
! Model size variables:
    integer(i4b) :: n_param, n_param_pre ! number of parameters
    integer(i4b) :: n_pen, n_pen_pre     ! number of penalty coefficients
    integer(i4b) :: n_g
    integer(i4b) :: n_conv, n_conv_pre   ! number of parameters for convergence

    integer(i4b) :: ierr, n_procs, iam
    integer(i4b) :: model_iters
    integer(i4b) :: nfeval         ! number of objective function evaluations
    logical :: restarting

    integer(i4b) :: opt_strategy, opt_strategy_pre

! Random seed variables:
    integer(i4b) :: rand_seed
    logical :: seeded

CONTAINS

!---------------------------------------------------------------------------
FUNCTION stringcat(string1, string2) RESULT(catstrings)
!  ==========================================================================
! |  This function does a simple concatenation of two strings.  It is        |
! |  included here because the C preprocessor treats Fortran's concatenation |
! |  operator as a comment and truncates the line at that point.  This is    |
! |  a way to work around that behavior.                                     |
!  ==========================================================================

    character(len=*) :: string1, string2
    character(len=len(string1)+len(string2)) :: catstrings

    catstrings = string1 // string2

END FUNCTION stringcat
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
SUBROUTINE skipline(iounit, record, lines)
!  ========================================================================
! |  Reads a number of lines in a formatted, sequentially accessed file    |
! |  to position the file marker.                                          |
!  ========================================================================

    implicit none

    integer, intent(IN) :: iounit
    character(len=*), intent(OUT) :: record
    integer, intent(IN) :: lines
    integer :: i, istat

    do i = 1, lines
      read(UNIT=iounit, FMT = "(a)", IOSTAT=istat) record
    end do

    RETURN

END SUBROUTINE skipline
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
SUBROUTINE error_msg(errnum, msgstr)
!  ==========================================================================
! |  This subroutine prints a friendly error message and stops the program.  |
!  ==========================================================================

    implicit none

    integer, intent(IN) :: errnum
    character(len=*), intent(IN) :: msgstr

    select case (errnum)
        case (1)
            print *, "Error: Problem allocating memory for ",trim(msgstr),"."
        case (2)
            print *, "Error: Problem deallocating memory for ",trim(msgstr),"."
    end select

    STOP

END SUBROUTINE error_msg
!---------------------------------------------------------------------------

END MODULE global
