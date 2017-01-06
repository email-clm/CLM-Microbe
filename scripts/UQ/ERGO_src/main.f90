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
PROGRAM Main

    USE global
    USE ERGO
    USE fileIO
    USE model

    implicit none

    real(DP), dimension(:), allocatable :: Prerequisite, Solution

! Initialize the global environment.
    call ergo_init()

! Initialize main program variables.
    call alloc_main()
    Prerequisite = 0.0d0
    Solution = 0.0d0

! Check for a restart file.
    call check_restart(restarting)

! Optimize for prerequisite solution.
  if (opt_strategy_pre > 0 .and. .not. restarting) then
    call optimize_pre(Prerequisite)
!    print *, "Prerequisite: ", Prerequisite
  end if

! Optimize for the main solution.
  if (opt_strategy > 0) then
    call optimize(Prerequisite, Solution)
  end if
!  print *, "proc", iam, "Done optimize"

! Report the solution.
    if (opt_strategy_pre == 0 .and. opt_strategy == 0) then
        call test_model()
    else if (opt_strategy == 0) then
        call report(Prerequisite, Prerequisite)
    else
        call report(Prerequisite, Solution)
    end if
  if (iam == MASTER) then
    call write_endtime()
  end if

! Clean up.
    call dealloc_main()
    call finalize()

CONTAINS
!---------------------------------------------------------------------------
SUBROUTINE alloc_main()
!  ========================================================================
! |  This subroutine allocates memory for the main program.
!  ========================================================================

    integer :: astat

    astat = 0
    allocate(Prerequisite(n_param_pre), STAT=astat)
    if (astat > 0) call error_msg(1,"Prerequisite")
    astat = 0
    allocate(Solution(n_param), STAT=astat)
    if (astat > 0) call error_msg(1,"Solution")

    RETURN

END SUBROUTINE alloc_main
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
SUBROUTINE dealloc_main()
!  ========================================================================
! |  This subroutine deallocates memory for the main program.
!  ========================================================================

    integer :: astat

    astat = 0
    deallocate(Prerequisite, STAT=astat)
    if (astat > 0) call error_msg(2,"Prerequisite")
    astat = 0
    deallocate(Solution, STAT=astat)
    if (astat > 0) call error_msg(2,"Solution")

    RETURN

END SUBROUTINE dealloc_main
!---------------------------------------------------------------------------
END PROGRAM Main
