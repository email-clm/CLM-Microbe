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

MODULE timer

    USE global

    implicit none

#ifdef HAVE_MPI
   include "mpif.h"

#endif
    private
    integer(i4b), dimension(8), public :: time_0
    integer(i4b) :: this_year

    integer(i4b), dimension(8), public :: start_time
    integer(i4b), dimension(8), public :: end_time, last_restart
    real(DP), public :: restart_intvl

    public :: init_timer, time_check, time_diff, parse_time, update_timer
    public :: reset_timer

CONTAINS
!---------------------------------------------------------------------------
SUBROUTINE init_timer()
!  ========================================================================
! |  Initialize timer variables for stopping and restarting.               |
!  ========================================================================

    implicit none

! Get the time before starting initializations.
    call date_and_time(VALUES=time_0)
    this_year = time_0(1)

    RETURN

END SUBROUTINE init_timer
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
SUBROUTINE update_timer(quitntime)
!  ========================================================================
! |  Subtract elapsed time for the next stage of optimization.             |
!  ========================================================================

    implicit none

    real(DP) :: stage1time
    real(DP), intent(INOUT) :: quitntime

    stage1time = 0.0d0

 if (iam == MASTER) then

    call date_and_time(VALUES=start_time)

    last_restart = start_time

    call time_diff(time_0, start_time, hours = stage1time)

! Subtract initialization time from the allotted time
    quitntime = quitntime - stage1time

! Give the rest of the program at least three minutes to run.
! (Program may time out if the requested run time is too short in the
!  submit file.)
    if (quitntime < 0.05d0) quitntime = 0.05d0

 end if 

#ifdef HAVE_MPI
    call MPI_Bcast(quitntime, 1, MPI_Double_Precision, MASTER, MPI_Comm_World, ierr)
#endif
    RETURN

END SUBROUTINE update_timer
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
SUBROUTINE reset_timer()
!  ========================================================================
!  ========================================================================

    implicit none

 if (iam == MASTER) then
    call date_and_time(VALUES=start_time)
 end if 

#ifdef HAVE_MPI
    call MPI_Bcast(start_time, 8, MPI_Integer, MASTER, MPI_Comm_World, ierr)
#endif
    last_restart = start_time

    RETURN

END SUBROUTINE reset_timer
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
SUBROUTINE time_check(iterations, quitntime, max_total, max_iters, restart_flag, stop_flag)
!  ========================================================================
! |  This subroutine checks the wall clock time, decides whether to write  |
! |  a restart file and/or stop the current run, and estimates how many    |
! |  more iterations to do before checking the time again.                 |
!  ========================================================================

    implicit none

    integer(i4b), intent(IN) :: iterations, max_total
    real(DP), intent(IN) :: quitntime
    integer(i4b), intent(OUT) :: max_iters
    logical, intent(OUT) :: restart_flag, stop_flag
    integer(i4b), dimension(8) :: timechk
    real(DP) :: time_since_start, time_since_rstout, iters_per_hr
    integer(i4b) :: iters_to_end, iters_to_rst

! Check the time.
    call date_and_time(values = timechk)

! If restart interval has passed, then signal for a restart file.
    call time_diff(last_restart, timechk, hours = time_since_rstout)
    if (time_since_rstout > restart_intvl) then
        last_restart = timechk
        restart_flag = .true.
    end if

! If quitntime has passed, then signal for a restart file and an exit.
    call time_diff(start_time, timechk, hours = time_since_start)
    if ( (time_since_start >= quitntime .and. max_total == 0) .or. &
         (iterations >= max_total .and. max_total > 0) ) then
        restart_flag = .true.
        stop_flag = .true.
        max_iters = 10
        RETURN
    end if

! Measure the average time for an iteration.
    iters_per_hr = real(iterations) / time_since_start

! Estimate how many iterations it will take to go 60% of the remaining time.
    iters_to_end = int(0.6*(quitntime - time_since_start)*iters_per_hr)

! Estimate how many iterations it will take to go 60% of the restart interval.
    iters_to_rst = int(0.6*(restart_intvl - time_since_rstout)*iters_per_hr)

! Reset max_iters to the lesser but not less than 10.
    max_iters = max(10, min(iters_to_end, iters_to_rst)) 

    RETURN

END SUBROUTINE time_check
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
SUBROUTINE time_diff(time1, time2, components, hours)
!  ========================================================================
! |  This subroutine returns the difference between two time arrays as     |
! |  reported by the date_and_time(values) subroutine [year, month, day,   |
! |  gmt_min, hr, min, sec, msec].  The output is expressed as an array    |
! |  of integer components {days, hrs, minutes, seconds} and as a real     |
! |  number of hours.  Leap years are not accounted for.                   |
!  ========================================================================

    implicit none

    integer(i4b), dimension(8), intent(IN)  :: time1, time2
    integer(i4b), dimension(4), optional, intent(OUT) :: components
    real(DP), optional, intent(OUT) :: hours
    integer(i4b), dimension(12), parameter :: mlen = (/ 0, 31, 28, 31, 30, &
                                           31, 30, 31, 31, 30, 31, 30 /)
    integer(i4b) :: days1, days2, days, hrs, mins
    real(DP) :: secs1, secs2, secs

! Convert time1 to days since Jan 1 of this year and seconds
    days1 = (time1(1) - this_year) * 365 + sum(mlen(1:time1(2))) + time1(3)
    secs1 = time1(5)*3600.0d0 + time1(6)*60.0d0 + time1(7) + 0.001d0*time1(8)

! Convert time2 to days since Jan 1 of this year and seconds
    days2 = (time2(1) - this_year) * 365 + sum(mlen(1:time2(2))) + time2(3)
    secs2 = time2(5)*3600.0d0 + time2(6)*60.0d0 + time2(7) + 0.001d0*time2(8)

! Find the elapsed days and seconds.
    days = days2 - days1
    if (secs2 < secs1) then
        secs2 = secs2 + 24.0d0*3600.0d0
        days = days - 1
    end if
    secs = secs2 - secs1
    hrs = int(secs/3600.0d0)
    secs = secs - hrs*3600.0d0
    mins = int(secs/60.0d0)
    secs = secs - mins*60.0d0

    if (present (components)) then
        components = (/ days, hrs, mins, int(secs) /)
    end if

    if (present (hours)) then
        hours = days*24 + hrs + real(mins)/60.0d0 + real(secs)/3600.0d0
    end if

    RETURN

END SUBROUTINE time_diff
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
SUBROUTINE parse_time(in_str, hours)
!  ========================================================================
! |  This subroutine assumes that the time string is in hhh:mm format, but |
! |  is flexible regarding the number of digits.  It returns hours as a    |
! |  real number.                                                          |
!  ========================================================================

    implicit none

    character(len=*), intent(IN) :: in_str
    real(DP), intent(OUT) :: hours
    character(len=len(in_str)) :: time_str
    integer, dimension(2) :: numbers
    integer :: colon_position

    time_str = in_str

! Change the colons into commas.  There should be only one.
    do
        colon_position = index(time_str, ":")
        if (colon_position == 0) exit
        time_str(colon_position:colon_position) = ","
    end do

! Convert the string into an array of integers via an internal file.
    read (time_str, fmt=*) numbers

    hours = numbers(1) + real(numbers(2))/60
    !hours = 120.

    RETURN

END SUBROUTINE parse_time
!---------------------------------------------------------------------------

END MODULE timer
