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
!===========================================================================
MODULE rndseed

    USE global

    implicit none

#ifdef HAVE_MPI
  include "mpif.h"
#endif

    private
    character(len=8) :: seedfname = 'rnd_seed'
    integer(i4b), parameter, public :: nRT = 10
    integer(i4b), parameter, public :: max_rnd_seed = 100000000

    public :: get_seed, check_rnd, rnd_dist

    interface rnd_dist
      module procedure rnd_d1, rnd_d2, rnd_d3
    end interface

CONTAINS
!---------------------------------------------------------------------------
SUBROUTINE get_seed(seed)
!  =========================================================================
! |  Read the random seed that is saved in the first line of the random     |
! |  seed file.  If the file does not exist, generate a seed from the date. |
!  =========================================================================
    
    implicit none

    integer(i4b), intent(OUT) :: seed
    integer, parameter :: IOU = 90
    integer :: istat

    if ( seedfile_exists() ) then
        open(UNIT=IOU, FILE=seedfname, STATUS='OLD', ACTION='READ', &
                FORM='FORMATTED', IOSTAT=istat)
        read(IOU, "(i11)", IOSTAT=istat) seed
        close(IOU)
    else
        seed = randate()
    end if

    RETURN

END SUBROUTINE get_seed
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
SUBROUTINE check_rnd(seed, test_distribution)
!  =========================================================================
! |  Compare a random distribution with one saved in a file to verify that  |
! |  the random seed is consistent across platforms.  If there is no file,  |
! |  write the input distribution and seed to a new file.                   |
!  =========================================================================
    
    implicit none

    integer(i4b), intent(IN) :: seed
    real(DP), dimension(nRT), intent(IN) :: test_distribution
    integer(i4b), dimension(nRT) :: int_dist, saved_dist, difference
    integer, parameter :: IOU = 95
    integer :: i, istat, skip

! Convert test_distribution to an array of nine-digit integers.
    int_dist = int(test_distribution*1000000000)

    if ( seedfile_exists() ) then    ! check it
    ! open the file
        open(UNIT=IOU, FILE=seedfname, STATUS='OLD', ACTION='READ', &
                FORM='FORMATTED', IOSTAT=istat)
    ! read lines 2 to 11
        read(IOU, "(i11)", IOSTAT=istat) skip
        do i=1,10
            read(IOU, "(i11)", IOSTAT=istat) saved_dist(i)
        end do
    ! close the file
        close(IOU)
    ! compare to the input test dist.
        difference = abs(int_dist - saved_dist)
    ! If they're not the same print a warning.
        if (sum(difference) .ne. 0) then
            print *, "Warning: Random number generator does not produce"
            print *, "numbers consistent with previous use of seed ", seed, "."
            print *, "The test distribution differs by ", difference
        end if
    else                        ! write one
    ! open the file
        open(UNIT=IOU, FILE=seedfname, STATUS='NEW', ACTION='WRITE', &
                FORM='FORMATTED', IOSTAT=istat)
    ! write the seed and the test dist.
        write(IOU,"(i11)") seed
        do i=1,10
            write(IOU,"(i11)") int_dist(i)
        end do
    ! close the file
        close(IOU)
    end if

    RETURN

END SUBROUTINE check_rnd
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
FUNCTION seedfile_exists() RESULT(yesorno)
!  =========================================================================
! |   Check to see if a random seed file exists.                            |
!  =========================================================================

    implicit none

    logical :: yesorno

!    if (iam == MASTER) then
        inquire(file=seedfname, exist=yesorno)
!    end if

!#ifdef HAVE_MPI
!    MPI_Bcast(yesorno, 1, MPI_LOGICAL, MASTER, MPI_COMM_WORLD, ierr)
!#endif

END FUNCTION seedfile_exists
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
FUNCTION randate() RESULT(rnd)
!===========================================================================
!   Generate a pseudo-random number out of the date_and_time() subroutine.
!   Brian Tuttle 25Feb05 <btuttle@psu.edu>
!===========================================================================

    implicit none

    integer :: rndint, maxrnd, rnd
    integer :: i, dig1, dig2, dig3, num1, num2, num3
    double precision :: nrnd
    character(len=10) :: timestr

    num1 = 0
    num2 = 0
    num3 = 0
! Get the time ("hhmmss.sss")
    call date_and_time(time=timestr)
! Convert the string into groups of three digits, avoiding the '.'
    do i=1,3
       dig1 =  ichar(timestr(i:i)) - ichar('0')
       dig2 =  ichar(timestr(i+3:i+3)) - ichar('0')
       dig3 =  ichar(timestr(i+7:i+7)) - ichar('0')
       num1 = 10*num1 + dig1
       num2 = 10*num2 + dig2
       num3 = 10*num3 + dig3
    end do
! Multiply them and add the last digit for good measure.
    rndint = num1*num2*num3 + dig3
! Normalize the result.
    maxrnd = 235 * 959 * 999 + 9
    nrnd = real(rndint)/real(maxrnd)
! Make the result a nine digit number max.
    rnd = int(max_rnd_seed * nrnd)

    RETURN

END FUNCTION randate
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
SUBROUTINE rnd_d1(randout)
!  =======================================================================
! |  Generates a one-dimensional array of random numbers                  |
!  =======================================================================

    USE mt19937

    implicit none

    double precision, dimension(:), intent(OUT) :: randout
    integer :: i

    do i = 1,size(randout,1)
        randout(i) = genrand_real1()
    end do

    RETURN

END SUBROUTINE rnd_d1
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
SUBROUTINE rnd_d2(randout)
!  =======================================================================
! |  Generates a two-dimensional array of random numbers                  |
!  =======================================================================

    USE mt19937

    implicit none

    double precision, dimension(:,:), intent(OUT) :: randout
    integer :: i, j

    do i = 1,size(randout,1)
      do j = 1,size(randout,2)
        randout(i,j) = genrand_real1()
      end do
    end do

    RETURN

END SUBROUTINE rnd_d2
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
SUBROUTINE rnd_d3(randout)
!  =======================================================================
! |  Generates a three-dimensional array of random numbers                |
!  =======================================================================

    USE mt19937

    implicit none

    double precision, dimension(:,:,:), intent(OUT) :: randout
    integer :: i, j, k

    do i = 1,size(randout,1)
      do j = 1,size(randout,2)
        do k = 1,size(randout,3)
          randout(i,j,k) = genrand_real1()
        end do
      end do
    end do

    RETURN

END SUBROUTINE rnd_d3
!---------------------------------------------------------------------------
END MODULE rndseed
