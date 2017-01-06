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

MODULE fileIO

    USE global

    implicit none

#ifdef HAVE_MPI
   include "mpif.h"

#endif

    !private

!    character(len=5) :: Coutfname = 'C_out'
    character(len=7) :: restfname = 'restart'
    character(len=13), parameter :: infname = './input'
    character(len=13), parameter :: inarchfname = './input.archive'
!    character(len=8) :: seedfname = 'rnd_seed'
    character(len=15) :: stpolfname = 'starting_policy'
    character(len=3), dimension(12), parameter :: month = (/ 'Jan','Feb', &
              'Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec' /)
    character(len=80) :: out_fbest, out_xbest, out_conv, out_restart !, out_dvar
    character(len=40) :: fbestname, xbestname, convfname !, dvarname
    logical :: writing_conv
    integer, parameter :: ioCout = 72

    
    character(len=50) :: out_path

    public :: init_outfiles, write_fbest, write_endtime, write_conv
    public :: read_input, print_input, write_xout, write_in_arch 
    public :: check_restart 

CONTAINS

!---------------------------------------------------------------------------
SUBROUTINE read_input()
!  ========================================================================
! |  This subroutine reads from a file called input, assigns values to     |
! |  the input structure in mod_global, and broadcasts the input data to   |
! |  the rest of the processors if any.                                    |
!  ========================================================================

    USE input_vars
    USE model

    implicit none

    integer, parameter :: tabstop = 31
    integer, parameter :: IOU = 10
    character(len=tabstop) :: tabover
    character(len=75) :: sec_head
    integer :: istat, i
    logical :: archivefile

if (iam == MASTER) then
! Check for input archive file.
    archivefile = .false.
    inquire(file=inarchfname, exist=archivefile)
    if (archivefile) then
    ! Open input archive file.
        print *, "Using ", inarchfname
        open(UNIT=IOU, FILE=inarchfname, STATUS='OLD', ACTION='READ', &
             FORM='FORMATTED', IOSTAT=istat)
    else
    ! Open input file.
        open(UNIT=IOU, FILE=infname, STATUS='OLD', ACTION='READ', &
             FORM='FORMATTED', IOSTAT=istat)
    end if

  headers: do
    read(IOU, "(a)", IOSTAT=istat) sec_head

    select case (trim(sec_head(1:35)))
    case ('Run time (wallclock time) [hhh:mm]:')
      rtime: do
        read(IOU, "(a)", ADVANCE="NO", IOSTAT=istat) tabover
        select case (tabover)
        case ('    Precalculation           = ')
           read(IOU, "(a)", IOSTAT=istat) rtimePRE_str
        case ('    Stage 1                  = ')
           read(IOU, "(a)", IOSTAT=istat) rtimeST1_str
        case ('    Stage 2                  = ')
           read(IOU, "(a)", IOSTAT=istat) rtimeST2_str
        case ('    Stage 3                  = ')
           read(IOU, "(a)", IOSTAT=istat) rtimeST3_str
        case ('    Maximum iterations pre   = ')
           read(IOU, "(i12)", IOSTAT=istat) maxiter_PRE
        case ('    Maximum iterations st1   = ')
           read(IOU, "(i12)", IOSTAT=istat) maxiter_ST1
        case ('    Maximum iterations st2   = ')
           read(IOU, "(i12)", IOSTAT=istat) maxiter_ST2
        case ('    Maximum iterations st3   = ')
           read(IOU, "(i12)", IOSTAT=istat) maxiter_ST3
        case ('                               ')
           exit rtime
        end select
      end do rtime
    case ('Optimization strategies:')
      optstrat: do
        read(IOU, "(a)", ADVANCE="NO", IOSTAT=istat) tabover
        select case (tabover(1:18))
        case ('    Precalculation')
           read(IOU, "(i1)", IOSTAT=istat) opt_strategy_pre
        case ('    Model optimiza')
           read(IOU, "(i1)", IOSTAT=istat) opt_strategy
        case ('                  ')
           exit optstrat
        end select
      end do optstrat
    case ('Convergence criteria:')
      ccrit: do
        read(IOU, "(a)", ADVANCE="NO", IOSTAT=istat) tabover
        select case (tabover)
        case ('    Target RMSE (BAU)   [%]  = ')
           read(IOU, "(f17.15)", IOSTAT=istat) RMStgt_pre
        case ('    Target RMSE (DE)    [%]  = ')
           read(IOU, "(f17.15)", IOSTAT=istat) RMStgt_de
        case ('    Target RMSE (SRES)  [%]  = ')
           read(IOU, "(f17.15)", IOSTAT=istat) RMStgt_sres
        case ('    Target RMSE (ISRES) [%]  = ')
           read(IOU, "(f17.15)", IOSTAT=istat) RMStgt_isres
        case ('                               ')
           exit ccrit
        end select
      end do ccrit
    case ('Thresholds for 3-chain strategy:')
      chain3: do
        read(IOU, "(a)", ADVANCE="NO", IOSTAT=istat) tabover
        select case (tabover)
        case ('    RMSE between Q2 & Q3 [%] = ')
           read(IOU, "(f17.15)", IOSTAT=istat) RMSrej_Q23
        case ('    Weighted Av. delta RMSE  = ')
           read(IOU, "(f17.15)", IOSTAT=istat) THdRMSE
        case ('    Weighted Av. delta Objf  = ')
           read(IOU, "(f17.15)", IOSTAT=istat) THdObjF
        case ('                               ')
           exit chain3
        end select
      end do chain3
    case ('ISRES parameters:')
      is_param: do
        read(IOU, "(a)", ADVANCE="NO", IOSTAT=istat) tabover
        select case (tabover)
        case ('    alpha [default 0.2]      = ')
           read(IOU, "(f7.4)", IOSTAT=istat) isres_alpha
        case ('                               ')
           exit is_param
        end select
      end do is_param
    case ('Output files (max 40 characters):')
      outfiles: do
        read(IOU, "(a)", ADVANCE="NO", IOSTAT=istat) tabover
        select case (tabover)
        case ('    Path                     = ')
           read(IOU, "(a)", IOSTAT=istat) out_path
        case ('    fbest                    = ')
           read(IOU, "(a)", IOSTAT=istat) fbestname
        case ('    xbest (set of files)     = ')
           read(IOU, "(a)", IOSTAT=istat) xbestname
        case ('    convergence (optional)   = ')
           read(IOU, "(a)", IOSTAT=istat) convfname
        case ('                               ')
           exit outfiles
        end select
      end do outfiles
    case ('------------------------ MODEL PARA')
      exit headers
    end select
  end do headers

!! Space through the header.
!    call skipline(IOU, tabover, 9)
!! --- Run time (wallclock) ---
!! --- Run time (iterations) ---
!    call skipline(IOU, tabover, 2)
!! --- Optimization strategies ---
!    call skipline(IOU, tabover, 2)
!    call skipline(IOU, tabover, 4)
!! --- Convergence criteria ---
!    call skipline(IOU, tabover, 2)
!! --- 3-chain strategy Thresholds ---
!    call skipline(IOU, tabover, 2)
!! --- ISRES parameters ---
!    call skipline(IOU, tabover, 2)
!#ifdef ERGO
!! --- DE parameters ---
!    call skipline(IOU, tabover, 2)
!    read(IOU, "(a)", ADVANCE='NO', IOSTAT=istat) tabover
!    read(IOU, "(i3)", IOSTAT=istat) DE_pop
!    read(IOU, "(a)", ADVANCE='NO', IOSTAT=istat) tabover
!    read(IOU, "(i1)", IOSTAT=istat) DE_rndscale
!    read(IOU, "(a)", ADVANCE='NO', IOSTAT=istat) tabover
!    read(IOU, "(i1)", IOSTAT=istat) DE_strat
!    read(IOU, "(a)", ADVANCE='NO', IOSTAT=istat) tabover
!    read(IOU, "(f7.4)", IOSTAT=istat) DE_F_CR
!#endif
!! --- Output files ---
!    call skipline(IOU, tabover, 2)
!    call skipline(IOU, tabover, 2)

end if

! --- Model parameters
    call read_model_input(IOU)

if (iam == MASTER) then

! Close input file.
    close(IOU)

! Construct file names.
    out_fbest = stringcat(trim(out_path), stringcat("/", trim(fbestname)))
    out_xbest = stringcat(trim(out_path), stringcat("/", trim(xbestname)))
    out_restart = stringcat(trim(out_path), stringcat("/", trim(restfname)))

    writing_conv = .false.
    if (convfname .ne. '') then
        out_conv = stringcat(trim(out_path), stringcat("/", trim(convfname)))
        writing_conv = .true.
    end if

end if

#ifdef HAVE_MPI
    call MPI_Bcast(maxiter_PRE, 1, MPI_Integer, MASTER, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(maxiter_ST1, 1, MPI_Integer, MASTER, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(maxiter_ST2, 1, MPI_Integer, MASTER, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(maxiter_ST3, 1, MPI_Integer, MASTER, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(opt_strategy_pre, 1, MPI_Integer, MASTER, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(opt_strategy, 1, MPI_Integer, MASTER, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(isres_alpha, 1, MPI_Double_Precision, MASTER, MPI_COMM_WORLD, ierr)
#endif

    RETURN

END SUBROUTINE read_input
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
SUBROUTINE check_restart(yesorno)
!  =========================================================================
! |   Check to see if a restart file exists.
!  =========================================================================

    implicit none

    logical, intent(OUT) :: yesorno

    if (iam == MASTER) then
        inquire(file=restfname, exist=yesorno)
    end if

#ifdef HAVE_MPI
    call MPI_Bcast(yesorno, 1, MPI_LOGICAL, MASTER, MPI_COMM_WORLD, ierr)
#endif

    RETURN

END SUBROUTINE check_restart
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
SUBROUTINE print_input()

    USE input_vars
    USE timer
    USE model

    implicit none

! Get the time.
    call date_and_time(VALUES=start_time)
    last_restart = start_time

! Print a disclaimer.
    print *, ""
    print *, "ERGO, Copyright (C) 2007 Klaus Keller and Brian Tuttle"
    print *, "ERGO comes with ABSOLUTELY NO WARRANTY.  This is free"
    print *, "software, and you are welcome to redistribute it under certain"
    print *, "conditions.  See file GPL.txt that came with this distribution."
    print *, ""
    print *, ""
    write(*,"(a,i2,a,a,i5,i3,a,i2,a,i2,a,i3)") "Start time: ", &
          start_time(3), " ", month(start_time(2)), start_time(1), &
          start_time(5), ":", start_time(6), ":", start_time(7), ".", &
          start_time(8)

! Print what was just read.
    print *, ""
    print *, "Run time (wallclock):"
    print *, "   Precalculation time             = ", rtimePRE_str
    print *, "   Stage 1 time                    = ", rtimeST1_str
    print *, "   Stage 2 time                    = ", rtimeST2_str
    print *, "   Stage 3 time                    = ", rtimeST3_str
    print *, "Run time (iterations):"
    print *, "   Maximum Precalc. iterations     = ", maxiter_PRE
    print *, "   Maximum Stage 1 iterations      = ", maxiter_ST1
    print *, "   Maximum Stage 2 iterations      = ", maxiter_ST2
    print *, "   Maximum Stage 3 iterations      = ", maxiter_ST3
    print *, "Optimization strategies:"
    print *, "   Precalculation                  = ", opt_strategy_pre
    print *, "   Model optimization              = ", opt_strategy
    print *, "Convergence criteria:"
    print *, "   Target RMSE (PRE)      [%]      = ", RMStgt_pre
    print *, "   Target RMSE (DE)       [%]      = ", RMStgt_de
    print *, "   Target RMSE (SRES)     [%]      = ", RMStgt_sres
    print *, "   Target RMSE (ISRES)    [%]      = ", RMStgt_isres
    print *, "Thresholds for 3-chain strategy:"
    print *, "   RMSE between Q2 & Q3 [%]        = ", RMSrej_Q23
    print *, "   Weighted Av. delta RMSE         = ", THdRMSE
    print *, "   Weighted Av. delta Objf         = ", THdObjF
#ifdef ERGO
    print *, "DE parameters:"
    print *, "   population factor               = ", DE_pop
    print *, "   Random scaling factor           = ", DE_rndscale
    print *, "   mutation strategy [1-6]         = ", DE_strat
#endif
    print *, "ISERS parameters:"
    print *, "   alpha                           = ", isres_alpha
    print *, "Output files:"
    print *, "   fbest                           = ", trim(out_fbest)
    print *, "   xbest                           = ", trim(out_xbest)
    if (restarting) then
       print *, "Restarting from file               = ", restfname
    else
       print *, "No restart file found."
    end if

    call print_model_input()

    RETURN

END SUBROUTINE print_input
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
SUBROUTINE init_outfiles()
!  ========================================================================
! |  Initialize the fbest file for objective function output.              |
!  ========================================================================

    USE input_vars
    USE timer
    USE model

    implicit none

    integer, parameter :: iof = 60
    integer, parameter :: ioc = 70
!    character(len=4) :: numstr
!    character(len=11) :: fmtstr_yr
!    integer(i4b), dimension(n_param) :: years
    integer :: istat !, tx

! Create fbest output file.
    istat = 0
!    if ( .not. restarting ) then
        print *, stringcat("Creating ", stringcat(trim(out_fbest), "."))
        open(UNIT=iof, FILE=trim(out_fbest), STATUS='REPLACE', &
                                ACCESS='SEQUENTIAL', IOSTAT=istat)
        write(iof,"(a)") "Iteration     Objf evals    Objf value"
!    else
!        print *, stringcat("Appending to ", stringcat(trim(out_fbest), "."))
!        open(UNIT=iof, FILE=trim(out_fbest), STATUS='OLD', ACCESS='SEQUENTIAL',&
!            ACTION='WRITE', POSITION='APPEND', IOSTAT=istat)
!        write(iof,"(a,i2,a,a,i5,i3,a,i2,a,i2,a,i3)") "Restart time: ", &
!          start_time(3), " ", month(start_time(2)), start_time(1), &
!          start_time(5), ":", start_time(6), ":", start_time(7), ".", &
!          start_time(8)
!    end if
    close(UNIT=iof)
    if (istat > 0) then
        print *, "Error: Problem opening file ", trim(out_fbest)
        STOP
    end if

! Create convergence output file.
!    if ( writing_conv .and. .not. restarting ) then
      istat = 0
      print *, stringcat("Creating ", stringcat(trim(out_conv), "."))
      open(UNIT=ioc, FILE=trim(out_conv), STATUS='REPLACE', &
                          ACCESS='SEQUENTIAL', IOSTAT=istat)
      write(ioc,"(a)") "RMS error"
      write(ioc,"(a)") "Iteration     Objf evals    Time(minutes)        RMSE"
      close(UNIT=ioc)
      if (istat > 0) then
          print *, "Error: Problem opening file ", trim(out_conv)
          STOP
      end if
!    end if

! Initialize any output files for the model.
    call model_init_outfiles()

    RETURN

END SUBROUTINE init_outfiles
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
SUBROUTINE write_fbest(iter, nfe, A, B, E)
!  ========================================================================
! |  This subroutine appends a record of updated fbest values to the file  |
! |  of objective function output.                                         |
!  ========================================================================

    implicit none

    integer, intent(IN) :: iter
    integer(i8b), intent(IN) :: nfe
    real(dp), intent(IN) :: A, B, E
    integer, parameter :: iou = 60
    integer :: istat

    open(UNIT=iou, FILE=trim(out_fbest), STATUS='OLD', &
         ACCESS='SEQUENTIAL', ACTION='WRITE', POSITION='APPEND', IOSTAT=istat)

!   write(iou,"(i12,3(a,F22.11))") iter, "  A:", A, "  B: ", B, " Err: ", E
!    write(iou,"(i12,2(a,F17.11),a,F15.11)") iter, "  A: ", -A, "  B: ", -B, &
!                                                 " Err:", E
!    write(iou,"(i12,es15.6,F22.11)") iter, real(nfe), A
    write(iou,"(i12,i15,3(F19.6,1x))") iter, nfe, A, B, E


    close(UNIT=iou)

    RETURN

END SUBROUTINE write_fbest
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
SUBROUTINE write_endtime()
!  ========================================================================
! |  This subroutine writes the time at the end of the file of objective   |
! |  function output.  Eventually it might write the elapsed time.         |
!  ========================================================================

    USE timer

    implicit none

    integer, dimension(4) :: elapsed_time
!    integer, parameter :: iou = 60
!    integer :: istat

! Get the time.
    call date_and_time(VALUES=end_time)

!! Open the file.
!    open(UNIT=iou, FILE=trim(out_fbest), STATUS='OLD', &
!         ACCESS='SEQUENTIAL', ACTION='WRITE', POSITION='APPEND', IOSTAT=istat)

! Write the ending time.
    print *, ""
    write(*,"(a,i2,a,a,i5,i3,a,i2,a,i2,a,i3)") "End time: ", end_time(3), &
     " ", month(end_time(2)), end_time(1), end_time(5), ":", end_time(6), &
     ":", end_time(7), ".", end_time(8)
! Write the elapsed time.
    call time_diff(time_0, end_time, COMPONENTS = elapsed_time)
    write(*, "(a,i3,a,i2,a,i2,a,i2,a)") "Total time elapsed: ", &
                elapsed_time(1), " days, ", elapsed_time(2), " hours, ", &
                elapsed_time(3), " minutes, and ", elapsed_time(4), " seconds."

!! Done.
!    close(UNIT=iou)

    RETURN

END SUBROUTINE write_endtime
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
SUBROUTINE write_xout(num, xout1, xout2, f1, f2, iters, rmserr)
!  ========================================================================
! |  This subroutine writes the arrays of real numbers "xout1" and "xout2" |
! |  to a file with the name constructed from the basename 'out_xbest' and |
! |  the input integer "num."  Xout1 and xout2 must be the same length.    |
!  ========================================================================

    implicit none

    integer(i4b), intent(IN) :: iters
    integer(i4b), intent(IN) :: num
    real(DP), dimension(:), intent(IN) :: xout1, xout2
    real(DP), intent(IN) :: f1, f2, rmserr
    character(len=len(out_xbest)) :: xfname
    character(len=10) :: timestr
    character(len=0) :: nullstr
    integer :: iounit = 40
    integer :: k, istat, datalen

! Check that the data lengths are the same.
    datalen = size(xout1,1)
    if (size(xout2,1) .ne. datalen) then
        print *, "write_xout: xout1 and xout2 are different sizes"
        RETURN
    end if

! Name the file.
!    xfname = filenamer(trim(out_xbest), num, 2, nullstr)
    xfname = trim(out_xbest)
! Open the file.
    open (UNIT=iounit, FILE=xfname, STATUS='REPLACE', &
                ACCESS='SEQUENTIAL', IOSTAT=istat)
! Write the header with time and iteration.
    call date_and_time(time=timestr)
    write (iounit,"(a)") stringcat("Wall time [hhmmss.sss]: ", timestr)
    write (iounit,"(a, i8)") "Iteration: ", iters
    write (iounit,"(a, f8.5)") "RMSerror: ", rmserr
    write (iounit,"(a, f20.11)") "FBest_A: ", f1
    write (iounit,"(a, f20.11)") "FBest_B: ", f2
    write (iounit,"(a)") "   XBestA:             XBestB:"
! Write the data.
    do k = 1, datalen
      write(iounit,'(2F20.11)') xout1(k), xout2(k)
    end do
! Close the file.
    close (UNIT=iounit)

! Write the RMS error to its own file.
!    if (writing_conv) call write_conv(iters, rmserr)

    RETURN

END SUBROUTINE write_xout
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
SUBROUTINE write_conv(iter, nfe, Err)        !, Av_dRMSe, Av_dobjf)
!  ========================================================================
! |  This subroutine appends a record of updated rmse values to the file   |
! |  of convergence output.                                                |
!  ========================================================================

    USE timer

    implicit none

    Integer, intent(IN) :: iter
    Integer(i8b), intent(IN) :: nfe
    Real(DP), dimension(:), intent(IN) :: Err
!    Real(DP), dimension(:), intent(IN) :: Av_dRMSe
!    Real(DP), dimension(:), intent(IN) :: Av_dobjf
    Integer(i4b), dimension(8) :: time_now
    Real(DP) :: hours_since_start, minutes
    Integer, parameter :: iou = 90
    Integer :: istat, errlen   !, oflen
    Character(len=19) :: fmtstr

! Find sizes.
    errlen = size(Err,1)
!    oflen = size(Av_dobjf,1)

! Get the elapsed time.
    call date_and_time(VALUES=time_now)
    call time_diff(start_time, time_now, HOURS = hours_since_start)
    minutes = (hours_since_start) * 60

    open(UNIT=iou, FILE=trim(out_conv), STATUS='OLD', &
         ACCESS='SEQUENTIAL', ACTION='WRITE', POSITION='APPEND', IOSTAT=istat)

!    write(fmtstr,"(a,i1,a)") "(i12,F10.2,",2*errlen+oflen,"F22.11)"
!    write(iou,fmtstr) iter, minutes, Err, Av_dRMSe, Av_dobjf
!    write(iou,"(i12,es15.6,F10.2,F22.11)") iter, nfe, minutes, Err(1)
    write(iou,"(i12,i15,F10.2,F22.11)") iter, nfe, minutes, Err(1)

    close(UNIT=iou)

    RETURN

END SUBROUTINE write_conv
!---------------------------------------------------------------------------

!!---------------------------------------------------------------------------
!SUBROUTINE write_restart(n, lambda, iter, filecnt, x_lower, x_upper, sigmaA, &
!                         sigmaB, indA, indB, xA, xB, xbestA, xbestB, fbestA, &
!                         fbestB, lastfbest, feasible)
!!  ========================================================================
!! |  This subroutine writes a restart file as Fortran unformatted binary.  |
!! |  The file was named restart above in the read_input subroutine and     |
!! |  be replaced each time a restart file is written, keeping only the     |
!! |  most recent state of the optimization.                                |
!!  ========================================================================
!
!    USE timer
!
!    implicit none
!
!    integer(i4b), intent(IN) :: n, lambda
!    integer(i4b), intent(IN) :: iter, filecnt, lastfbest
!    real(DP), dimension(n), intent(IN) :: x_lower
!    real(DP), dimension(n), intent(IN) :: x_upper
!    real(DP), dimension(n,lambda), intent(IN) :: sigmaA, sigmaB, xA, xB
!    integer(i4b), dimension(lambda), intent(IN) :: indA, indB
!    real(DP), dimension(n), intent(IN) :: xbestA, xbestB
!    real(DP), intent(IN) :: fbestA, fbestB
!    logical, intent(IN) :: feasible
!    integer, parameter :: iou = 60
!    integer :: istat, fnumlen, i, j
!    integer(i4b), dimension(8) :: timenow
!    real(DP) :: time_since_start, total_run_time
!    character(len=len(out_restart)) :: fname
!
!    print *, "Writing a restart file... "
!! Get the time.
!    call date_and_time(values = timenow)
!    call time_diff(start_time, timenow, hours = time_since_start)
!    total_run_time = time_since_start
!! Make a filename.
!    fname = trim(out_restart)
!! Open the file.
!    open (UNIT=iou, FILE=fname, STATUS='REPLACE', ACTION='WRITE', &
!            FORM='UNFORMATTED', IOSTAT=istat)
!! Write data.
!    write(IOU) iter
!    write(IOU) filecnt
!    do i = 1, n
!        write(IOU) x_lower(i)
!        write(IOU) x_upper(i)
!    end do
!    do i = 1, n
!      do j = 1, lambda
!        write(IOU) xA(i,j)
!        write(IOU) xB(i,j)
!        write(IOU) sigmaA(i,j)
!        write(IOU) sigmaB(i,j)
!      end do
!    end do
!    do j = 1, lambda
!        write(IOU) indA(j)
!        write(IOU) indB(j)
!    end do
!    do i = 1, n
!        write(IOU) xbestA
!        write(IOU) xbestB
!    end do
!    write(IOU) fbestA
!    write(IOU) fbestB
!    write(IOU) lastfbest
!    write(IOU) feasible
!    write(IOU) total_run_time
! 
!! Close the file.
!    close (UNIT=iou)
!
!    print *, "                            done."
!
!    RETURN
!
!END SUBROUTINE write_restart
!!---------------------------------------------------------------------------

!!---------------------------------------------------------------------------
!SUBROUTINE read_restart(n, lambda, iter, filecnt, x_lower, x_upper, sigmaA, &
!                        sigmaB, indA, indB, xA, xB, xbestA, xbestB, fbestA, &
!                        fbestB, lastfbest, feasible)
!!  ========================================================================
!! |  This subroutine checks to see that a restart file was specified and   |
!! |  proceeds to read the Fortran unformatted binary file, assigns values, |
!! |  and transmits them to the other processors.                           |
!!  ========================================================================
!
!    USE timer
!
!    implicit none
!
!    integer(i4b), intent(IN) :: n, lambda
!    integer(i4b), intent(OUT) :: iter, filecnt, lastfbest
!    real(DP), dimension(n), intent(OUT) :: x_lower
!    real(DP), dimension(n), intent(OUT) :: x_upper
!    real(DP), dimension(n,lambda), intent(OUT) :: sigmaA, sigmaB, xA, xB
!    integer(i4b), dimension(lambda), intent(OUT) :: indA, indB
!    real(DP), dimension(n), intent(OUT) :: xbestA, xbestB
!    real(DP), intent(OUT) :: fbestA, fbestB
!    logical, intent(OUT) :: feasible
!    integer, parameter :: IOU = 70
!    integer :: istat, i, j
!
!! Open the file.
!    open (UNIT=IOU, FILE=restfname, STATUS='OLD', ACTION='READ', &
!            FORM='UNFORMATTED', IOSTAT=istat)
!! Read data.
!    read (IOU) iter
!    read (IOU) filecnt
!    do i = 1, n
!        read(IOU) x_lower(i)
!        read(IOU) x_upper(i)
!    end do
!    do i = 1, n
!      do j = 1, lambda
!        read(IOU) xA(i,j)
!        read(IOU) xB(i,j)
!        read(IOU) sigmaA(i,j)
!        read(IOU) sigmaB(i,j)
!      end do
!    end do
!    do j = 1, lambda
!        read(IOU) indA(j)
!        read(IOU) indB(j)
!    end do
!    do i = 1, n
!        read(IOU) xbestA
!        read(IOU) xbestB
!    end do
!    read(IOU) fbestA
!    read(IOU) fbestB
!    read(IOU) lastfbest
!    read(IOU) feasible
! 
!
!! Close the file.
!    close (UNIT=IOU)
!
!    RETURN
!
!END SUBROUTINE read_restart
!!---------------------------------------------------------------------------

!!---------------------------------------------------------------------------
!SUBROUTINE read_startpol(startpol, nop, n_sow)
!!  ========================================================================
!! |  This routine opens a starting policy text file containing the         |
!! |  abatement (mu) variable formatted the same as dvar_mu output file.    |
!! |  Output is an array nop by n_sow containing the data.                  |
!!  ========================================================================
!
!    implicit none
!
!    integer(i4b), intent(IN) :: nop, n_sow
!    real(DP), dimension(nop,n_sow), intent(OUT) :: startpol
!    integer, parameter :: IOU = 75
!    character(len=4) :: numstr, header
!    character(len=14) :: fmtstr
!    integer :: istat, tx, year
!
!! Make a format string.
!    write(numstr,"(i4)") n_sow
!    fmtstr = stringcat("(i4,", stringcat(trim(numstr), "F12.6)"))
!! Open the file.
!    open (UNIT=IOU, FILE=stpolfname, STATUS='OLD', ACTION='READ', &
!            FORM='FORMATTED', IOSTAT=istat)
!! Read data.
!    call skipline(IOU, header, 1)
!    do tx = 1, nop
!        read (IOU, trim(fmtstr)) year, startpol(tx,:)
!    end do
!
!! Close the file.
!    close (UNIT=IOU)
!
!    RETURN
!
!END SUBROUTINE read_startpol
!!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
SUBROUTINE write_in_arch()
!  ========================================================================
! |  Writes an input file that will recreate the current model run.        |
!  ========================================================================

    USE input_vars
    USE model

    implicit none

    integer(i4b), dimension(8) :: timestamp
    integer(i4b), parameter :: IOU = 80
    integer(i4b) :: istat

! Get the date and time.
    call date_and_time(VALUES=timestamp)

! Open the file.
    open (UNIT=IOU, FILE=inarchfname, STATUS='REPLACE', ACTION='WRITE', &
            FORM='FORMATTED', IOSTAT=istat)
! Write.
write(IOU,"(a)") &
" ==========================================================================="
write(IOU,"(a)") &
"|  Input archive file for FRANC optimization using ERGO                     |"
write(IOU,"(a)") &
"|  automatically generated by mod_fileIO.f90: SUBROUTINE write_in_arch()    |"
write(IOU,"(a,i2,3a,i4,a,i2,a,i2,a)")  "|  ", timestamp(3), &
    " ", month(timestamp(2)), " ", timestamp(1), " at ", timestamp(5), &
    ":", timestamp(6), "                                                     |"
write(IOU,"(a)") &
"|  ** Changes to the format of this file must be reflected in               |"
write(IOU,"(a)") &
"|  ** mod_fileIO.f90: read_input()                                          |"
write(IOU,"(a)") &
" ==========================================================================="
write(IOU,"(a)") &
"------------------------ ERGO PARAMETERS ------------------------------------"
write(IOU,"(a)") "Run time (wallclock time) [hhh:mm]:"
write(IOU,"(2a)") "    Precalculation           = ", rtimePRE_str
write(IOU,"(2a)") "    Stage 1                  = ", rtimeST1_str
write(IOU,"(2a)") "    Stage 2                  = ", rtimeST2_str
write(IOU,"(2a)") "    Stage 3                  = ", rtimeST3_str
write(IOU,"(a)") "Run time (iterations):"
write(IOU,"(a)") " A value greater than 0 overrides wallclock spec."
write(IOU,"(a,i12)") "    Maximum iterations pre   = ", maxiter_pre
write(IOU,"(a,i12)") "    Maximum iterations st1   = ", maxiter_ST1
write(IOU,"(a,i12)") "    Maximum iterations st2   = ", maxiter_ST2
write(IOU,"(a,i12)") "    Maximum iterations st3   = ", maxiter_ST3
write(IOU,"(a)") " "
write(IOU,"(a)") "Optimization strategies:"
write(IOU,"(a,i1)") "    Precalculation (0-2)     = ", opt_strategy_pre
write(IOU,"(a,i1)") "    Model optimization (0-8) = ", opt_strategy
write(IOU,"(a)") &
" Strategy choices (number of chains) (Stage 1 [-> Stage 2 [-> Stage 3]]):"
write(IOU,"(a)") &
"    1. DE(2)    2. SRES(2)  3. DE(3)    4. SRES(3)  5. DE(2) -> SRES(2)"
write(IOU,"(a)") "    6. SRES(2) -> ISRES(2)      7. SRES(3) -> ISRES(2)"
write(IOU,"(a)") &
"    8. DE(3) -> SRES(3) -> ISRES(2)     0. Skip optimization step"
write(IOU,"(a)") " "
write(IOU,"(a)") "Convergence criteria:"
write(IOU,"(a,f17.15)") "    Target RMSE (PRE)   [%]  = ", RMStgt_pre
write(IOU,"(a,f17.15)") "    Target RMSE (DE)    [%]  = ", RMStgt_de
write(IOU,"(a,f17.15)") "    Target RMSE (SRES)  [%]  = ", RMStgt_sres
write(IOU,"(a,f17.15)") "    Target RMSE (ISRES) [%]  = ", RMStgt_isres
write(IOU,"(a)") " "
write(IOU,"(a)") "Thresholds for 3-chain strategy:"
write(IOU,"(a,f17.15)") "    RMSE between Q2 & Q3 [%] = ", RMSrej_Q23
write(IOU,"(a,f17.15)") "    Weighted Av. delta RMSE  = ", THdRMSE
write(IOU,"(a,f17.15)") "    Weighted Av. delta Objf  = ", THdObjF
write(IOU,"(a)") " "
#ifdef ERGO
write(IOU,"(a)") "DE parameters:"
write(IOU,"(a,i3)") "    population factor        = ", DE_pop
write(IOU,"(a,i1)") "    Random scaling factor    = ", DE_rndscale
write(IOU,"(a,i1)") "    mutation strategy [1-6]  = ", DE_strat
write(IOU,"(a)") " "
#endif
write(IOU,"(a)") "ISRES parameters:"
write(IOU,"(a,f7.4)") "    alpha [default 0.2]      = ", isres_alpha
write(IOU,"(a)") " "
write(IOU,"(a)") "Output files (max 40 characters):"
write(IOU,"(2a)") "    Path                     = ", out_path
write(IOU,"(2a)") "    fbest                    = ", fbestname
write(IOU,"(2a)") "    xbest (set of files)     = ", xbestname
write(IOU,"(2a)") "    convergence (optional)   = ", convfname
write(IOU,"(a)") " "
write(IOU,"(a)") &
"------------------------ MODEL PARAMETERS ---------------------------------"

call model_write_arch(IOU)

write(IOU,"(a)") " "
write(IOU,"(a)") &
"---------------------- END MODEL PARAMETERS -------------------------------"

! Close the file.
    close (UNIT=IOU)

    RETURN

END SUBROUTINE write_in_arch
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
FUNCTION filenamer(fileset, indx, ind_len, separator) RESULT(filename)
!  ========================================================================
! |  This function generates a file name from a root name called           |
! |  "fileset" and an two-digit integer.                                   |
! |  It translates the integer into a character string and                 |
! |  appends it to the end of the name.                                    |
!  ========================================================================

    implicit none

    character(len=*) :: fileset
    character(len=*) :: separator
    integer :: ind_len
    character(len=len(fileset)+ind_len+len(separator)) :: filename
    character(len=ind_len) :: numstr
    integer :: indx

    numstr = int2str(indx,ind_len)

    filename = stringcat(trim(fileset), stringcat(separator, numstr))

END FUNCTION filenamer
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
FUNCTION int2str(intnumber, length) RESULT(strnumber)
!  ========================================================================
! |  This function converts an integer into a string representing          |
! |  that same number, but preceded with zeros to fill the                 |
! |  length specified in the input.                                        |
!  ========================================================================

    implicit none

    integer :: intnumber, length, numlen
    character(len=length) :: strnumber
    character(len=4) :: lenstr

    if (intnumber > 1000000000) then
        print *, "Input range of function exceeded: int2str(intnumber, length)"
        write (*,"(a,i10)") "intnumber = ", intnumber
        STOP
    end if

    if (intnumber == 0) then
        numlen = 1
    else
        numlen = int(log10(real(intnumber))) + 1
    end if
    write(lenstr,"(a,i1,a)") "(i", numlen, ")"
    write(strnumber,lenstr) intnumber
    do
        if (len(trim(strnumber)) == length) then
            exit
        else
            strnumber = stringcat('0', trim(strnumber))
        end if
    end do

END FUNCTION int2str
!---------------------------------------------------------------------------

END MODULE fileIO
