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

MODULE converge

    USE global

    implicit none

    Private

    integer(i4b), public :: NQ      ! number of concurrent solutions
    integer(i8b), public :: nfe_total ! number of function evals to report
    integer(i4b), public :: cx_len  ! length of convergence index vector
    integer(i4b), dimension(:), allocatable, public :: conv_indx
    real(DP) :: RMScrit, rms_norm

!    integer :: iter_newQ3 
!    real(DP) :: rmse_newQ3
!    real(DP) :: bestval_newQ3
    integer, dimension(:), allocatable :: last_Qrank, iter_newQ
    real(DP), dimension(:), allocatable :: last_rmse, prev_best
    real(DP), dimension(:), allocatable :: WArmsnum_prev, WArmsden_prev
    real(DP), dimension(:), allocatable :: WAobfnum_prev, WAobfden_prev

    public :: init_conv_check, fin_conv_check, conv_check !, init_conv 
    public :: sort_val, dealloc_conv_indx, alloc_conv_indx

CONTAINS
!---------------------------------------------------------------------------
SUBROUTINE init_conv_check(rmstgt, numberQs, lower_limit, upper_limit)
!  =========================================================================
! |  This subroutine sets up global variables and allocates size NQ         |
! |  variables to accomodate running with 1, 2 or 3 convergence queues.     |
!  =========================================================================

    implicit none

    real(DP), intent(IN) :: rmstgt
    integer(i4b), intent(IN) :: numberQs
    real(DP), dimension(:), intent(IN) :: lower_limit, upper_limit
    integer :: astat

    NQ = min(numberQs, 3)       ! Maximum number of queues is 3.

! Allocate NQ-dependent global variables.
    astat = 0
    allocate(last_Qrank(NQ), STAT=astat)
    if (astat > 0) call error_msg(1,"last_Qrank")

    astat = 0
    allocate(last_rmse(NQ), STAT=astat)
    if (astat > 0) call error_msg(1,"last_rmse")

    astat = 0
    allocate(prev_best(NQ), STAT=astat)
    if (astat > 0) call error_msg(1,"prev_best")

    astat = 0
    allocate(WArmsnum_prev(NQ), STAT=astat)
    if (astat > 0) call error_msg(1,"WArmsnum_prev")

    astat = 0
    allocate(WArmsden_prev(NQ), STAT=astat)
    if (astat > 0) call error_msg(1,"WArmsden_prev")

    astat = 0
    allocate(WAobfnum_prev(NQ), STAT=astat)
    if (astat > 0) call error_msg(1,"WAobfnum_prev")

    astat = 0
    allocate(WAobfden_prev(NQ), STAT=astat)
    if (astat > 0) call error_msg(1,"WAobfden_prev")

    astat = 0
    allocate(iter_newQ(NQ), STAT=astat)
    if (astat > 0) call error_msg(1,"iter_newQ")

! Define convergence criterion.
    RMScrit = rmstgt

! Define normalizing factor for RMS error measurement.
    rms_norm = sqrt(sum((upper_limit - lower_limit)**2)/real(size(upper_limit)))

! Initialize convergence records.
!    bestval_newQ3 = 0.0d0
!    rmse_newQ3 = 100.0d0
!    iter_newQ3 = 0
    iter_newQ = 0
    last_Qrank = 1    
    last_rmse = 100.0d0
    prev_best = 0.0d0
    WArmsnum_prev = 0
    WArmsnum_prev = 0
    WAobfnum_prev = 0
    WAobfnum_prev = 0

END SUBROUTINE init_conv_check
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
SUBROUTINE conv_check(iters,intvl,bestval,best_array, Q_rank, stop_flg, &
                        action_flg)
!  =========================================================================
! |  This routine decides what actions to take by measuring the convergence |
! |  of the three concurrent solutions.  It returns flags for ending the    |
! |  optimization and for restarting the worst solution as well as the      |
! |  indices of Q1 and Q3.                                                  |
!  =========================================================================

    USE fileIO
    USE input_vars
    USE model

    implicit none

    integer i
    integer(i4b), intent(IN) :: iters, intvl
    real(DP), dimension(NQ), intent(IN) :: bestval
    real(DP), dimension(:,:), intent(IN) :: best_array
    integer(i4b), dimension(NQ), intent(OUT) :: Q_rank
    logical, intent(OUT) :: stop_flg, action_flg
    real(DP) :: Q3_metric
    real(DP), dimension(NQ) :: rmserr, delta_rmse, delta_objf
    real(DP), dimension(NQ) :: WAnum, WAdenom, WAdrmse, WAdobjf

! Initialize flags.
    stop_flg = .false.
    action_flg = .false.

  if (NQ > 2) then      !! 3 Queues
! Sort the queues in order of objective function values.
    Q_rank = sort_val(bestval)
!    print *, "Iter", iters, "bestval", bestval, "Q_rank", Q_rank

! Measure the RMS error between the 1st and 2nd, and the 1st and 3rd queues.
    rmserr(1) = rms_error(best_array(:,Q_rank(1)),best_array(:,Q_rank(2)))
    rmserr(2) = rms_error(best_array(:,Q_rank(2)),best_array(:,Q_rank(3)))
    rmserr(3) = rms_error(best_array(:,Q_rank(1)),best_array(:,Q_rank(3)))
!    print *, "rmserr: ", rmserr

! Find the change in rmserr.
    delta_rmse = (rmserr - last_rmse) / intvl
! Find the change in objective function.
    delta_objf(1) = (bestval(Q_rank(1)) - prev_best(1)) / intvl
    delta_objf(2) = (bestval(Q_rank(2)) - prev_best(2)) / intvl
    delta_objf(3) = (bestval(Q_rank(3)) - prev_best(3)) / intvl
! Calculate the weighted running average of delta_rmse.
    WAnum = WArmsnum_prev + (iters - iter_newQ)*delta_rmse
    WAdenom = WArmsden_prev + (iters - iter_newQ)
    WAdrmse = WAnum / WAdenom
    WArmsnum_prev = WAnum
    WArmsden_prev = WAdenom
! Calculate the weighted running average of delta_objf.
    WAnum = WAobfnum_prev + (iters - iter_newQ)*delta_objf
    WAdenom = WAobfden_prev + (iters - iter_newQ)
    WAdobjf = WAnum / WAdenom
    WAobfnum_prev = WAnum
    WAobfden_prev = WAdenom
! Determine convergence of 2 best queues.
!    print *, "iters = ", iters
!    print *, "rmserr(1) = ", rmserr(1)
!    print *, "RMScrit = ", RMScrit
!    print *, "Weighted average drmse for stopflag ", abs(WAdrmse(1))
!    print *, "Weighted average drmse for restart Q3", abs(WAdrmse(3))
!    print *, "Threshold for dRMSE = ", THdRMSE
!    print *, "Weighted average dobjf for restart Q3", abs(WAdobjf(3))
!    print *, "Threshold for dObjF = ", THdObjF
!    print *, "WA dobjf(3) / rmserr(3) = ", abs(WAdobjf(3)/rmserr(3))
    if (rmserr(1) < RMScrit .and. abs(WAdrmse(1)) < THdRMSE) then
        stop_flg = .true.
! Determine whether to restart the 3rd best queue.
    else if (rmserr(2) < min(RMSrej_Q23, rmserr(1))) then 
!        print *, "Solution 3 is near solution 2."
        action_flg = .true.
!    else if (abs(WAdrmse(3)) < THdRMSE .and. abs(WAdobjf(3)) < THdObjF) then
!    else if (abs(WAdrmse(3) * WAdobjf(3)) < 1.0) then
    else if (abs(WAdobjf(3)/rmserr(3)) < 1.0) then
!        print *, "Solution 3 is stalled."
        action_flg = .true.
    end if
! Check whether the RMS error has changed since last time.
    if (abs(sum(rmserr - last_rmse)) > 0.00000000001 .or. rmserr(1) .ne. rmserr(1))  then
        call write_fbest(iters, nfe_total, bestval(Q_rank(1)), &
                        bestval(Q_rank(2)), rmserr(1))
!        call write_fbest(iters, bestval(Q_rank(2)), &
!                            bestval(Q_rank(3)), rmserr(2))
!        call write_fbest(iters, bestval(Q_rank(1)), &
!                            bestval(Q_rank(3)), rmserr(3))
        !call model_out(best_array(:,Q_rank(1)), out_path)
        open(UNIT=43, FILE='./solution_out', STATUS='REPLACE', ACTION='WRITE', &
             ACCESS='SEQUENTIAL')

        write(43,"(a,g12.5)") '#  f(x) = ', bestval(Q_rank(1))
        write(43,"(a)") '# i    x(i)'
        do i = 1,size(best_array(:,Q_rank(1)))
           write(43, "(i3,1x,g17.10)") i, best_array(i,Q_rank(1))
        end do
        close(43)
        call write_conv(iters, nfe_total, rmserr)      !, WAdrmse, WAdobjf)


!    ! Determine convergence of 2 best queues.
!        print *, "iters = ", iters
!        print *, "rmserr(1) = ", rmserr(1)
!        print *, "RMScrit = ", RMScrit
!        print *, "abs(WAdrmse(1)) = ", abs(WAdrmse(1))
!        print *, "THdRMSE = ", THdRMSE
!        if (rmserr(1) < RMScrit .and. abs(WAdrmse(1)) < THdRMSE) then
!            stop_flg = .true.
!    ! Determine whether to restart the 3rd best queue.
!        else if (rmserr(2) < min(RMSrej_Q23, rmserr(1))) then 
!!            print *, "Solution 3 is near solution 2."
!            action_flg = .true.
!        else if (abs(WAdrmse(3)) < THdRMSE .and. abs(WAdobjf(3)) < THdObjF) then
!!            print *, "Solution 3 is stalled."
!            action_flg = .true.
!        end if
    ! Save the error for next time.
        last_rmse = rmserr
        prev_best(1) = bestval(Q_rank(1))
        prev_best(2) = bestval(Q_rank(2))
        prev_best(3) = bestval(Q_rank(3))
!        if (action_flg) then
!            iter_newQ(3) = iters
!            WArmsnum_prev(2:3) = 0
!            WArmsden_prev(2:3) = 0
!            last_rmse(3) = 100.0d0
!            WAobfnum_prev(2:3) = 0
!            WAobfden_prev(2:3) = 0
!            prev_best(3) = 0.0d0
!        end if
    end if
    if (action_flg) then
        iter_newQ(3) = iters
        WArmsnum_prev(2:3) = 0
        WArmsden_prev(2:3) = 0
        last_rmse(3) = 100.0d0
        WAobfnum_prev(2:3) = 0
        WAobfden_prev(2:3) = 0
        prev_best(3) = 0.0d0
    end if

  else      !! 2 Queues
! Sort the queues in order of objective function values.
    Q_rank = sort_val(bestval)
!    print *, "Iter", iters, "bestval", bestval, "Q_rank", Q_rank

! Measure the RMS error between the 1st and 2nd queues.
    rmserr(1) = rms_error(best_array(:,Q_rank(1)),best_array(:,Q_rank(2)))
!    print *, "rmserr(1) = ", rmserr(1)

! Find the change in rmserr.
    delta_rmse = (rmserr - last_rmse) / intvl
! Find the change in objective function.
    delta_objf(1) = (bestval(Q_rank(1)) - prev_best(1)) / intvl
! Calculate the weighted running average of delta_rmse.
    WAnum = WArmsnum_prev + (iters - iter_newQ)*delta_rmse
    WAdenom = WArmsden_prev + (iters - iter_newQ)
    WAdrmse = WAnum / WAdenom
    WArmsnum_prev = WAnum
    WArmsden_prev = WAdenom
!! Calculate the weighted running average of delta_objf.
!    WAnum = WAobfnum_prev + (iters - iter_newQ)*delta_objf
!    WAdenom = WAobfden_prev + (iters - iter_newQ)
!    WAdobjf = WAnum / WAdenom
!    WAobfnum_prev = WAnum
!    WAobfden_prev = WAdenom
! Determine convergence of 2 queues.
    !print*, 'test1'
    if (rmserr(1) < RMScrit .and. abs(WAdrmse(1)) < THdRMSE) then
        stop_flg = .true.
    end if
! Check whether the RMS error has changed since last time.
    if (sum(rmserr - last_rmse) .ne. 0 ) then
        call write_fbest(iters, nfe_total, bestval(Q_rank(1)), &
                        bestval(Q_rank(2)), rmserr(1))
        call model_out(best_array(:,Q_rank(1)), out_path)
        call write_conv(iters, nfe_total, rmserr)  !, WAdrmse, WAdobjf)

!    ! Determine convergence of 2 queues.
!        if (rmserr(1) < RMScrit .and. abs(WAdrmse(1)) < THdRMSE) then
!            stop_flg = .true.
!        end if
!   ! Determine whether convergence has stalled and is ready to switch to isres.
!        if (rmserr(1) < 1.0d0 .and. abs(WAdrmse(1)) < 0.005 &
!                              .and. abs(WAdobjf(1)) < 0.005) then
!            action_flg = .true.
!            iter_newQ(1) = iters
!            WArmsnum_prev(1) = 0
!            WArmsden_prev(1) = 0
!            WAobfnum_prev(1) = 0
!            WAobfden_prev(1) = 0
!        end if
    ! Save the error for next time.
        last_rmse = rmserr
        prev_best(1) = bestval(Q_rank(1))
    end if
    print*, 'test2'

  end if

! Save the ranking for next time.
    last_Qrank = Q_rank

    RETURN

END SUBROUTINE conv_check
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
SUBROUTINE fin_conv_check()
!  ========================================================================
! |  Deallocate global NQ-dependent variables.                             |
!  ========================================================================

    implicit none

    integer :: astat

    astat = 0
    deallocate(last_Qrank, STAT=astat)
    if (astat > 0) call error_msg(2,"last_Qrank")

    astat = 0
    deallocate(last_rmse, STAT=astat)
    if (astat > 0) call error_msg(2,"last_rmse")

    astat = 0
    deallocate(prev_best, STAT=astat)
    if (astat > 0) call error_msg(2,"prev_best")

    astat = 0
    deallocate(WArmsnum_prev, STAT=astat)
    if (astat > 0) call error_msg(2,"WArmsnum_prev")

    astat = 0
    deallocate(WArmsden_prev, STAT=astat)
    if (astat > 0) call error_msg(2,"WArmsden_prev")

    astat = 0
    deallocate(WAobfnum_prev, STAT=astat)
    if (astat > 0) call error_msg(2,"WAobfnum_prev")

    astat = 0
    deallocate(WAobfden_prev, STAT=astat)
    if (astat > 0) call error_msg(2,"WAobfden_prev")

    astat = 0
    deallocate(iter_newQ, STAT=astat)
    if (astat > 0) call error_msg(2,"iter_newQ")

    RETURN

END SUBROUTINE fin_conv_check
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
FUNCTION sort_val(values) RESULT(order)
!  ========================================================================
! |  Sort values in descending order and return their indices.  Copied     |
! |  and adapted from randperm() in mod_DE.f90.                            |
!  ========================================================================

    implicit none

    real(DP), dimension(:) :: values
    integer, dimension(size(values)) :: order
    integer :: i, j, k, indx, length

    length = size(values,1)
    do i=1,length
        indx = 1
        do j=1,length
            if (values(i) > values(j)) then
                indx = indx + 1
            end if
        end do
        do k=1,i-1
            if (values(i) <= values(k) .and. values(i) >= values(k)) then
                indx = indx + 1
            end if
        end do
        order(indx) = i
    end do


END FUNCTION sort_val
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
FUNCTION rms_error(vectorA, vectorB) RESULT(rmse)
!  ========================================================================
! |  Calculate the RMS error between mu and I/Q policy1 and policy2.       |
! |  Error is reported in percent.                                         |
!  ========================================================================

    real(DP), dimension(:) :: vectorA, vectorB
    real(DP) :: rmse
    real(DP), dimension(cx_len) :: diff
    integer(i4b) :: i

! Find the error between the two policies for the chosen elements.
    diff = vectorA(conv_indx) - vectorB(conv_indx)

! Compute the root mean square of the error.
    rmse = 100.0 * sqrt(sum(diff*diff)/real(cx_len)) / rms_norm

    RETURN

END FUNCTION rms_error
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
SUBROUTINE alloc_conv_indx()
!  ========================================================================
! |  Allocate memory for the convergence index vector.                     |
!  ========================================================================

    implicit none

    integer :: astat

    astat = 0
    allocate(conv_indx(cx_len), STAT=astat)
    if (astat > 0) call error_msg(1,"conv_indx")

    RETURN

END SUBROUTINE alloc_conv_indx
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
SUBROUTINE dealloc_conv_indx()
!  ========================================================================
! |  Deallocate memory for the convergence index vector.                   |
!  ========================================================================

    implicit none

    integer :: astat

    astat = 0
    deallocate(conv_indx, STAT=astat)
    if (astat > 0) call error_msg(2,"conv_indx")

    cx_len = 0

    RETURN

END SUBROUTINE dealloc_conv_indx
!---------------------------------------------------------------------------
END MODULE converge
