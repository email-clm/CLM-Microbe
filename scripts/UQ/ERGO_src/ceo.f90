!  CEO Constrained Evolutionary Optimization
!  Copyright (C) 2005 Thomas Philip Runarsson 
!
!  Fortran 90 implementation by Gudlaugur Johannesson and Thomas Runarsson
!
!  This program is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2 of the License, or
!  (at your option) any later version.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program; if not, write to the Free Software
!  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!  Please properly cite this paper for algorithm SRES:
!
!  Thomas Philip Runarsson and Xin Yao, Stochastic Ranking for 
!  Constrained Evolutionary Optimization. IEEE Transactions on 
!  Evolutionary Computation, Vol. 4, No. 3, September pp 274-283, 2000.
!
!  Please properly cite this paper for algorithm ISRES:
!
!  Thomas Philip Runarsson and Xin Yao, Search Biases in
!  Constrained Evolutionary Optimization. IEEE Transactions on
!  Systems, Man and Cybernetics Special Issue on Knowledge
!  Extraction and Incorporation in Evolutionary Computation (in press)
!---------------------------------------------------------------------------
!  Modified for use with DICE objective function by Brian Tuttle.
!  25 April 2005  btuttle@psu.edu
!---------------------------------------------------------------------------

MODULE ceo
!==============================================================================
! This module contains the subroutines sres and isres which are public
! Additionally there is a subroutine for stochastic ranking 
! There are also two functions randn() for generating Gaussian distributed
! random variables and a function that produces a random seed by reading
! the Linux random device.
!==============================================================================

    USE global
    USE model
    USE converge

    implicit none

    Private

! Define logical TRUE
    Integer, Parameter, Public :: lgt = Kind(.true.)

    Public :: run_sres, dealloc_sres
!    Public :: sres, run_isres, init_sres, dealloc_sres, extend

! Probability of ranking w.r.t. objective
    Real(dp), Parameter :: pf = 0.45d0    

    Integer(i4b), parameter :: refresh = 3 ! Iterations between progress chks
    Integer(i4b) :: n_x                      ! Number of variables
    Integer(i4b), parameter :: mu = 15  !15       ! Number of parents
    Integer(i4b), parameter :: lambda = 128   ! Number of individuals in 
                                             !   the population
    Real(dp) :: alpha                        ! ISRES randomization range 
                                             !   (originally 0.2)
    Real(dp), Parameter :: gamma = 0.85d0
    Real(dp) :: tau, taud, varphi 
    Real(dp), Dimension(:), Allocatable :: sigmaub
    Integer(i4b) :: total_iters
    Logical(lgt) :: isres

CONTAINS

!---------------------------------------------------------------------------
SUBROUTINE run_sres(solution, solution_sigma, max_iters, max_time, &
                    rms_target, numQs, INIT_POP, INIT_SIGMA, IN_ALPHA)
!  =========================================================================
! |
!  =========================================================================

    USE mt19937
    USE rndseed

    implicit none

#ifdef HAVE_MPI
  include "mpif.h"
#endif

    real(DP), dimension(:,:), intent(OUT) :: solution
    real(DP), dimension(:,:), intent(OUT) :: solution_sigma
    integer(i4b), intent(INOUT) :: max_iters
    real(DP), intent(IN) :: max_time
    real(DP), intent(IN) :: rms_target
    integer(i4b), intent(IN) :: numQs
    real(DP), dimension(:,:), intent(IN), optional :: INIT_POP
    real(DP), dimension(:,:), intent(IN), optional :: INIT_SIGMA
    real(DP), intent(IN), optional :: IN_ALPHA
    real(DP), dimension(size(solution,1), size(solution,2)) :: x_out
    real(DP), dimension(size(solution,1), size(solution,2)) :: sig_out
    real(DP), dimension(nRT) :: rand_test
    integer(i4b) :: i, q, seed
!    real(DP), dimension(size(solution,1), numQs) :: initial_pop

    n_x = size(solution,1)

! Initialize alpha for isres.
    alpha = 0.2d0

! Set up conv_check variables.
  if (iam == 0) then
    call init_conv_check(rms_target, numQs, x_lower, x_upper)
  end if
#ifdef HAVE_MPI
        call MPI_Bcast(NQ, 1, MPI_Integer, 0, MPI_Comm_World, ierr)
#endif

! Initialize the random seed.
    if (.not. seeded) then
        if (iam == 0) then
           seed = rand_seed
        end if
#ifdef HAVE_MPI
        Call MPI_Bcast(seed, 1, MPI_Integer4, 0, MPI_Comm_World, ierr)
        seed = mod(seed + iam, max_rnd_seed)
#endif
        Print*, "processor", iam, "SEED = ", seed
        Call init_genrand(seed)
        if (iam == 0) then
            call rnd_dist(rand_test)
            call check_rnd(seed, rand_test)
        end if
        seeded = .true.
    end if

! Allocate memory for sres variables.
    call alloc_sres()

! Do the optimization.
!    call isres(max_iters, x_out, initial_pop)
    if (present(IN_ALPHA)) then
        alpha = IN_ALPHA
        isres = .true.
    else
        isres = .false.
    end if
    if (present(INIT_POP)) then
        if (present(INIT_SIGMA)) then
            call sres(max_iters, max_time, x_out, sig_out, INIT_POP, INIT_SIGMA)
        else
            call sres(max_iters, max_time, x_out, sig_out, INIT_POP)
        end if
    else
        call sres(max_iters, max_time, x_out, sig_out)
    end if

! Return the solution.
    solution = x_out
    solution_sigma = sig_out

! Clean up sres variables.
    call dealloc_sres()

! Clean up conv_check variables.
  if (iam == 0) then
    call fin_conv_check()
  end if

    RETURN

END SUBROUTINE run_sres
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
SUBROUTINE sres(max_total_iter, time_limit, xout, sigma_out, initial_pop, &
                in_sigma) 
!  ======================================================================== 
! |  Stochastic Ranking using an Evolution Strategy algorithm              |
!  ========================================================================

  USE mt19937
  USE timer

  Implicit None
  
#ifdef HAVE_MPI
  include "mpif.h"
#endif

  Real(dp), Dimension(:,:), Intent(Out) :: xout
  Real(dp), Dimension(n_x,size(xout,2)), Intent(Out) :: sigma_out
  Integer(i4b), Intent(INOUT) :: max_total_iter
  Real(dp), Intent(IN) :: time_limit
  Real(dp), Dimension(:,:), Intent(IN), optional :: initial_pop
  Real(dp), Dimension(:,:), Intent(IN), optional :: in_sigma
!  Real(dp), Dimension(n,NQ), Intent(IN) :: x1_pop
!  Real(dp), Intent(IN) :: exitcrit          ! Exit criterion
!  Integer(i4b), Intent(IN) :: mu      ! Number of parents
!  Integer(i4b), Intent(IN) :: lambda  ! Number of individuals in the population
  Real(dp) :: fnout
  !Real(dp), Dimension(n_x,lambda) :: temp_x
  !Real(dp), Dimension(n_x,lambda,NQ) :: x
  double precision x(n_x,lambda,NQ)
  Real(dp), Dimension(n_x,lambda,NQ) :: sigma
  Real(dp), Dimension(n_x,lambda,NQ) :: xpr
  Real(dp), Dimension(n_x,lambda,NQ) :: sigmapr
  Real(dp), Dimension(lambda,NQ) :: f, phi
  Real(dp), Dimension(n_g) :: g
  Real(dp), Dimension(n_x,NQ) :: xbest
  Real(dp) :: chi, evout
  Real(dp), Dimension(NQ) :: fbest
  Integer(i4b), Dimension(NQ) :: Qrank
  Integer(i4b) :: miter
  Integer(i4b), Dimension(lambda,NQ) :: ind
  Integer(i4b) :: i, j, n, q, iters, filecnt, nfeval_sum
  Logical(lgt) :: feasible, stopping, writing_rst
  Logical(lgt) :: restartingQ3 = .false.
  integer, parameter :: iounit = 30
  character(len=10) :: timestr
  double precision :: temp
  integer :: pp, count, tempint, num_loops, starting_iter

  varphi = 1.0d0

! Initialize sigmaub.
  sigmaub = (x_upper - x_lower) / sqrt (dble(n_x))

  ! Now we are ready to initialize parameters
  ! note that we may want to use an initial search point
  if (isres) then
    chi = (1.0d0 / (2.0d0 * dble(n_x)) + 1.0d0 / (2.0d0 * sqrt(dble(n_x))))
    varphi = Sqrt((2d0/chi)*Log((1d0/alpha)*(Exp(varphi**2*chi/2.0d0) - &
                                                (1d0-alpha))))
  end if
  tau  = varphi / Sqrt(2d0*Sqrt(dble(n_x)))
  taud = varphi / Sqrt(2d0*Real(n_x))

  print*, n_x, lambda, NQ

! Initialize population.
!  Start with random population.
    !if (iam .eq. 0) then 
    do q = 1,NQ
      xpr(:,:,q)=0
      do j = 1+iam, lambda, n_procs
      !do j = 1, lambda
        do i = 1,n_x
          !print*, x_lower(i), x_upper(i)
          xpr(i,j,q) = x_lower(i)+(x_upper(i)-x_lower(i))*genrand_real1()
        end do
        sigmapr(:,j,q) = sigmaub
      end do
    end do
    !writE(*,'(14(g13.6,1x))'), xpr(:,:,2)
    !end if

    fbest = 10000000000000000.0d0
    xbest = 0.0d0
!  Overwrite the first population member and best values with starting dist.
!    if (present(initial_pop)) then
   
    !if (present(initial_pop) .and. iam == 0) then
    !  do q = 1, size(initial_pop,2)
    !    xpr(:,1,q) = initial_pop(:,q)
    !
    !    xbest(:,q) = initial_pop(:,q)
    !  ! Evaluate the best one.
    !    call model_objf(xbest(:,q), fbest(q), g)
    !    do i = 1, n_g
    !       if (g(i) > 0d0) then
    !           phi(1,q) = phi(1,q) + g(i)**2
    !       end if
    !    end do
    !
    !  end do
    !  if (present(in_sigma)) then
    !    do q = 1, size(initial_pop,2)
    !      sigmapr(:,1,q) = in_sigma(:,q)
    !    end do
    !  end if
    !end if

!! Overwrite the initial population if there is a starting policy file, this is
!! not a restart case, and damages are being considered (no BAU).
!  if (startpol .and. .not. restarting .and. damages) then
!    call get_startpol(xpr(:,:,1))
!  endif

  do q = 1,NQ
      ind(1:lambda,q) = (/ (i, i=1,lambda) /)
  end do

#ifdef HAVE_MPI
  Call Mpi_Allreduce(xpr, x, lambda*n_x*NQ, MPI_DOUBLE_PRECISION, &
                    MPI_SUM, MPI_COMM_WORLD, ierr)
  Call Mpi_Allreduce(sigmapr, sigma, lambda*n_x*NQ, Mpi_Double_Precision, &
                        Mpi_Sum, Mpi_Comm_World, ierr)
  !x = xpr
  !sigma = sigmapr
  !broadcast
  !Call Mpi_Bcast(x,     lambda*n_x*NQ, MPI_DOUBLE_PRECISION, 0, Mpi_Comm_World, ierr)
  !Call Mpi_Bcast(sigma, lambda*n_x*NQ, MPI_DOUBLE_PRECISION, 0, Mpi_Comm_World, ierr)
#else
  x = xpr
  sigma = sigmapr
#endif
    !writE(*,'(14(g13.6,1x))'), sigma(:,:,1)

!  xbest = 0.0d0
  feasible = .false.
!  new_fbest = .false.
  filecnt = 0
!  lastcheck = 0

! Start the main iteration loop
! Get ready.  Set miter to 10 unless the input file specifies more.
  if (max_total_iter > 0) then
!    miter = max(10, max_total_iter)
    miter = max_total_iter
  else
    miter = 10
  end if
  starting_iter = 0
! Get restart values and overwrite variables accordingly.
  if (restarting) then
!    call read_restart(n_x, lambda, starting_iter, filecnt, x_lower, x_upper, &
!                      sigma, ind, x, xbest, fbest, lastfbest, feasible)
  end if
  stopping = .false.
  writing_rst = .false.
  total_iters = 0
! Go.
  timeloop: do
  do iters = 1, miter
    model_iters = model_iters + 1
    total_iters = total_iters + 1

  ! Do the optimization step for each queue.
    do q = 1, NQ
      if (isres) then
        call inner_isres(mu, lambda, sigma(:,:,q), ind(:,q), &
                         x(:,:,q), f(:,q), phi(:,q))
      else
        call inner_sres(mu, lambda, sigma(:,:,q), ind(:,q), &
                         x(:,:,q), f(:,q), phi(:,q))
      end if
    !  print*, 'test1'
    ! Store the best feasible solution for run q.
      print*, ind(1,q)
      !print*, phi(ind(1,q))
      print*, f(ind(1,q),q)
      print*, size(fbest)
      print*, fbest(q)
        if ((phi(ind(1,q),q) == 0.0d0) .and. (f(ind(1,q),q)) < fbest(q)) then
            xbest(:,q) = x(:,ind(1,q),q)
            fbest(q) = f(ind(1,q),q)
!            evout = dice_rmse(xbest) ! <- CHECK THIS
!            new_fbest = .true.
            feasible = .true.
!            lastfbest = total_iters + starting_iter
!            if (iam == 0) then
!                call write_fbest(lastfbest, fbest(1), fbest(2), evout) 
!            end if
        end if
    end do
 
    ! Write out xbest from runs A and B, but not until the flurry of updates
    ! is over.
!    if (new_fbest .and. (starting_iter+total_iters-lastfbest > 99) ) then
!         ((iters-lastfbest > 99).or.(iters-lastcheck > 499)))  then
    print*, 'test2'
    if ( (refresh > 0) .and. (mod(total_iters,refresh)==0) ) then
#ifdef HAVE_MPI
          call MPI_Reduce(nfeval, nfeval_sum, 1, MPI_Integer, MPI_SUM, 0, &
                            MPI_COMM_WORLD, ierr)
          nfe_total = nfe_total + nfeval_sum
#else
          nfe_total = nfe_total + nfeval
#endif
        nfeval = 0
!        reset_nfe = .false.
        print*, 'test3'
        if (iam == 0) then
            call conv_check(total_iters, refresh, fbest, xbest, Qrank, &
                            stopping, restartingQ3) !, reset_nfe )
        end if
#ifdef HAVE_MPI
          call MPI_Bcast(Qrank, NQ, MPI_Integer, 0, MPI_Comm_World, ierr)
          call MPI_Bcast(stopping, 1, MPI_Logical, 0, MPI_Comm_World, ierr)
          call MPI_Bcast(restartingQ3, 1, MPI_Logical, 0, MPI_Comm_World, ierr)
!          call MPI_Bcast(reset_nfe, 1, MPI_Logical, 0, MPI_Comm_World, ierr)
#endif
!        if (reset_nfe) nfeval = 0
        print*, 'test4'
        if (restartingQ3) then
!          print *, "Restarting Q", Qrank(NQ)
          do j = 1+iam, lambda, n_procs
            do i = 1,n_x
              xpr(i,j,Qrank(NQ)) = x_lower(i)+(x_upper(i)-x_lower(i))*genrand_real1()
            end do
            sigmapr(:,j,Qrank(NQ)) = sigmaub
          end do
          fbest(Qrank(NQ)) = 100000000.0d0
          xbest(:,Qrank(NQ)) = 0.0d0
#ifdef HAVE_MPI
          Call Mpi_Allreduce(xpr(:,:,Qrank(NQ)), x(:,:,Qrank(NQ)), &
            lambda*n_x, Mpi_Double_Precision, Mpi_Sum, Mpi_Comm_World, ierr)
          Call Mpi_Allreduce(sigmapr(:,:,Qrank(NQ)), sigma(:,:,Qrank(NQ)), &
                lambda*n_x, Mpi_Double_Precision, Mpi_Sum, Mpi_Comm_World, ierr)
#else
          x(:,:,Qrank(NQ)) = xpr(:,:,Qrank(NQ))
          sigma(:,:,Qrank(NQ)) = sigmapr(:,:,Qrank(NQ))
#endif
        end if
       
       filecnt = filecnt + 1
    end if
    if (stopping) exit
  end do
  
  if ( iam == 0 ) then
  ! Check the time and signal to write restart files and/or stop.
    call time_check(total_iters, time_limit, max_total_iter, miter, &
                    writing_rst, stopping)
    if (writing_rst) then
!        call write_restart(n_x, lambda, total_iters+starting_iter, filecnt, &
!                           x_lower, x_upper, sigma, ind, x, xbest, fbest, &
!                           lastfbest, feasible)
        writing_rst = .false.
    end if
  end if
#ifdef HAVE_MPI
  Call Mpi_Bcast(stopping, 1, Mpi_Logical, 0, Mpi_Comm_World, ierr)
  Call Mpi_Bcast(miter, 1, Mpi_Integer4, 0, Mpi_Comm_World, ierr)
#endif
  if (stopping) exit
  end do timeloop       ! End of the main loop.

! Report the number of iterations completed.
    max_total_iter = total_iters

! In case conv_check has not yet been called, find Qbest.
    if (iam == 0) then
        Qrank = sort_val(fbest)
        do q = 1, size(xout,2)
            xout(:,q) = xbest(:,Qrank(q))
            sigma_out(:,q) = sigma(:,ind(1,q),Qrank(q))
        end do
        fnout = fbest(Qrank(1))

! Report the answer or lack thereof.
      if (feasible) then
         Print*, "feasible solution: F(X) = ", fnout
         Print*, "X = ", xout(:,1)
      else
         Print*, "No feasible solution found, increase Miter or decrease pf"
      end if
   end if

   RETURN

END SUBROUTINE sres
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
SUBROUTINE inner_sres(mu, lambda, sigma, ind, x, f, phi)
!===========================================================================
! This subroutine handles the inner operations of the main iteration loop
! of the SRES scheme.  It is rewritten here as a subroutine to facilitate
! running two simultaneous optimization queues for automatic detection of
! convergence.  BCT 4 Apr 2005 <btuttle@psu.edu>
!===========================================================================

  USE mt19937

#ifdef HAVE_MPI
  include "mpif.h"
#endif

  Integer(i4b), Intent(IN) :: mu      ! Number of parents
  Integer(i4b), Intent(IN) :: lambda  ! Number of individuals in the population
  Real(dp), Dimension(n_x,lambda), intent(INOUT) :: sigma
  Integer(i4b), Dimension(lambda), intent(INOUT) :: ind
  Real(dp), Dimension(n_x,lambda), intent(INOUT) :: x
  Real(dp), Dimension(lambda), intent(OUT) :: f, phi

! Internal variables: g, fpr, phipr, sigmapr, xpr, i, j, k, normal, m, cnt
  Real(dp), Dimension(lambda) :: fpr, phipr
  Real(dp), Dimension(n_x,lambda) :: xpr, sigmapr
  Real(dp), Dimension(n_g) :: g
  Integer(i4b) :: i, j, k, m, cnt
  Real(dp) :: normal

    ! Zero some variables, because we will use Mpi_Sum in Mpi_Allreduce
    fpr = 0d0
    phipr = 0d0
    sigmapr = 0d0
    xpr = 0d0

    ! Generate new search points
    Do j = iam+1, lambda, n_procs
      i = Mod(j-1,mu) + 1
      normal = nrand()
      Do k = 1, n_x
        ! strategy space: arithmetic crossover and lognormal distibution
        m = Mod(genrand_int31(),mu) + 1
        sigmapr(k,j) = (sigma(k,ind(i)) + sigma(k,ind(m))) / 2.0d0
        sigmapr(k,j) = sigmapr(k,j) * Exp(taud * normal + tau * nrand())
        If (sigmapr(k,j) > sigmaub(k)) sigmapr(k,j) = sigmaub(k)
        ! object space: Gaussian search distribution (rotation sensitive)
        cnt = 0;
        Do  
          xpr(k,j) = x(k,ind(i)) + sigmapr(k,j) * nrand()
          cnt = cnt + 1
          
         If ((xpr(k,j)>=x_lower(k)) .and. (xpr(k,j)<=x_upper(k))) Then
            exit
          End If
          If (cnt > 10) Then
             xpr(k,j) = x(k,ind(i))
             exit
          End If
        End Do

      End Do

! Call the objective function.
      call model_objf(xpr(:,j), fpr(j), g)
      !print*, j, fpr(j), xpr(1,j)
      do k = 1, n_g
        if (g(k) > 0d0) then
          phipr(j) = phipr(j) + g(k)**2
        end if
      end do

    End Do

    ! Synchonize with all other nodes. 

#ifdef HAVE_MPI
! Combine fpr from all processes to f on the root process.  Lambda elements.
    Call Mpi_Reduce(fpr, f, lambda, Mpi_Double_Precision, Mpi_Sum, 0, &
                    Mpi_Comm_World, ierr)

! Combine phipr from all processes to phi on the root process.  Lambda elements.
    Call Mpi_Reduce(phipr, phi, lambda, Mpi_Double_Precision, Mpi_Sum, 0, &
                    Mpi_Comm_World, ierr)

! Combine xpr from all processes and distribute to x on all processes.  
! Lambda*n is the number of elements.
    Call Mpi_Allreduce(xpr, x, lambda*n_x, Mpi_Double_Precision, Mpi_Sum, &
                    Mpi_Comm_World, ierr)

! Combine sigmapr from all processes and distribute to sigma on all processes.  
! Lambda*n is the number of elements.
    Call Mpi_Allreduce(sigmapr, sigma, lambda*n_x, Mpi_Double_Precision, &
                    Mpi_Sum, Mpi_Comm_World, ierr)
#else
   f = fpr
   phi = phipr
   x = xpr
   sigma = sigmapr
#endif
    ! Only the first node does the sorting and the Bcasts the ordering
   If (iam == 0) Then
      Call srsort(f,phi,ind)

   End If

#ifdef HAVE_MPI
    Call Mpi_Bcast(ind, lambda, Mpi_Integer4, 0, Mpi_Comm_World, ierr)
    Call Mpi_Bcast(f, lambda, Mpi_Double_Precision, 0, Mpi_Comm_World, ierr)
    Call Mpi_Bcast(phi, lambda, Mpi_Double_Precision, 0, Mpi_Comm_World, ierr)
#endif

    RETURN

END SUBROUTINE inner_sres
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
SUBROUTINE inner_isres(mu, lambda, sigma, ind, x, f, phi)
!===========================================================================
! This subroutine handles the inner operations of the main iteration loop
! of the iSRES scheme.  It is rewritten here as a subroutine to facilitate
! running two simultaneous optimization queues for automatic detection of
! convergence.  BCT 4 Apr 2005
!===========================================================================

  USE mt19937
!  USE mod_sow

#ifdef HAVE_MPI
  include "mpif.h"
#endif

  Integer(i4b), Intent(IN) :: mu      ! Number of parents
  Integer(i4b), Intent(IN) :: lambda  ! Number of individuals in the population
  Real(dp), Dimension(n_x,lambda), intent(INOUT) :: sigma
  Integer(i4b), Dimension(lambda), intent(INOUT) :: ind
  Real(dp), Dimension(n_x,lambda), intent(INOUT) :: x
  Real(dp), Dimension(lambda), intent(OUT) :: f, phi

! Internal variables: g, fpr, phipr, sigmapr, xpr, i, j, k, normal, m, cnt
  Real(dp), Dimension(lambda) :: fpr, phipr
  Real(dp), Dimension(n_x,lambda) :: xpr, sigmapr
  Real(dp), Dimension(n_g) :: g
  Integer(i4b) :: i, j, k, m, cnt
  Real(dp) :: normal
!  Integer(i4b) :: stat(Mpi_Status_Size)
  Real(dp) :: ftemp
!  Integer(i4b) :: sow
  integer pp, qq, jj, kk, tempint, num_loops
  double precision temp

    ! Zero some variables, because we will use Mpi_Sum in Mpi_Allreduce
    fpr = 0d0
    phipr = 0d0
    sigmapr = 0d0
    xpr = 0d0

    ! Generate new search points
    do j = iam + 1, lambda, n_procs
      i = Mod(j-1,mu) + 1
      normal = nrand()
      do k = 1, n_x
        ! object space: differential variation
        if (j < mu) then
          sigmapr(k,j) = sigma(k,ind(i))
          xpr(k,j) = x(k,ind(i)) + gamma * (x(k,ind(1)) - x(k,ind(i+1)))
        Else
        ! strategy space: arithmetic crossover and lognormal distibution
          sigmapr(k,j) = sigma(k,ind(i)) * Exp(taud * normal + tau * nrand())
          if (sigmapr(k,j) > sigmaub(k)) sigmapr(k,j) = sigmaub(k)
        ! object space: Gaussian search distribution
          xpr(k,j) = x(k,ind(i)) + sigmapr(k,j) * nrand()
        end if
        cnt = 0;
        do 
          if ((xpr(k,j)>=x_lower(k)) .and. (xpr(k,j)<=x_upper(k))) then
            exit
          end if
          cnt = cnt + 1
          if (cnt > 10) then
             xpr(k,j) = x(k,ind(i))
             exit
          end if
          xpr(k,j) = x(k,ind(i)) + sigmapr(k,j) * nrand()
        end do
        if (j >= mu) then
          sigmapr(k,j) = sigma(k,ind(i))+alpha*(sigmapr(k,j)-sigma(k,ind(i)))
        end if

      end do

! Call the objective function.
      call model_objf(xpr(:,j), fpr(j), g)
      do k = 1, n_g
        if (g(k) > 0d0) then
           phipr(j) = phipr(j) + g(k)**2
        end if
      end do
    end do

#ifdef HAVE_MPI
  ! Synchonize with all other nodes
    Call Mpi_Reduce(fpr, f, lambda, MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
                    MPI_COMM_WORLD, ierr)
    Call Mpi_Reduce(phipr, phi, lambda, MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
                    MPI_COMM_WORLD, ierr)
    Call Mpi_Allreduce(xpr, x, lambda*n_x, MPI_DOUBLE_PRECISION, MPI_SUM, &
                        MPI_COMM_WORLD, ierr)
    Call Mpi_Allreduce(sigmapr, sigma, lambda*n_x, MPI_DOUBLE_PRECISION, &
                        MPI_SUM, MPI_COMM_WORLD, ierr)
#else
   f = fpr
   phi = phipr
   x = xpr
   sigma = sigmapr
#endif

    ! Only the first node does the sorting and the Bcasts the ordering
    if (iam == 0) then
      Call srsort(f,phi,ind)

    end if

#ifdef HAVE_MPI
    Call Mpi_Bcast(ind, lambda, Mpi_Integer4, 0, Mpi_Comm_World, ierr)
    Call Mpi_Bcast(f, lambda, Mpi_Double_Precision, 0, Mpi_Comm_World, ierr)
    Call Mpi_Bcast(phi, lambda, Mpi_Double_Precision, 0, Mpi_Comm_World, ierr)
#endif

    RETURN

END SUBROUTINE inner_isres
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
FUNCTION nrand()
  USE mt19937
  Implicit None
  Real(dp) :: nrand, v, u
  Do
    v = genrand_real1()
    Do
      u = genrand_real1()
      If (u /= 0) Exit
    End Do
    ! the constant 1.715 ... = sqrt(8/e)
    nrand = 1.71552776992141359295 * (v-0.5)/u
    If (nrand**2 <= -4.0*Log(u)) Exit
  End Do
END FUNCTION nrand
!---------------------------------------------------------------------------
  
!---------------------------------------------------------------------------
SUBROUTINE srsort(f,phi,ind)
  USE mt19937
  Implicit None
  
  Real(dp), Intent(In), Dimension(:)  :: phi, f
  Integer(i4b), Intent(Out), Dimension(:) :: ind
  integer counter
  double precision  temp
  ! input parameters f,phi - vectors to be sorted
  !                    ind - vector of length >= Size(f)
  ! f,phi are not altered by this routine

  ! output parameter - ind - sequence of indices 1,...,N permuted in the same
  !                          fashion as x would be. Thus, the ordering on 
  !                          x is defined by y(i) = x(ind(i))
  
  Integer(i4b) :: i, j, k
  ! local parameters -
  ! i,j : simple counters
  ! k   : used as a stopping criteria for complete sweep
  
  ! initialize
  Do i = 1, Size(f)
    ind(i) = i
  End Do
  
  ! Perform stochastic bubble sort
  Do i = 1, Size(f)
    k = 0
    Do j = 1, (Size(f) - 1)

      If (((phi(ind(j)) == phi(ind(j+1))) .and. (phi(ind(j)) == 0)) .or. (genrand_real1() < pf)) Then
        If (f(ind(j+1)) < f(ind(j))) Then
          k = ind(j)
          ind(j) = ind(j+1)
          ind(j+1) = k 
        End If
      Else
        If (phi(ind(j+1)) < phi(ind(j))) Then
          k = ind(j)
          ind(j) = ind(j+1)
          ind(j+1) = k 
        End If
      End If
    End Do
    If (k == 0) Then
      Exit
    End If
  End Do
  
END SUBROUTINE srsort
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
FUNCTION randate() RESULT(rnd)
!===========================================================================
!   Generate a pseudo-random number out of the date_and_time() subroutine.
!   Brian Tuttle 25Feb05 <btuttle@psu.edu>
!===========================================================================

    integer(i4b) :: rndint, maxrnd, rnd
    integer(i4b) :: i, dig1, dig2, dig3, num1, num2, num3
    real(DP) :: nrnd
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
    rnd = int(100000000 * nrnd)

END FUNCTION randate
!---------------------------------------------------------------------------

!!---------------------------------------------------------------------------
!FUNCTION getseed()
!
!  Integer(i4b) :: getseed
!
!  ! Use the random device in Linux (if this fails we need to use the wall clock or user input)
!  ! The alternative Linux device is /dev/random 
!  Open (Unit=76, File='/dev/urandom', Form='Unformatted',Access='Direct', Status='Old', Recl=1)
!  Read (76, Rec=1) getseed
!  Close(76)
!
!END FUNCTION getseed
!!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
SUBROUTINE alloc_sres()
!  =========================================================================
! |  This subroutine allocates memory for the variables in the ceo module.  |
!  =========================================================================

    integer :: astat

    astat = 0
    allocate(sigmaub(n_x), STAT=astat)
    if (astat > 0) call error_msg(1,"sigmaub")

    RETURN
    
END SUBROUTINE alloc_sres
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
SUBROUTINE dealloc_sres()
!  ===========================================================================
! |  This subroutine deallocates memory for the variables in the ceo module.  |
!  ===========================================================================

    integer :: astat

    astat = 0
    deallocate(sigmaub, STAT=astat)
    if (astat > 0) call error_msg(2,"sigmaub")

    RETURN
    
END SUBROUTINE dealloc_sres
!---------------------------------------------------------------------------
END MODULE ceo
