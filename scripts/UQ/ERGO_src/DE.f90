MODULE DE

  USE global
  USE model
  USE converge
  USE rndseed

  implicit none

#ifdef HAVE_MPI
  include "mpif.h"
#endif

  integer, parameter :: refresh=1
  integer, parameter :: interval = 99
  integer, parameter :: iwrite=7
#ifdef ERGO
  integer :: strategy 
  integer, dimension(3) :: method
#else
  integer, parameter :: strategy=6
  integer, dimension(3), parameter :: method=(/0, 1, 0/)
#endif
  double precision, parameter :: VTR=-100000000d0, CR_XC=0.5d0
  integer :: NP
  integer :: n_x
!  double precision, dimension(:,:), allocatable :: bestmem_XC
!  double precision :: prob_thresh = 0.0d0
!  double precision :: prob_coll = 0.0d0
#ifdef ERGO
! Put these here publicly because input_vars is not included as a module in DE.
  integer, public :: inDE_pop       ! factor to multiply n_x to get NP
  integer, public :: inDE_rndscale  ! method(1)
  integer, public :: inDE_strat     ! strategy 1 to 6
  double precision, public :: inDE_F_CR ! Crossover factor: 0 implies random
                                ! and method(2) = 1.  Otherwise method(2) = 0
#endif

CONTAINS
!---------------------------------------------------------------------------
SUBROUTINE run_DE(solution, max_iters, max_time, rms_target, numQs)
!  =====================================================================
! |
!  =====================================================================

    USE mt19937

    implicit none

    double precision, dimension(:,:), intent(OUT) :: solution
    integer, intent(INOUT) :: max_iters
    double precision, intent(IN) :: max_time
    double precision, intent(IN) :: rms_target
    integer, intent(IN) :: numQs
    integer :: seed
    double precision, dimension(size(solution,1),size(solution,2)) :: bestmember
    double precision, dimension(nRT) :: rand_test
    

    n_x = size(solution,1)
#ifdef ERGO
    NP = inDE_pop * n_x
    method(1) = inDE_rndscale
    method(2) = 0
    if (inDE_F_CR == 0.0d0) method(2) = 1
    strategy = inDE_strat
#else
    NP = 10*n_x
#endif

! Set up conv_check variables.
  if (iam == 0) then
    call init_conv_check(rms_target, numQs, x_lower, x_upper)
  end if
#ifdef HAVE_MPI
        call MPI_Bcast(NQ, 1, MPI_Integer, 0, MPI_Comm_World, ierr)
#endif

! Initialize random seed.
    if (.not. seeded) then
        if (iam == 0) then
          seed = rand_seed
        end if
#ifdef HAVE_MPI
        call MPI_Bcast(seed, 1, MPI_Integer4, 0, MPI_Comm_World, ierr)
        seed = mod(seed + iam, max_rnd_seed)
#endif
        print *, "processor", iam, "seed = ", seed
        call init_genrand(rand_seed)
        if (iam == 0) then
            call rnd_dist(rand_test)
            call check_rnd(seed, rand_test)
        end if
        seeded = .true.
    end if

! Do the optimization.
    call DE_Fortran90(max_iters, max_time, bestmember)

! Return the best solution.
    solution = bestmember

! Clean up conv_check variables.
  if (iam == 0) then
    call fin_conv_check()
  end if

    RETURN

END SUBROUTINE run_DE
!---------------------------------------------------------------------------

!!---------------------------------------------------------------------------
!SUBROUTINE Rosen(inv_choice, max_iterations)
!!  =====================================================================
!! |
!!  =====================================================================
!
!    implicit none
!
!    character(len=3), intent(IN) :: inv_choice
!    integer, intent(IN) :: max_iterations
!    double precision, allocatable :: x(:)
!    double precision :: f
!    integer :: i
!    integer , dimension(8) :: time 
!    integer :: alloc_stat
!    integer :: seed(2)
!    character(len=6) :: fname
!    intrinsic date_and_time
!
!    itermax = max_iterations
!
!    if(inv_choice .eq. 'BAU') then
!        fname = 'BAUout'
!    else
!        fname = 'OPTout'
!    end if
!    open(iwrite,file=fname)
!    call date_and_time(values=time)
!    write(unit=iwrite, FMT=11) time(1:3), time(5:7)
!
!!    call setup_go_problem
!    call init_dice(inv_choice)
!
!    allocate(x(n_x))
!    if (inv_choice .eq. 'BAU') then
!        alloc_stat = 0
!        allocate(bau_inv(nop), STAT=alloc_stat)
!        if (alloc_stat > 0) call error_msg(1,"bau_inv")
!    end if
!
!       seed(1) = randate()  
!!       seed(1) = 1 
!       seed(2) = seed(1)**0.5  
!       print*,"seed =",seed
!
!call random_seed(put=seed)
!
!    alloc_stat = 0
!    allocate(bestmem_XC(n_x),stat=alloc_stat)
!    if (alloc_stat > 0) then
!       print *, "Error allocating memory for bestmem_XC"
!       stop
!    end if
!
!print*,"NP =",NP
!print*,"itermax =",itermax
!print*,"method(1:3) =",method(1:3)
!write(*,201) F_XC, CR_XC, F_CR
!
!    call DE_Fortran90
!              
!    write(iwrite,205) NP, nfeval, method(1:3)
!    write(iwrite,FMT=201) F_XC, CR_XC, F_CR
!    write(iwrite,FMT=200) bestval
!    do i=1,n_x
!       write(iwrite,FMT=202) i,bestmem_XC(i)
!    end do
!200 format(/2x, 'Bestval=', ES14.7)
!201 format(2x, 'F_XC =',F6.3, 2x, 'CR_XC =', F6.3, 2x, 'F_CR =', F6.3)
!!202 format(2x, 'best_XC(',I3,') =',ES14.7)
!202 format(2x, I3, ES16.7)
!205 format(2x, 'NP=', I4, 4x, 'No. function call =', I9, &
!         /2x, 'method(1:3) =',3I3)
!    call date_and_time(values=time)
!    write(unit=iwrite, FMT=10)time(1:3), time(5:7)
!10  format(/1x, 'End of Program. Date:', I4, '/', I2,'/', I2, ', Time: ', I2,':',I2,':',I2) 
!11  format(1x, 'Beginning of Program. Date:', I4, '/', I2,'/', I2, ', Time: ', I2,':',I2,':',I2)
!    
!    if (inv_choice .eq. 'BAU') then
!        bau_inv = bestmem_XC(nop+1:n_x)
!    end if
!    deallocate(x)
!    deallocate(bestmem_XC)
!    call dealloc_dvar
!
!END SUBROUTINE Rosen
  
!---------------------------------------------------------------------------
SUBROUTINE DE_Fortran90(max_total_iter, time_limit, bestmem)
!=======================================================================
!  
!   Differential Evolution for Optimal Control Problems
!
!-----------------------------------------------------------------------
!    This Fortran 90 program translates from the original MATLAB 
!    version of differential evolution (DE). This FORTRAN 90 code 
!    has been tested on Compaq Visual Fortran v6.1. 
!    Any users new to the DE are encouraged to read the article of 
!    Storn and Price. 
!
!  Refences:
!    Storn, R., and Price, K.V., (1996). Minimizing the double precision 
!    function of the ICEC96 contest by differential evolution. IEEE conf. 
!    on Evolutionary Comutation, 842-844.
!
!  This Fortran 90 program written by Dr. Feng-Sheng Wang 
!  Department of Chemical Engineering, National Chung Cheng University, 
!  Chia-Yi 621, Taiwan, e-mail: chmfsw@ccunix.ccu.edu.tw
!=======================================================================
!             obj : The user provided file for evlauting the objective function.
!                   subroutine obj(xc,fitness)
!                   where "xc" is the double precision decision parameter 
!                   vector.(input)
!                   "fitness" is the fitness value.(output)
!             n_x : Dimension of the double precision decision parameters.
!    x_lower(n_x) : The lower bound of the double precision decision parameters.
!    x_upper(n_x) : The upper bound of the double precision decision parameters.
!             VTR : The expected fitness value to reach.
!              NP : Population size.
!         itermax : The maximum number of iteration.
!            F_XC : Mutation scaling factor for double precision decision 
!                   parameters.
!           CR_XC : Crossover factor for double precision decision parameters.
!        strategy : The strategy of the mutation operations is used in HDE.
!         refresh : The intermediate output will be produced after "refresh"
!                   iterations. No intermediate output will be produced if
!                   "refresh < 1".
!        interval : The intermediate output will be produced after "interval"
!                   iterations following the last update to bestval.
!          iwrite : The unit specfier for writing to an external data file.
! bestmem_XC(n_x) : The best double precision decision parameters.
!         bestval : The best objective function.
!          nfeval : The number of function call.
!         method(1) = 0, Fixed mutation scaling factors (F_XC)
!                   = 1, Random mutation scaling factors F_XC=[0, 1]
!                   = 2, Random mutation scaling factors F_XC=[-1, 1] 
!         method(2) = 1, Random combined factor (F_CR) used for strategy = 6
!                        in the mutation operation 
!                   = other, fixed combined factor provided by the user 
!         method(3) = 1, Saving results in a data file.
!                   = other, displaying results only.
! ========================================================================
    
    USE timer

    implicit none

    integer, intent(INOUT) :: max_total_iter
    double precision, intent(IN) :: time_limit
    double precision, dimension(:,:), intent(OUT) :: bestmem
    double precision, dimension(NP,n_x,NQ) :: pop_XC, bm_XC, mui_XC, &
         mpo_XC, popold_XC, rand_XC, ui_XC
    double precision, dimension(NQ) :: F_XC, F_CR
    integer, dimension(NQ) :: Qrank
    integer, dimension(NP) :: rot, rt
    integer, dimension(NP,NQ) :: a1, a2, a3, a4, a5
    integer, dimension(4) :: ind
    double precision, dimension(NP) :: tempval, tempval_pr, phi, phi_pr
    double precision, dimension(NP,NQ) :: val
    double precision, dimension(NQ) :: bestval
!    double precision, dimension(n_x,NQ) :: bestmemit_XC
    double precision, dimension(n_x,NQ) :: bestmem_XC
    intrinsic max, min, random_number, mod, abs, any, all, maxloc
    double precision :: fbest, phibest
!    double precision, dimension(n_g) :: gbest
!    logical :: feasible = .true.
    logical :: done !, reset_nfe
    logical :: writing_rst = .false. , restartingQ3 = .false.
    integer :: i, ibest, total_iters, miter, iter, q, nfeval_sum

! Initialize mutation factors.
    F_XC = 0.8d0
#ifdef ERGO
    F_CR = inDE_F_CR
#else
    F_CR = 0.8d0
#endif
  
! Initialize the population and evaluate fitness functions to find the best 
! member.
    pop_XC=0.0d0
    val=0.0d0
    bestval = 0.0d0
    do q=1,NQ
        call init_pop(pop_XC(:,:,q), bestmem_XC(:,q), val(:,q), bestval(q))
    end do
 
    if (iam == 0) then
        bm_XC=0.0d0
        rot=(/(i,i=0,NP-1)/)
    end if

    iter=1 

! Prepare for main iteration loop.
    if (max_total_iter > 0) then
!        miter = max(10, max_total_iter)
        miter = max_total_iter
    else
        miter = 10
    end if
    total_iters = 0

! Perform evolutionary computation: Start TIME LOOP.

 timeloop: do
  
  do iter = 1, miter

    model_iters = model_iters + 1
    total_iters = total_iters + 1

    if (iam == 0) then      ! MASTER processor task
       popold_XC=pop_XC

! Mutation operation:
       do q=1,NQ
          ind=randperm(4)
          a1(:,q)=randperm(NP)
          rt=mod(rot+ind(1),NP)
          a2(:,q)=a1(rt+1,q)
          rt=mod(rot+ind(2),NP)
          a3(:,q)=a2(rt+1,q)
          rt=mod(rot+ind(3),NP)
          a4(:,q)=a3(rt+1,q)
          rt=mod(rot+ind(4),NP)
          a5(:,q)=a4(rt+1,q)
       end do
!       bm_XC=spread(bestmemit_XC, DIM=1, NCOPIES=NP)
       bm_XC=spread(bestmem_XC, DIM=1, NCOPIES=NP)

! Generating a random scaling factor
       select case (method(1))
       case (1)
          call rnd_dist(F_XC)
       case(2)
          call rnd_dist(F_XC)
          F_XC=2.0d0*F_XC-1.0d0
       end select

! Select a mutation strategy.
       select case (strategy)
       case (1)
          do q=1,NQ
             ui_XC(:,:,q) = bm_XC(:,:,q) + &
                       F_XC(q)*(popold_XC(a1(:,q),:,q)-popold_XC(a2(:,q),:,q))
          end do
          
       case default
          do q=1,NQ
             ui_XC(:,:,q) = popold_XC(a3(:,q),:,q) + &
                       F_XC(q)*(popold_XC(a1(:,q),:,q)-popold_XC(a2(:,q),:,q))
          end do
          
       case (3)
          do q=1,NQ
             ui_XC(:,:,q) = popold_XC(:,:,q) + &
                       F_XC(q)*(bm_XC(:,:,q)-popold_XC(:,:,q) + &
                                popold_XC(a1(:,q),:,q)-popold_XC(a2(:,q),:,q))
          end do
          
       case (4)
          do q=1,NQ
             ui_XC(:,:,q) = bm_XC(:,:,q) + &
                     F_XC(q)*(popold_XC(a1(:,q),:,q)-popold_XC(a2(:,q),:,q) + &
                              popold_XC(a3(:,q),:,q)-popold_XC(a4(:,q),:,q))
          end do
          
       case (5)
          do q=1,NQ
             ui_XC(:,:,q) = popold_XC(a5(:,q),:,q) + &
                     F_XC(q)*(popold_XC(a1(:,q),:,q)-popold_XC(a2(:,q),:,q) + &
                              popold_XC(a3(:,q),:,q)-popold_XC(a4(:,q),:,q))
          end do
       case (6) ! A linear crossover combination of bm_XC and popold_XC
          if (method(2) == 1) call rnd_dist(F_CR) 
          do q=1,NQ
             ui_XC(:,:,q) = popold_XC(:,:,q) + &
                       F_CR(q)*(bm_XC(:,:,q)-popold_XC(:,:,q)) + &
                       F_XC(q)*(popold_XC(a1(:,q),:,q)-popold_XC(a2(:,q),:,q))
          end do
          
       end select

! Crossover operation:
       call rnd_dist(rand_XC)
       mui_XC=0.0d0
       mpo_XC=0.0d0
       where (rand_XC < CR_XC)
          mui_XC=1.0d0
       elsewhere
          mpo_XC=1.0d0
       end where

       ui_XC=popold_XC*mpo_XC+ui_XC*mui_XC

! Evaluate fitness functions and find the best member.
       do q=1,NQ
         do i=1,NP
! Confine each of feasible individuals in the lower-upper bound.
           ui_XC(i,:,q)=max(min(ui_XC(i,:,q),x_upper),x_lower)
         end do
       end do

    end if      !! End MASTER task

#ifdef HAVE_MPI
    call MPI_Bcast(ui_XC, NP*n_x*NQ, MPI_Double_Precision, 0, &
                    MPI_Comm_World, ierr)
#endif

       do q=1,NQ

         tempval_pr = 0.0d0
         phi_pr = 0.0d0
         do i=1+iam,NP,n_procs
           call FTN(ui_XC(i,:,q), tempval_pr(i), phi_pr(i))
         end do

#ifdef HAVE_MPI
         Call MPI_Reduce(tempval_pr, tempval, NP, MPI_Double_Precision, &
                        MPI_Sum, 0, MPI_Comm_World, ierr)
         Call MPI_Reduce(phi_pr, phi, NP, MPI_Double_Precision, &
                        MPI_Sum, 0, MPI_Comm_World, ierr)
#else
         tempval = tempval_pr
         phi = phi_pr
#endif

         if (iam == 0) then     !** MASTER task

           do i=1,NP
             if (tempval(i) < val(i,q) .and. phi(i) == 0.0d0) then
               pop_XC(i,:,q)=ui_XC(i,:,q)
               val(i,q) = tempval(i)

               if (tempval(i) < bestval(q)) then
                 bestval(q) = tempval(i)
                 bestmem_XC(:,q) = ui_XC(i,:,q)
               end if

             end if
           end do

         end if

       end do                   ! end MASTER task **
            
!       bestmemit_XC = bestmem_XC
       done = .false.

       if( (refresh > 0) .and. (mod(total_iters,refresh)==0)) then
#ifdef HAVE_MPI
          call MPI_Reduce(nfeval, nfeval_sum, 1, MPI_Integer, MPI_SUM, 0, &
                            MPI_COMM_WORLD, ierr)
          nfe_total = nfe_total + nfeval_sum
#else
          nfe_total = nfe_total + nfeval
#endif
          nfeval = 0
!          reset_nfe = .false.

          if (iam == 0) then
              call conv_check(total_iters, refresh, bestval, bestmem_XC, &
                              Qrank, done, restartingQ3) !, reset_nfe)
          end if
#ifdef HAVE_MPI
          call MPI_Bcast(Qrank, NQ, MPI_Integer, 0, MPI_Comm_World, ierr)
          call MPI_Bcast(done, 1, MPI_Logical, 0, MPI_Comm_World, ierr)
          call MPI_Bcast(restartingQ3, 1, MPI_Logical, 0, MPI_Comm_World, ierr)
!          call MPI_Bcast(reset_nfe, 1, MPI_Logical, 0, MPI_Comm_World, ierr)
#endif
!          if (reset_nfe) nfeval = 0

          if (restartingQ3) then
!              print *, "** Restarting Q:", Qrank(NQ)
              call init_pop(pop_XC(:,:,Qrank(NQ)), bestmem_XC(:,Qrank(NQ)), &
                            val(:,Qrank(NQ)), bestval(Qrank(NQ)))
              restartingQ3 = .false.
          end if

       end if
            
!       if ( minval(bestval) <= VTR .and. refresh > 0) then
!          write(unit=iwrite, FMT=*) ' The best fitness is smaller than VTR'
!          write(unit=*, FMT=*) 'The best fitness is smaller than VTR' 
!          exit
!       endif
       if (done) exit

  end do      ! End miter loop.

  if (iam  == 0) then
    call time_check(total_iters, time_limit, max_total_iter, miter, &
                  writing_rst, done)
  end if
#ifdef HAVE_MPI
  call MPI_Bcast(miter, 1, MPI_Integer, 0, MPI_Comm_World, ierr)
  call MPI_Bcast(writing_rst, 1, MPI_Logical, 0, MPI_Comm_World, ierr)
  call MPI_Bcast(done, 1, MPI_Logical, 0, MPI_Comm_World, ierr)
#endif

  if (done) exit

 end do timeloop     
! End main timeloop. ***

! Report the number of iterations completed.
    max_total_iter = total_iters

    if (iam == 0) then
    ! In case conv_check has not been called yet, find Qbest.
        Qrank = sort_val(bestval) 
    end if
#ifdef HAVE_MPI
    call MPI_Bcast(Qrank, NQ, MPI_Integer, 0, MPI_Comm_World, ierr)
    call MPI_Bcast(bestmem_XC, n_x*NQ, MPI_Double_Precision, 0, &
                   MPI_Comm_World, ierr)
#endif
    call FTN(bestmem_XC(:,Qrank(1)), fbest, phibest)

    if (iam == 0 .and. phibest > 0.0d0 ) then
       print *, "Solution is not feasible."
       print *, "Increase pen_val in subroutine FTN in DE.f90" 
    end if

    do q=1,size(bestmem,2)
        bestmem(:,q) = bestmem_XC(:,Qrank(q))
    end do

    RETURN

END SUBROUTINE DE_Fortran90
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
SUBROUTINE init_pop(poparray, best_mem, obfval, best_obfval)
!  =======================================================================
! |
!  =======================================================================

    implicit none

    double precision, dimension(NP,n_x), intent(OUT) :: poparray
    double precision, dimension(n_x), intent(OUT) :: best_mem
    double precision, dimension(NP), intent(OUT) :: obfval
    double precision, dimension(NP) :: obfval_pr, phi_pr, phi
    double precision, intent(OUT) :: best_obfval
    double precision, dimension(n_x) :: rand_C1
    integer :: i, ibest, j

! Initialize a random population.
    if (iam == 0) then
        poparray = 0.0d0
        do i=1,NP
            call rnd_dist(rand_C1)
            poparray(i,:) = x_lower + rand_C1*(x_upper-x_lower)
        end do
    end if

#ifdef HAVE_MPI
! Distribute population array to the rest of the processors.
    call MPI_Bcast(poparray, NP*n_x, MPI_Double_Precision, 0, &
                    MPI_Comm_World, ierr)
    call MPI_Barrier(MPI_Comm_World, ierr)
#endif

! Evaluate fitness function to find the best member.
    obfval = 0.0d0
    obfval_pr = 0.0d0
    phi = 0.0d0
    phi_pr = 0.0d0
    do i=1+iam,NP,n_procs
        call FTN(poparray(i,:), obfval_pr(i), phi_pr(i))
    end do

#ifdef HAVE_MPI
    call MPI_Reduce(obfval_pr, obfval, NP, MPI_Double_Precision, &
                    MPI_Sum, 0, MPI_Comm_World, ierr)
    call MPI_Reduce(phi_pr, phi, NP, MPI_Double_Precision, &
                    MPI_Sum, 0, MPI_Comm_World, ierr)
#else
    obfval = obfval_pr
    phi = phi_pr
#endif

    if (iam == 0) then
        ibest=1
        best_obfval = obfval(1)
        do i=2,NP
          if (obfval(i)<best_obfval .and. phi(i) == 0.0d0) then
            ibest = i
            best_obfval = obfval(i)
          end if
        end do

        best_mem = poparray(ibest,:)
    end if

    RETURN

END SUBROUTINE init_pop
!---------------------------------------------------------------------------
       
!---------------------------------------------------------------------------
FUNCTION randperm(num)
!  =======================================================================
! |
!  =======================================================================

    implicit none

    integer, intent(in) :: num
    integer :: number, i, j, k
    integer, dimension(num) :: randperm
    double precision, dimension(num) :: rand2
    intrinsic random_number

    call rnd_dist(rand2)
    do i=1,num
       number=1
       do j=1,num
          if (rand2(i) > rand2(j)) then
             number=number+1
          end if
       end do
       do k=1,i-1
          if (rand2(i) <= rand2(k) .and. rand2(i) >= rand2(k)) then
             number=number+1
          end if
       end do
       randperm(i)=number
    end do

    RETURN

END FUNCTION randperm
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
SUBROUTINE FTN(X, objval, phi)
!  =======================================================================
! |  This subroutine handles penalties for DE.  It is the generic         |
! |  objective function routine which calls the specific model function.  |
!  =======================================================================

    implicit none

    double precision, dimension(n_x), intent(IN) :: X
    double precision, intent(OUT) :: objval, phi
    double precision :: f
    double precision, dimension(n_g) :: g
!    double precision, parameter :: pen_val = 10000d0 
!   DE_pen_val is now declared in model.f90 of the ERGO package.
    integer :: i

    phi = 0.0d0

    call model_objf(x,f,g)

! Lower objval is better, so penalties are positive when added.
    objval = f
    do i=1,n_g
       if (g(i).gt.0d0) then
        objval = objval + DE_pen_val*g(i)
        phi = phi + g(i)
       end if
    end do
    
END SUBROUTINE FTN
!---------------------------------------------------------------------------

!!---------------------------------------------------------------------------
!SUBROUTINE alloc_DE()
!!  =======================================================================
!! |  Allocate memory for DE variables.                                    |
!!  =======================================================================
!
!    implicit none
!
!    integer :: astat
!
!    astat = 0
!    allocate(bestmem_XC(n_x,NQ), STAT=astat)
!    if (astat > 0) call error_msg(1,"bestmem_XC")
!
!    RETURN
!
!END SUBROUTINE alloc_DE
!!---------------------------------------------------------------------------
!
!!---------------------------------------------------------------------------
!SUBROUTINE dealloc_DE()
!!  =======================================================================
!! |  Deallocate memory for DE variables.                                  |
!!  =======================================================================
!
!    implicit none
!
!    integer :: astat
!
!    astat = 0
!    deallocate(bestmem_XC, STAT=astat)
!    if (astat > 0) call error_msg(2,"bestmem_XC")
!
!    RETURN
!
!END SUBROUTINE dealloc_DE
!!---------------------------------------------------------------------------
END MODULE DE
