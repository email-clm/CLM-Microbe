program qpso

implicit none
include "mpif.h"

integer maxiter, maxpop, maxparms
parameter (maxiter = 10000)
parameter (maxpop = 2048)
parameter (maxparms = 512)

integer i, j,k, npop, nparms, niter, gInx, nfunc(maxpop)
integer nfuncall(maxpop)
integer np, myid, ierr, start_iter
integer pft(maxparms)
double precision gbest(maxparms)
double precision mbest(maxparms)
double precision feval, beta_l, beta_u, beta 
double precision betapro(maxparms), pupdate(maxparms)
double precision pbest(maxpop, maxparms)
double precision pbestall(maxpop, maxparms)
double precision f_x(maxpop), x(maxpop, maxparms)
double precision xall(maxpop, maxparms)
double precision f_pbest(maxpop), f_pbestall(maxpop), f_gbest
double precision gpar(maxiter,maxparms)
double precision gobj(maxiter)
double precision pmin(maxparms), pmax(maxparms)
double precision fi(maxparms), u(maxparms), v(maxparms)
logical isvalid, restart
character(len=4) popst
character(len=100) mymachine, dummy, parm_name(maxparms), case_name, thisfmt

!------- user-tunable QPSO algorithm parameters ----------------

npop = 64         !number of particles
niter = 200       !number of iterations
beta_l = 0.4d0
beta_u = 0.7d0
restart = .false.
!mymachine = 'oic'   !Fill in correct machine and uncomment this line

!---------------------------------------------------------------

call mpi_init(ierr)
call mpi_comm_size(mpi_comm_world, np, ierr)
call mpi_comm_rank(mpi_comm_world, myid, ierr)

!get parameter information from the parm_list file
if (myid .eq. 0) then 
  open(unit = 8, status='old', file = './parm_list')
  read(8,*) case_name 
  nparms=0
  do i=1,maxparms
     read(8,*, end=10) parm_name(i), pft(i), pmin(i), pmax(i)
     nparms = nparms+1
  end do
end if

10 continue
if (myid .eq. 0) then 
  close(8)
  print*, nparms, ' Parameters optimized'
end if

!broadcast parameter info to other procs
call mpi_bcast(nparms, 1, mpi_integer, 0, mpi_comm_world, ierr)
call mpi_bcast(pmin, maxparms, mpi_double, 0, mpi_comm_world, ierr)
call mpi_bcast(pmax, maxparms, mpi_double, 0, mpi_comm_world, ierr)


nfunc(:) = 0   !keep track of total function evaluations 

x(:,:) = 0d0
f_x(:) = 0d0
f_pbest(:) = 0d0

if (restart .eqv. .false.) then 
   do i=myid+1,npop,np      
      !randomize starting locations
      call init_random_seed
      call random_number(u)
      x(i,:) = pmin + (pmax-pmin) * u
      f_x(i) = feval(x(i,:), nparms, i, mymachine)
      nfunc(i) = nfunc(i)+1
      f_pbest(i) = f_x(i)
   end do
   
   call mpi_allreduce(x, xall, maxparms*maxpop, mpi_double, mpi_sum, &
        mpi_comm_world, ierr)
   call mpi_allreduce(f_pbest, f_pbestall, maxpop, mpi_double, mpi_sum, &
        mpi_comm_world, ierr)
   pbestall = xall

   !initialize pbest and gbest
   if (myid .eq. 0) then 
      gInx = 1
      do i=2,npop
         if (f_pbestall(i) .lt. f_pbestall(gInx)) gInx = i
      end do
      gbest = pbestall(gInx,:)
      f_gbest = f_pbestall(gInx)
   end if
   call mpi_bcast(gbest, maxparms, mpi_double, 0, mpi_comm_world, ierr)
   call mpi_bcast(f_gbest, 1, mpi_double, 0, mpi_comm_world, ierr)
   start_iter = 1
else
   !load restart information 
   if (myid .eq. 0) then 
      open(unit=8, file='./qpso_restart.txt')
      read(8,*) start_iter
      do j=1,npop
         read(8,*) xall(j,1:nparms)
         read(8,*) pbestall(j,1:nparms)
         read(8,*) f_pbestall(j)
      end do
      read(8,*) gbest(1:nparms)
      read(8,*) f_gbest
   end if
   xall=pbestall   
   call mpi_bcast(xall, maxparms*maxpop, mpi_double, 0, mpi_comm_world, ierr) 
   call mpi_bcast(pbestall, maxparms*maxpop, mpi_double, 0, mpi_comm_world, ierr)
   call mpi_bcast(f_pbestall, maxpop, mpi_double, 0, mpi_comm_world, ierr)
   call mpi_bcast(gbest, maxparms, mpi_double, 0, mpi_comm_world, ierr)
   call mpi_bcast(f_gbest, 1, mpi_double, 0, mpi_comm_world, ierr)
end if



!QPSO algorithm
do i=start_iter,niter
   if (myid .eq. 0) print*, 'Iteration', i
   beta = beta_u - (beta_u-beta_l)*i/niter
   !compute mean of best parameters (all procs)
   do k=1, nparms
      mbest(k) = sum(pbestall(1:npop,k))/npop
   end do
   !print*, mbest(1:nparms)
   !MPI over population
   x(:,:) = 0d0
   pbest(:,:) = 0d0
   f_pbest(:) = 0d0
   do j = myid+1,npop,np
      isvalid = .false.
      do while (isvalid .eqv. .false.)   !Force a set of parameters within the bounds
         call random_number(fi)
         call random_number(u)
         call random_number(v)

         isvalid=.true.
         do k=1,nparms
            pupdate = fi(k)*pbestall(j,k) + (1-fi(k))*gbest(k)
            betapro = beta * abs(mbest(k)-xall(j,k))

            x(j,k) = pupdate(k)+((-1d0)**ceiling(0.5+v(k)))*betapro(k)*(-log(u(k)))

            if (x(j,k) .lt. pmin(k) .or. x(j,k) .gt. pmax(k)) isvalid=.false.
            !DMR for better load balancing, instead of not running, run at the boundary.  Testing only.
            !if (x(j,k) .lt. pmin(k)) x(j,k) = pmin(k)
            !if (x(j,k) .gt. pmax(k)) x(j,k) = pmax(k)
         end do
      end do

      !run the model to get the cost function
      f_x(j) = feval(x(j,:),nparms, j, mymachine)
      nfunc(j) = nfunc(j)+1
    
20 continue
      if (f_x(j) .lt. f_pbestall(j)) then 
         pbest(j,:) = x(j,:)
         f_pbest(j) = f_x(j)
      else 
         pbest(j,:) = pbestall(j,:)
         f_pbest(j) = f_pbestall(j)
      end if
   end do 
      
   call mpi_allreduce(pbest, pbestall, maxparms*maxpop, mpi_double, mpi_sum, &
        mpi_comm_world, ierr)
   call mpi_allreduce(x, xall, maxpop, mpi_double, mpi_sum, &
        mpi_comm_world, ierr)
   call mpi_allreduce(f_pbest, f_pbestall, maxpop, mpi_double, mpi_sum, &
        mpi_comm_world, ierr)

   !update overall best (all procs)
   do j=1,npop
      if (f_pbestall(j) .lt. f_gbest) then 
         gbest = pbestall(j,:)
         f_gbest = f_pbestall(j)
      end if
   end do

   !save info from this iteration
   gpar(i,:) = gbest
   gobj(i) = f_gbest
   call mpi_allreduce(nfunc,nfuncall, maxpop, mpi_integer, mpi_sum, &
        mpi_comm_world, ierr)
   if (myid .eq. 0) then 
     open(unit=8, file='qpso_best.txt')
     write(8,*), 'Iteration', i
     write(8,*), 'Objective function:', gobj(i)
     do k=1,nparms
        write(8,'(A,1x,I2,1x,g13.6)') trim(parm_name(k)), pft(k), gpar(i,k)
     end do
     close(8)
     if (i .eq. 1) then
        open(unit=10, file='qpso_costfunc.txt')
     else
        open(unit=10, file='qpso_costfunc.txt', status='old', position='append', action='write')
     end if
     write(10,*) i, sum(nfuncall) , gobj(i)
     close(10)

     !write the restart file
     write(popst, '(I4)') nparms
     thisfmt = '(' // trim(popst) // '(g13.6,1x))'
     open(unit=11, file = 'qpso_restart.txt')
     write(11,'(I4)') i        !current iteration number
     do j=1,npop
        write(11,fmt=trim(thisfmt)) xall(j,1:nparms)      !current parameters for each population
        write(11,fmt=trim(thisfmt)) pbestall(j,1:nparms)  !best parameters for each population 
        write(11,'(g13.6)') f_pbestall(j)                 !best objective function for each population
     end do
     write(11,fmt=trim(thisfmt)) gbest(1:nparms)          !overall best parameters
     write(11,'(g13.6)') f_gbest                          !overall best objectivefunction
     close(11)
  end if
end do
call mpi_finalize(ierr)

end program qpso


!Function to evaluate the CLM/ALM model
double precision function feval(parms, nparms, thispop, mymachine)

integer nparms, i, thispop
double precision parms(500), trueparms(4)
double precision mydata(1000), model(1000), sse(1000)
double precision temp(1000), par(1000)
character(len=6) thispopst
character(len=100) mymachine


write(thispopst, '(I6)') 100000+thispop

!write the parameters to file 
open(unit=9, file='./parm_data_files/parm_data_' // thispopst(2:6))
do i=1,nparms
   write(9,*) parms(i)
end do
close(9)

!Call python workflow to set up and launch model simulation
call system('python UQ_runens.py --ens_num ' // thispopst(2:6) // &
     ' --parm_list ' // 'parm_list --parm_data ./parm_data_files/' // &
     'parm_data_' // thispopst(2:6) // ' --constraints constraints' // &
     ' --machine ' // trim(mymachine))

!get the sum of squared errors
open(unit=9, file='./ssedata/mysse_' // thispopst(2:6) // '.txt')
read(9,*) feval
close(9)

return 

end function feval
