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
    logical, dimension(:), allocatable :: islog

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
    
   integer i
   character(len=200) dummy, fname
   real input(4)

   !get the nubmer of CLM parameters
   fname = './parm_list'
   open(unit=8, file=fname, status='old')
   read(8,*) dummy

   n_param=0
   do i=1,500
      read(8,*,end=10) dummy
      n_param=n_param+1
   end do

10 continue
   close(8)

   print*, n_param, 'parameters optimized'

   n_param_pre = 0
   n_pen = 0
   n_pen_pre = 0
   n_conv = n_param
   n_conv_pre = n_param_pre
     
   
   RETURN

 END SUBROUTINE model_size

SUBROUTINE model_init(eval_index, prereq_in)
  !  ========================================================================
  ! |  This routine does any initialization that the model may require.  Any |
  ! |  prerequisite solution array may be included in the second argument.   |
  ! |  The output array is a list of parameter indices to be included by the |
  ! |  optimization program.
  ! |  The default is to use all parameters for convergence evaluation.      |
  !  ========================================================================

  implicit none

  integer(i4b), dimension(:), intent(OUT) :: eval_index
  real(DP), dimension(:), intent(IN), optional :: prereq_in

  !real parm_ranges(2,200)
  real(DP) parm_def(200)
  integer useparm(200)
 
  integer :: i, j, pct, s, astat, thissite, flag !, n_param

  !character(len=6) mysite
  character(len=200) dummy, fname !, CLM_base, basecase, pftfile_in
  !integer lst, npd, mypft, startyr, nyears

  !get the nubmer of CLM parameters
  fname = './parm_list'
  open(unit=8, file=fname)
  read(8,*) dummy
  n_param=0
  do i=1,200
     read(8,*,end=10) dummy
     n_param=n_param+1
  end do
10 continue
  close(8)
 
  ! Default value for pen_val in DE.
  DE_pen_val = 10000d0

  ! Initialize evaluation index with default values.
  eval_index = (/ (i, i=1, size(eval_index)) /)

  if (present(prereq_in)) then
     nop = n_param
     N_G = n_pen
     ! Allocate bounds.
     call model_allocate(nop)
     ! Initialize the model_siodel for the main optimization.
  else
     nop = n_param_pre
     n_g = n_pen_pre
     ! Allocate bounds.
     call model_allocate(nop)
     ! Initialize the model for the pre-calculation.
  end if

  open(unit=8, file=fname)
  read(8,*) dummy
  do i=1,nop
     islog(i) = .false.
     read(8,*) dummy, useparm(i), x_lower(i), x_upper(i)
     if (x_lower(i) .gt. 0 .and. x_upper(i) .gt. 0) then 
       if (log10(x_upper(i)) - log10(x_lower(i)) .gt. 2) then 
         !if range spans more than two orders of magnitude and is
         !postive/nonzero, use log distribution
         x_lower(i) = log10(x_lower(i))
         x_upper(i) = log10(x_upper(i))
         islog(i) = .true.
       end if
     end if
  end do
  close(8)


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
  
    integer i, j, pct
    real(DP), dimension(nop), intent(IN) :: x
    real(DP), intent(OUT) :: f
    real(DP), dimension(n_g), intent(OUT) :: g
    real(DP) :: posterior, tmp
    integer :: nfe
    integer :: flag, try
    character(len=800) dummy, fname, pftfile, surffile
    character(len=200) myrundir, case
    character(len=6) grouptag, runtag_temp
    !integer lst, npd, mypft, startyr, nyears, ident
    double precision sse

    ! Initialize output.
    f = 0d0
    g = 0d0

    nfe = 1  !number of function evaluations
    pct = 1  !current parameter index

    !get tags for each run directory
    write(grouptag, '(I6)') 100000+iam+1
    call system('pwd')
    print*, 'IAM', iam, grouptag    

   !Default number of function evaluations is 1 for each call to model_objf.
    nfe = 1
   
    !write parameters to parm_list file (specific for each proc)
    open(unit=9, file='./parm_data_files/parm_data_' // grouptag(2:6))
    do i=1,nop
      if (islog(i)) then 
         write(9,*) 10**(x(i))
      else
         write(9,*) x(i)
      end if
    end do
    close(9)

    !Call python script to manage simulation
    call system('python UQ_runens.py --ens_num ' // grouptag(2:6) // &
         ' --parm_list ' // 'parm_list --parm_data ./parm_data_files/' // &
         'parm_data_' // grouptag(2:6) // ' --constraints ' // &
         '/home/zdr/models/CLM_SPRUCE/scripts/UQ/constraints')

    !get the sum of squared errors
    open(unit=9, file='./ssedata/mysse_' // grouptag(2:6) // '.txt')
    read(9,*) f
    close(9)

    !keep track of number of function evals
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
SUBROUTINE model_out(Solution, out_path, obf_out)
!  ========================================================================
! |  This output routine is called by the optimization program as it       |
! |  progresses toward convergence on a solution.  It allows the model to  |
! |  report the state of its various parameters given the latest best      |
! |  solution.  It is also called at the end of the optimization process.  |
!  ========================================================================

    implicit none

    real(DP), dimension(nop), intent(IN) :: Solution
    real(DP), intent(OUT), OPTIONAL :: obf_out
    character(len=50) :: fname, out_path 
    real(DP), dimension(n_g) :: g
    real(DP) :: fout
    integer :: FID = 10
    integer :: istat, i

    fname = trim(out_path) // "/solution_out"
    print*, fname
    call model_objf(Solution,fout,g)

    open(UNIT=FID, FILE=fname, STATUS='REPLACE', ACTION='WRITE', &
        ACCESS='SEQUENTIAL', IOSTAT=istat)

    write(FID,"(a,g12.5)") '# g04:  f(x) = ', fout
    write(FID,"(a)") '# i    x(i)'

    do i = 1,size(Solution)
        write(FID, "(i3,g17.10)") i, Solution(i)
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
    astat = 0
    allocate(islog(n_bounds), STAT=astat)
    if (astat > 0) call error_msg(1,"islog") 

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
    astat = 0
    deallocate(islog, STAT=astat)
    if (astat > 0) call error_msg(2,"islog")

    RETURN

END SUBROUTINE model_finalize
!---------------------------------------------------------------------------
END MODULE model
