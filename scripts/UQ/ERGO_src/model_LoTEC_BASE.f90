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
    use LoTEC_assim

    implicit none

    integer(i4b) :: nop
    real(DP), dimension(:), allocatable :: x_lower, x_upper
    real(DP) :: DE_pen_val

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
    
    integer i, j, option, usesite, thissite
    character(len=6) site_code(nsites_max)
    character(len=80) fname
    logical heterosk
    double precision step(N_parms_max)
    character(len=4) nstr4
    character(len=200) dummy

    !------------ User-specified model run options ------------------
    steadystate  =  1         !LoTEC steady state mode
    zobler_code  =  1         
    nsoillayers(:)  =  14     !number of specified soil layers
    varyparms    = .false.
    writemodel   = .false.   
    !constraint options
    option       = #OPTION#   !1 = NEE/LE, 2 = NEE/LE/nobio 3 = NEE/nobio
    usenee = .true.
    usegpp = .false.
    usele  = .true.
    usebio = .false.
    autocorr     = .true.    !use autocorrelation in likelihood
    heterosk     = .false.    !use heteroskedastic errors in likelihood
    distresid    =  1         !0 = Gaussian error, 1 = Laplace error
    if (option .eq. 1) then   !1. NEE+LE+biometric
       usele  = .true.
       usebio = .true.
    else if (option .eq. 2) then   !2. NEE+LE only
       usele  = .true.
       usebio = .false.
    else if (option .eq. 3) then   !3. NEE only
       usele = .false.
       usebio = .false.
    else if (option .eq. 4) then  
       distresid = 1
       heterosk = .true.
       autocorr = .true.
    end if
    useerrparms  = .true.     !estimate error (stddev) as parameter
    oneacerr     = .true.     !same AC and stdevs for all sites
    nobtypes     = -1*(usenee+usegpp+usele)
    !other options
    dynlai       = .false.    !dynamic LAI model (based on stem carbon)
    usephen      = .true.     !use phenology model (rather than LAI forcing)
    
    !-----------------------------------------------------------

    ! Default value for model sizes.
    !load parameter range file to get number of parameters

    nsites = 1 !change number of sites HERE

    !get site info from siteinfo.txt
    write(nstr4, '(I4)'), thissite+1000
    open(unit=10, file = '/home/zdr/models/LoTEC/assim/input/siteinfo.txt')
    read(10,*), dummy
    j=0

    optsites(1) = #SITE1#
    optsites(2) = 13
    optsites(3) = 18
    optsites(4) = 24
    optsites(5) = 15 

    do i=1,30  !nsites_max
       read(10,30), usesite,siteno(j+1),site_code(j+1),sitename(j+1),startyear(j+1), nyears(j+1), LAI(j+1), cstem(j+1), croot(j+1), csoil(j+1), fsoil(j+1), LMA(j+1), sand(j+1), clay(j+1), fracd(j+1), fracc(j+1), fracg(j+1)
       if (i .eq. optsites(j+1) .and. j .lt. nsites) then 
          !if (usesite .eq. 1) then 
          !   print*, 'reading data for ', sitename(j+1)
          !   j=j+1
          !end if
          print*, 'reading data for ', sitename(j+1)
          j=j+1
       end if
    end do

30  format(I1,3x,I2,2x,A6,2x,A25,I4,3x,I2,2x,f3.1,1x,f5.1,1x,f5.2,1x,f5.1,1x,f5.0,1x,f3.0,1x,2(f5.2,1x),2(f5.3,2x),f5.3)
    print*, nsites, ' total sites'

    !determine conifer or deciduous/grass
    if (fracd(1) .gt. 0.5 .or. fracg(1) .gt. 0.5) then 
       isconifer = .false.
       if (fracg(1) .gt. 0.5) isgrass = .true.
    else
       isconifer = .true.
    end if
    
    fname = '/home/zdr/models/LoTEC/assim/input/parm_ranges.dat'
    open(unit=9, FILE = fname, STATUS = 'OLD')
 
    !load parameter ranges and defaults
    do i=1,N_parms_LoTEC
       read(9,90), parm_ranges(1:2,i), parms(1,i), step(i), priortype(i), useparms(i), pnames(i)
       parms(:,i) = parms(1,i)           
    end do
90  format (4(e8.2,2X),A1,2x,I1,A200)    
    close(unit = 9) 

    !get number of parameters to optimize from parm_ranges.txt
    N_parms_opt=0

    if (autocorr) then 
       useparms(93:93+nobtypes-1)=1
    else
       useparms(93:93+nobtypes-1)=0
    end if
    if (heterosk) then
       useparms(96:96+nobtypes-1)=1
       useparms(99:99+nobtypes-1)=1
       useparms(102)=1
       useparms(105)=1
    else
       useparms(96:96+nobtypes-1)=1
       useparms(99:107)=0
       parms(:,99:107)=0.0
    end if
    !useparms(25)=2
    if (usegpp) then 
       useparms(71:86)=0
    end if
    if (isconifer) then 
       useparms(10)=1      !leaf time constant
       useparms(85)=1      !Ts threshold
       useparms(65:67) = 0 !deciduous phenology parameters
       useparms(90:92) = 0
    else
       useparms(10)=0
       useparms(85)=0
       useparms(65:67)=1
       if (nsites .gt. 1) useparms(90:92)=1
       if (isgrass) then 
          if (optsites(1) .eq. 30) parms(:,1) = 4.0  !C4 grass/crop
          parms(:,6)     = 0.30  !allocation
          parms(:,7:8)   = 0.00
          parms(:,9)     = 0.70
          parms(:,72:73) = 0.00  !coarse/stem carbon
          useparms(72:73) = 0
          useparms(10)    = 1    !leaf turnover
          useparms(13)    = 1    !fine root turnover
       end if
    end if

    do i=1,N_parms_LoTEC
       if (useparms(i) .eq. 1) N_parms_opt=N_parms_opt+1   
       if (useparms(i) .eq. 2) N_parms_opt=N_parms_opt+nsites
    end do
    
    print*, N_parms_opt, 'parameters optimized'

    n_param = N_parms_opt
    n_param_pre = 0
    n_pen = 0
    n_pen_pre = 0
    n_conv = n_param
    n_conv_pre = n_param_pre
    

    RETURN

END SUBROUTINE model_size
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
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
  real ptemp, parmstemp(N_parms_max), post
  double precision vector(2), std
  integer :: i, j, s, astat, thissite, flag
  logical isreal
  character(len=80) fname

  character(len=6) site_code(nsites_max)
  character(len=4) nstr4
  character(len=200) dummy
  integer usesite, pct, ntoadd, ylower, yupper
  logical dopseudo

  dopseudo = .false.

  print*, s, startyear(1)
  ylower = startyear(1)
  yupper = startyear(1)+nyears(1)-1
  call get_LoTEC_data(1, ylower, yupper)

  !psuedodata experiment:  overwrite constraints with model output, add random error
  !(heteroskedastic)  Default parameters are used.
  if (dopseudo .eq. .true.) then 
     call getposterior(post, flag)
     constraints(3,:) = flux_out(1,:)
     constraints(7,:) = flux_out(3,:)
     constraints(:,1:8784)=-999
     do i=1,nhours(1)
        call RANDOM_NUMBER (vector)
        !print*, vector
        std = sqrt(-2*log(vector(1)))*cos(2*3.1416*vector(2))
        if (constraints(3,i) .ne. -999) then 
           if (constraints(3,i) .lt. 0) then 
              constraints(3,i) = constraints(3,i)+std*(1.42+0.19*constraints(3,i))
           else
              constraints(3,i) = constraints(3,i)+std*(0.62+0.63*constraints(3,i))
           end if
        end if
        if (constraints(7,i) .ne. -999) then 
           constraints(7,i) = constraints(7,i)+std*(15.3+0.23*constraints(7,i))
        end if
     end do
     csoil(1) = sum(cpools_out(5:9,12000))
     cstem(1) = cpools_out(2,12000)
     croot(1) = cpools_out(3,12000)
     fsoil(1) = sum(flux_out(6,8761:17520))*12*3600/1e6
  end if

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

  J=1
  !load parameter ranges

  pct=1

  do i=1,N_parms_LoTEC
     ntoadd=0
     if (useparms(i) .eq. 1) then
        ntoadd=1
     endif
     if (useparms(i) .eq. 2) then 
        ntoadd=nsites
     endif
     if (ntoadd .gt. 0) then 
        if (priortype(i) .eq. 'U') then 
           x_lower(pct:pct+ntoadd-1)=parm_ranges(1,i)
           x_upper(pct:pct+ntoadd-1)=parm_ranges(2,i)
        end if
        if (priortype(i) .eq. 'N') then 
           x_lower(pct:pct+ntoadd-1)=parm_ranges(1,i)-5.0*parm_ranges(2,i)
           x_upper(pct:pct+ntoadd-1)=parm_ranges(1,i)+5.0*parm_ranges(2,i)
        end if
     end if
     pct=pct+ntoadd
  end do

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
  
    integer I, J, pct
    real(DP), dimension(nop), intent(IN) :: x
    real(DP), intent(OUT) :: f
    real(DP), dimension(n_g), intent(OUT) :: g
    real(DP) :: posterior
    integer :: nfe
    integer flag

! Initialize output.
    f = 0d0
    g = 0d0

! Default number of function evaluations is 1 for each call to model_objf.
    nfe = 1
    pct=1
 
    do i=1,N_parms_LoTEC
       if (useparms(i) .eq. 1) then
          parms(:,i) = x(pct)
          pct=pct+1
       end if
       if (useparms(i) .eq. 2) then 
          parms(:,i) = x(pct:pct+nsites-1)
          pct=pct+nsites
       end if 
    end do
    
    call getposterior(posterior, flag)

    !if (ISNAN(post)) then
    !   f=-99999999999
    !end if
    f = -1.0*posterior     !use negative log likelihood
    nfeval = nfeval + nfe

    print*, nfeval, f

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

    real(DP), dimension(:), intent(IN) :: Solution
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

    RETURN

END SUBROUTINE model_finalize
!---------------------------------------------------------------------------
END MODULE model
