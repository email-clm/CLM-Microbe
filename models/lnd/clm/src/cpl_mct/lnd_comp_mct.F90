module lnd_comp_mct
  
!---------------------------------------------------------------------------
!BOP
!
! !MODULE: lnd_comp_mct
!
!  Interface of the active land model component of CESM the CLM (Community Land Model)
!  with the main CESM driver. This is a thin interface taking CESM driver information
!  in MCT (Model Coupling Toolkit) format and converting it to use by CLM.
!
! !DESCRIPTION:
!
! !USES:
  use shr_kind_mod     , only : r8 => shr_kind_r8
  use shr_sys_mod      , only : shr_sys_flush
  use mct_mod          , only : mct_aVect, mct_gsmap
!
! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  SAVE
  private                              ! By default make data private
!
! !PUBLIC MEMBER FUNCTIONS:
!
  public :: lnd_init_mct               ! clm initialization
  public :: lnd_run_mct                ! clm run phase
  public :: lnd_final_mct              ! clm finalization/cleanup
!
! !PUBLIC DATA MEMBERS: None
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
! Dec/18/2009 Make sure subroutines have documentation. Erik Kluzek
!
! !PRIVATE MEMBER FUNCTIONS:
  private :: lnd_SetgsMap_mct         ! Set the land model MCT GS map
  private :: lnd_domain_mct           ! Set the land model domain information
  private :: lnd_export_mct           ! export land data to CESM coupler
  private :: lnd_import_mct           ! import data from the CESM coupler to the land model
  private :: sno_export_mct
  private :: sno_import_mct
!
! !PRIVATE DATA MEMBERS:
!
! Time averaged flux fields
!  
  type(mct_aVect)   :: l2x_l_SNAP     ! Snapshot of land to coupler data on the land grid
  type(mct_aVect)   :: l2x_l_SUM      ! Summation of land to coupler data on the land grid

  type(mct_aVect)   :: s2x_s_SNAP     ! Snapshot of sno to coupler data on the land grid
  type(mct_aVect)   :: s2x_s_SUM      ! Summation of sno to coupler data on the land grid 
!
! Time averaged counter for flux fields
!
  integer :: avg_count                ! Number of times snapshots of above flux data summed together
  integer :: avg_count_sno
!
! Atmospheric mode  
!
  logical :: atm_prognostic           ! Flag if active atmosphere component or not

!EOP
!===============================================================
contains
!===============================================================

!---------------------------------------------------------------------------
!BOP
!
! !IROUTINE: lnd_init_mct
!
! !INTERFACE:
  subroutine lnd_init_mct( EClock, cdata_l, x2l_l, l2x_l, &
                                   cdata_s, x2s_s, s2x_s, &
                                   NLFilename )
!
! !DESCRIPTION:
! Initialize land surface model and obtain relevant atmospheric model arrays
! back from (i.e. albedos, surface temperature and snow cover over land).
!
! !USES:
    use abortutils       , only : endrun
    use clm_time_manager , only : get_nstep, get_step_size, set_timemgr_init, &
                                  set_nextsw_cday
    use clm_atmlnd       , only : clm_l2a
    use clm_glclnd       , only : clm_s2x
    use clm_initializeMod, only : initialize1, initialize2
    use clm_varctl       , only : finidat,single_column, set_clmvarctl, iulog, noland, &
                                  inst_index, inst_suffix, inst_name, &
                                  create_glacier_mec_landunit 
    use clm_varorb       , only : eccen, obliqr, lambm0, mvelpp
    use controlMod       , only : control_setNL
    use decompMod        , only : get_proc_bounds
    use domainMod        , only : ldomain
    use shr_kind_mod     , only : r8 => shr_kind_r8
    use shr_file_mod     , only : shr_file_setLogUnit, shr_file_setLogLevel, &
                                  shr_file_getLogUnit, shr_file_getLogLevel, &
                                  shr_file_getUnit, shr_file_setIO
    use seq_cdata_mod    , only : seq_cdata, seq_cdata_setptrs
    use seq_timemgr_mod  , only : seq_timemgr_EClockGetData
    use seq_infodata_mod , only : seq_infodata_type, seq_infodata_GetData, seq_infodata_PutData, &
                                  seq_infodata_start_type_start, seq_infodata_start_type_cont,   &
                                  seq_infodata_start_type_brnch
    use seq_comm_mct     , only : seq_comm_suffix, seq_comm_inst, seq_comm_name
    use spmdMod          , only : masterproc, spmd_init
    use clm_varctl       , only : nsrStartup, nsrContinue, nsrBranch
    use clm_cpl_indices  , only : clm_cpl_indices_set, nflds_l2x
    use seq_flds_mod
    use mct_mod
    use ESMF
    implicit none
!
! !ARGUMENTS:
    type(ESMF_Clock),           intent(in)    :: EClock           ! Input synchronization clock
    type(seq_cdata),            intent(inout) :: cdata_l          ! Input land-model driver data
    type(mct_aVect),            intent(inout) :: x2l_l, l2x_l     ! land model import and export states
    type(seq_cdata),            intent(inout) :: cdata_s          ! Input snow-model (land-ice) driver data
    type(mct_aVect),            intent(inout) :: x2s_s, s2x_s     ! Snow-model import and export states
    character(len=*), optional, intent(in)    :: NLFilename       ! Namelist filename to read
!
! !LOCAL VARIABLES:
    integer                          :: LNDID	     ! Land identifyer
    integer                          :: mpicom_lnd   ! MPI communicator
    type(mct_gsMap),         pointer :: GSMap_lnd    ! Land model MCT GS map
    type(mct_gGrid),         pointer :: dom_l        ! Land model domain
    type(mct_gsMap),         pointer :: GSMap_sno
    type(mct_gGrid),         pointer :: dom_s
    type(seq_infodata_type), pointer :: infodata     ! CESM driver level info data
    integer  :: lsize                                ! size of attribute vector
    integer  :: g,i,j                                ! indices
    integer  :: dtime_sync                           ! coupling time-step from the input synchronization clock
    integer  :: dtime_clm                            ! clm time-step
    logical  :: exists                               ! true if file exists
    logical  :: atm_aero                             ! Flag if aerosol data sent from atm model
    logical  :: samegrid_al                          ! true if atmosphere and land are on the same grid
    integer  :: nstep
    real(r8) :: scmlat                               ! single-column latitude
    real(r8) :: scmlon                               ! single-column longitude
    real(r8) :: nextsw_cday                          ! calday from clock of next radiation computation
    character(len=SHR_KIND_CL) :: caseid             ! case identifier name
    character(len=SHR_KIND_CL) :: ctitle             ! case description title
    character(len=SHR_KIND_CL) :: starttype          ! start-type (startup, continue, branch, hybrid)
    character(len=SHR_KIND_CL) :: calendar           ! calendar type name
    character(len=SHR_KIND_CL) :: hostname           ! hostname of machine running on
    character(len=SHR_KIND_CL) :: version            ! Model version
    character(len=SHR_KIND_CL) :: username           ! user running the model
    integer :: nsrest                                ! clm restart type
    integer :: perpetual_ymd                         ! perpetual date
    integer :: ref_ymd                               ! reference date (YYYYMMDD)
    integer :: ref_tod                               ! reference time of day (sec)
    integer :: start_ymd                             ! start date (YYYYMMDD)
    integer :: start_tod                             ! start time of day (sec)
    integer :: stop_ymd                              ! stop date (YYYYMMDD)
    integer :: stop_tod                              ! stop time of day (sec)
    logical :: brnch_retain_casename                 ! flag if should retain the case name on a branch start type
    logical :: perpetual_run                         ! flag if should cycle over a perpetual date or not
    integer :: lbnum                                 ! input to memory diagnostic
    integer :: shrlogunit,shrloglev                  ! old values for log unit and log level
    integer :: begg, endg
    character(len=32), parameter :: sub = 'lnd_init_mct'
    character(len=*),  parameter :: format = "('("//trim(sub)//") :',A)"
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

    ! Set cdata data

    call seq_cdata_setptrs(cdata_l, ID=LNDID, mpicom=mpicom_lnd, &
         gsMap=GSMap_lnd, dom=dom_l, infodata=infodata)

    ! Determine attriute vector indices

    call clm_cpl_indices_set()

    ! Initialize clm MPI communicator 

    call spmd_init( mpicom_lnd, LNDID )

#if (defined _MEMTRACE)
    if(masterproc) then
       lbnum=1
       call memmon_dump_fort('memmon.out','lnd_init_mct:start::',lbnum)
    endif
#endif                      

    inst_name   = seq_comm_name(LNDID)
    inst_index  = seq_comm_inst(LNDID)
    inst_suffix = seq_comm_suffix(LNDID)

    ! Initialize io log unit

    call shr_file_getLogUnit (shrlogunit)
    if (masterproc) then
       inquire(file='lnd_modelio.nml'//trim(inst_suffix),exist=exists)
       if (exists) then
          iulog = shr_file_getUnit()
          call shr_file_setIO('lnd_modelio.nml'//trim(inst_suffix),iulog)
       end if
       write(iulog,format) "CLM land model initialization"
    else
       iulog = shrlogunit
    end if

    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (iulog)
    
    ! Use infodata to set orbital values

    call seq_infodata_GetData( infodata, orb_eccen=eccen, orb_mvelpp=mvelpp, &
         orb_lambm0=lambm0, orb_obliqr=obliqr )

    ! Consistency check on namelist filename	

    call control_setNL("lnd_in"//trim(inst_suffix))

    ! Initialize clm
    ! initialize1 reads namelist, grid and surface data (need this to initialize gsmap) 
    ! initialize2 performs rest of initialization	

    call seq_timemgr_EClockGetData(EClock,                               &
                                   start_ymd=start_ymd,                  &
                                   start_tod=start_tod, ref_ymd=ref_ymd, &
                                   ref_tod=ref_tod, stop_ymd=stop_ymd,   &
                                   stop_tod=stop_tod,                    &
                                   calendar=calendar )
    call seq_infodata_GetData(infodata, perpetual=perpetual_run,                &
                              perpetual_ymd=perpetual_ymd, case_name=caseid,    &
                              case_desc=ctitle, single_column=single_column,    &
                              scmlat=scmlat, scmlon=scmlon,                     &
                              brnch_retain_casename=brnch_retain_casename,      &
                              start_type=starttype, model_version=version,      &
                              hostname=hostname, username=username,             &
                              samegrid_al=samegrid_al                           &
                                )
    call set_timemgr_init( calendar_in=calendar, start_ymd_in=start_ymd, start_tod_in=start_tod, &
                           ref_ymd_in=ref_ymd, ref_tod_in=ref_tod, stop_ymd_in=stop_ymd,         &
                           stop_tod_in=stop_tod,  perpetual_run_in=perpetual_run,                &
                           perpetual_ymd_in=perpetual_ymd )
    if (     trim(starttype) == trim(seq_infodata_start_type_start)) then
       nsrest = nsrStartup
    else if (trim(starttype) == trim(seq_infodata_start_type_cont) ) then
       nsrest = nsrContinue
    else if (trim(starttype) == trim(seq_infodata_start_type_brnch)) then
       nsrest = nsrBranch
    else
       call endrun( sub//' ERROR: unknown starttype' )
    end if

    call set_clmvarctl(    caseid_in=caseid, ctitle_in=ctitle,                     &
                           brnch_retain_casename_in=brnch_retain_casename,         &
                           single_column_in=single_column, scmlat_in=scmlat,       &
                           scmlon_in=scmlon, nsrest_in=nsrest, version_in=version, &
                           hostname_in=hostname, username_in=username)

    ! Read namelist, grid and surface data

    call initialize1( )

    ! If no land then exit out of initialization

    if ( noland ) then
           call seq_infodata_PutData( infodata, sno_present   =.false.)
           call seq_infodata_PutData( infodata, lnd_present   =.false.)
           call seq_infodata_PutData( infodata, lnd_prognostic=.false.)
       return
    end if

    ! Determine if aerosol and dust deposition come from atmosphere component

    call seq_infodata_GetData(infodata, atm_aero=atm_aero )
    if ( .not. atm_aero )then
       call endrun( sub//' ERROR: atmosphere model MUST send aerosols to CLM' )
    end if

    ! Initialize lnd gsMap and domain

    call lnd_SetgsMap_mct( mpicom_lnd, LNDID, gsMap_lnd ) 	
    lsize = mct_gsMap_lsize(gsMap_lnd, mpicom_lnd)
    call lnd_domain_mct( lsize, gsMap_lnd, dom_l )

    ! Initialize lnd attribute vectors coming from driver

    call mct_aVect_init(x2l_l, rList=seq_flds_x2l_fields, lsize=lsize)
    call mct_aVect_zero(x2l_l)

    call mct_aVect_init(l2x_l, rList=seq_flds_l2x_fields, lsize=lsize)
    call mct_aVect_zero(l2x_l)

    call mct_aVect_init(l2x_l_SNAP, rList=seq_flds_l2x_fluxes, lsize=lsize)
    call mct_aVect_zero(l2x_l_SNAP)

    call mct_aVect_init(l2x_l_SUM , rList=seq_flds_l2x_fluxes, lsize=lsize)
    call mct_aVect_zero(l2x_l_SUM )

    if (masterproc) then
       write(iulog,format)'time averaging the following flux fields over the coupling interval'
       write(iulog,format) trim(seq_flds_l2x_fluxes)
    end if

    ! Finish initializing clm

    call initialize2()

    ! Check that clm internal dtime aligns with clm coupling interval

    call seq_timemgr_EClockGetData(EClock, dtime=dtime_sync )
    dtime_clm = get_step_size()
    if (masterproc) write(iulog,*)'dtime_sync= ',dtime_sync,&
         ' dtime_clm= ',dtime_clm,' mod = ',mod(dtime_sync,dtime_clm)
    if (mod(dtime_sync,dtime_clm) /= 0) then
       write(iulog,*)'clm dtime ',dtime_clm,' and Eclock dtime ',&
            dtime_sync,' never align'
       call endrun( sub//' ERROR: time out of sync' )
    end if

    ! Create land export state 

    call get_proc_bounds(begg, endg) 

#ifndef CPL_BYPASS
    call lnd_export_mct( clm_l2a, l2x_l, begg, endg )
#endif

    if (create_glacier_mec_landunit) then
       call seq_cdata_setptrs(cdata_s, gsMap=gsMap_sno, dom=dom_s)

       ! Initialize sno gsMap (same as gsMap_lnd)
       call lnd_SetgsMap_mct( mpicom_lnd, LNDID, gsMap_sno )
       lsize = mct_gsMap_lsize(gsMap_sno, mpicom_lnd)

       ! Initialize sno domain (same as lnd domain)
       call lnd_domain_mct( lsize, gsMap_sno, dom_s )

       ! Initialize sno attribute vectors
       call mct_aVect_init(x2s_s, rList=seq_flds_x2s_fields, lsize=lsize)
       call mct_aVect_zero(x2s_s)

       call mct_aVect_init(s2x_s, rList=seq_flds_s2x_fields, lsize=lsize)
       call mct_aVect_zero(s2x_s)

       ! In contrast to l2x_l_SNAP / l2x_l_SUM, for s2x we accumulate/average all fields,
       ! not just fluxes. This is because glc wants the time-averaged tsrf field (and the
       ! other state field, topo, is not time-varying, so it doesn't matter what we do
       ! with that field)
       call mct_aVect_init(s2x_s_SUM , rList=seq_flds_s2x_fields, lsize=lsize)
       call mct_aVect_zero(s2x_s_SUM )

       call mct_aVect_init(s2x_s_SNAP , rList=seq_flds_s2x_fields, lsize=lsize)
       call mct_aVect_zero(s2x_s_SNAP )

       ! Create mct sno export state
       call sno_export_mct(clm_s2x, s2x_s)
    endif   ! create_glacier_mec_landunit

    ! Initialize averaging counter

    avg_count = 0

    ! Fill in infodata

    call seq_infodata_PutData( infodata, lnd_prognostic=.true.)
    call seq_infodata_PutData( infodata, lnd_nx=ldomain%ni, lnd_ny=ldomain%nj)
    if (create_glacier_mec_landunit) then
       call seq_infodata_PutData( infodata, sno_present=.true.)
       call seq_infodata_PutData( infodata, sno_prognostic=.false.)
       call seq_infodata_PutData( infodata, sno_nx=ldomain%ni, sno_ny=ldomain%nj)
    else
       call seq_infodata_PutData( infodata, sno_present=.false.)
       call seq_infodata_PutData( infodata, sno_prognostic=.false.)
    endif

    call seq_infodata_GetData(infodata, nextsw_cday=nextsw_cday )

    call set_nextsw_cday( nextsw_cday )

#ifdef CPL_BYPASS
    !Calcualte next radiation calendar day (since atm model did not run to set this)
    !DMR:  NOTE this assumes a no-leap calendar and equal input/model timesteps
    nstep = get_nstep()
    nextsw_cday = mod((nstep/(86400._r8/dtime_clm))*1.0_r8,365._r8)+1._r8 
    call set_nextsw_cday( nextsw_cday )
#endif

    ! Determine atmosphere modes

    call seq_infodata_GetData(infodata, atm_prognostic=atm_prognostic)
    if (masterproc) then
       if ( atm_prognostic )then
          write(iulog,format) 'Atmospheric input is from a prognostic model'
       else
          write(iulog,format) 'Atmospheric input is from a data model'
       end if
    end if

    ! Reset shr logging to original values

    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)

#if (defined _MEMTRACE)
    if(masterproc) then
       write(iulog,*) TRIM(Sub) // ':end::'
       lbnum=1
       call memmon_dump_fort('memmon.out','lnd_int_mct:end::',lbnum)
       call memmon_reset_addr()
    endif
#endif

  end subroutine lnd_init_mct

!---------------------------------------------------------------------------
!BOP
!
! !IROUTINE: lnd_run_mct
!
! !INTERFACE:
  subroutine lnd_run_mct( EClock, cdata_l, x2l_l, l2x_l, &
                                  cdata_s, x2s_s, s2x_s)
!
! !DESCRIPTION:
! Run clm model
!
! !USES:
    use shr_kind_mod    ,only : r8 => shr_kind_r8
    use clmtype
    use clm_atmlnd      ,only : clm_l2a, clm_a2l
    use clm_driver      ,only : clm_drv
    use clm_time_manager,only : get_curr_date, get_nstep, get_curr_calday, get_step_size, &
                                advance_timestep, set_nextsw_cday,update_rad_dtime
    use decompMod       ,only : get_proc_bounds
    use abortutils      ,only : endrun
    use clm_varctl      ,only : iulog, create_glacier_mec_landunit 
    use clm_varorb      ,only : eccen, obliqr, lambm0, mvelpp
    use shr_file_mod    ,only : shr_file_setLogUnit, shr_file_setLogLevel, &
                                shr_file_getLogUnit, shr_file_getLogLevel
    use seq_cdata_mod   ,only : seq_cdata, seq_cdata_setptrs
    use seq_timemgr_mod ,only : seq_timemgr_EClockGetData, seq_timemgr_StopAlarmIsOn, &
                                seq_timemgr_RestartAlarmIsOn, seq_timemgr_EClockDateInSync
    use seq_infodata_mod,only : seq_infodata_type, seq_infodata_GetData
    use spmdMod         ,only : masterproc, mpicom
    use perf_mod        ,only : t_startf, t_stopf, t_barrierf
    use clm_glclnd      ,only : clm_s2x, clm_x2s, unpack_clm_x2s
    use shr_orb_mod     ,only : shr_orb_decl
    use clm_varorb      ,only : eccen, mvelpp, lambm0, obliqr
    use clm_cpl_indices ,only : nflds_l2x, nflds_x2l
    use mct_mod
    use ESMF
    implicit none
!
! !ARGUMENTS:
    type(ESMF_Clock) , intent(in)    :: EClock    ! Input synchronization clock from driver
    type(seq_cdata)  , intent(inout) :: cdata_l   ! Input driver data for land model
    type(mct_aVect)  , intent(inout) :: x2l_l     ! Import state to land model
    type(mct_aVect)  , intent(inout) :: l2x_l     ! Export state from land model
    type(seq_cdata)  , intent(in)    :: cdata_s   ! Input driver data for snow model (land-ice)
    type(mct_aVect)  , intent(inout) :: x2s_s     ! Import state for snow model
    type(mct_aVect)  , intent(inout) :: s2x_s     ! Export state for snow model
!
! !LOCAL VARIABLES:
    integer :: ymd_sync                   ! Sync date (YYYYMMDD)
    integer :: yr_sync                    ! Sync current year
    integer :: mon_sync                   ! Sync current month
    integer :: day_sync                   ! Sync current day
    integer :: tod_sync                   ! Sync current time of day (sec)
    integer :: ymd                        ! CLM current date (YYYYMMDD)
    integer :: yr                         ! CLM current year
    integer :: mon                        ! CLM current month
    integer :: day                        ! CLM current day
    integer :: tod                        ! CLM current time of day (sec)
    integer :: dtime                      ! time step increment (sec)
    integer :: nstep                      ! time step index
    logical :: rstwr_sync                 ! .true. ==> write restart file before returning
    logical :: rstwr                      ! .true. ==> write restart file before returning
    logical :: nlend_sync                 ! Flag signaling last time-step
    logical :: nlend                      ! .true. ==> last time-step
    logical :: dosend                     ! true => send data back to driver
    logical :: doalb                      ! .true. ==> do albedo calculation on this time step
    real(r8):: nextsw_cday                ! calday from clock of next radiation computation
    real(r8):: caldayp1                   ! clm calday plus dtime offset
    integer :: shrlogunit,shrloglev       ! old values for share log unit and log level
    integer :: begg, endg                 ! Beginning and ending gridcell index numbers
    integer :: lbnum                      ! input to memory diagnostic
    type(seq_infodata_type),pointer :: infodata ! CESM information from the driver
    type(mct_gGrid),        pointer :: dom_l    ! Land model domain data
    integer  :: g,i,lsize                       ! counters
    logical,save :: first_call = .true.         ! first call work
    logical  :: glcrun_alarm          ! if true, sno data is averaged and sent to glc this step
    logical  :: update_glc2sno_fields ! if true, update glacier_mec fields
    real(r8) :: calday                ! calendar day for nstep
    real(r8) :: declin                ! solar declination angle in radians for nstep
    real(r8) :: declinp1              ! solar declination angle in radians for nstep+1
    real(r8) :: eccf                  ! earth orbit eccentricity factor
    real(r8) :: recip                 ! reciprical
    character(len=32)            :: rdate       ! date char string for restart file names
    character(len=32), parameter :: sub = "lnd_run_mct"
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!---------------------------------------------------------------------------

#if (defined _MEMTRACE)
    if(masterproc) then
       lbnum=1
       call memmon_dump_fort('memmon.out','lnd_run_mct:start::',lbnum)
    endif
#endif

    ! Reset shr logging to my log file
    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (iulog)

    ! Determine time of next atmospheric shortwave calculation
    call seq_cdata_setptrs(cdata_l, infodata=infodata, dom=dom_l)
    call seq_timemgr_EClockGetData(EClock, &
         curr_ymd=ymd, curr_tod=tod_sync,  &
         curr_yr=yr_sync, curr_mon=mon_sync, curr_day=day_sync)
    call seq_infodata_GetData(infodata, nextsw_cday=nextsw_cday )

    call set_nextsw_cday( nextsw_cday )
    dtime = get_step_size()
#ifdef CPL_BYPASS
    !Calcualte next radiation calendar day (since atm model did not run to set this)
    !DMR:  NOTE this assumes a no-leap calendar and equal input/model timesteps
    nstep = get_nstep()
    nextsw_cday = mod((nstep/(86400._r8/dtime))*1.0_r8,365._r8)+1._r8 
    call set_nextsw_cday( nextsw_cday )
#endif

    write(rdate,'(i4.4,"-",i2.2,"-",i2.2,"-",i5.5)') yr_sync,mon_sync,day_sync,tod_sync
    nlend_sync = seq_timemgr_StopAlarmIsOn( EClock )
    rstwr_sync = seq_timemgr_RestartAlarmIsOn( EClock )

    call get_proc_bounds(begg, endg)

    ! Map MCT to land data type
    ! Perform downscaling if appropriate

    call t_startf ('lc_lnd_import')
    call lnd_import_mct( x2l_l, clm_a2l, begg, endg )
    
    ! Map to clm (only when state and/or fluxes need to be updated)

    if (create_glacier_mec_landunit) then
       update_glc2sno_fields  = .false.
       call seq_infodata_GetData(infodata, glc_g2supdate = update_glc2sno_fields)
       if (update_glc2sno_fields) then
          call sno_import_mct( x2s_s, clm_x2s )
          call unpack_clm_x2s(clm_x2s)
       endif ! update_glc2sno
    endif ! create_glacier_mec_landunit
    call t_stopf ('lc_lnd_import')

    ! Use infodata to set orbital values if updated mid-run

    call seq_infodata_GetData( infodata, orb_eccen=eccen, orb_mvelpp=mvelpp, &
         orb_lambm0=lambm0, orb_obliqr=obliqr )

    ! Loop over time steps in coupling interval

    dosend = .false.
    do while(.not. dosend)

       ! Determine if dosend
       ! When time is not updated at the beginning of the loop - then return only if
       ! are in sync with clock before time is updated

       call get_curr_date( yr, mon, day, tod )
       ymd = yr*10000 + mon*100 + day
       tod = tod
       dosend = (seq_timemgr_EClockDateInSync( EClock, ymd, tod))

       ! Determine doalb based on nextsw_cday sent from atm model

       nstep = get_nstep()
       caldayp1 = get_curr_calday(offset=dtime)
       if (nstep == 0) then
	  doalb = .false. 	
       else if (nstep == 1) then 
          doalb = (abs(nextsw_cday- caldayp1) < 1.e-10_r8) 
       else
          doalb = (nextsw_cday >= -0.5_r8) 
       end if
       call update_rad_dtime(doalb)

       ! Determine if time to write cam restart and stop

       rstwr = .false.
       if (rstwr_sync .and. dosend) rstwr = .true.
       nlend = .false.
       if (nlend_sync .and. dosend) nlend = .true.

       ! Run clm 

       call t_barrierf('sync_clm_run1', mpicom)
       call t_startf ('clm_run')
       call t_startf ('shr_orb_decl')
       calday = get_curr_calday()
       call shr_orb_decl( calday     , eccen, mvelpp, lambm0, obliqr, declin  , eccf )
       call shr_orb_decl( nextsw_cday, eccen, mvelpp, lambm0, obliqr, declinp1, eccf )
       call t_stopf ('shr_orb_decl')
       call clm_drv(doalb, nextsw_cday, declinp1, declin, rstwr, nlend, rdate)
       call t_stopf ('clm_run')

       ! Create l2x_l export state - add river runoff input to l2x_l if appropriate
       
       call t_startf ('lc_lnd_export')
#ifndef CPL_BYPASS
       call lnd_export_mct( clm_l2a, l2x_l, begg, endg )
#endif
       call t_stopf ('lc_lnd_export')

       ! Do not accumulate on first coupling freq - consistency with ccsm3

       nstep = get_nstep()
       if (nstep <= 1) then
          call mct_aVect_copy( l2x_l, l2x_l_SUM )
          avg_count = 1
       else
          call mct_aVect_copy( l2x_l, l2x_l_SNAP )
          call mct_aVect_accum( aVin=l2x_l_SNAP, aVout=l2x_l_SUM )
          avg_count = avg_count + 1
       endif
       
       ! Map sno data type to MCT

       if (create_glacier_mec_landunit) then
          call sno_export_mct(clm_s2x, s2x_s)
          if (nstep <= 1) then
             call mct_aVect_copy( s2x_s, s2x_s_SUM )
             avg_count_sno = 1
          else
             call mct_aVect_copy( s2x_s, s2x_s_SNAP )
             call mct_aVect_accum( aVin=s2x_s_SNAP, aVout=s2x_s_SUM )
             avg_count_sno = avg_count_sno + 1
          endif
       endif    ! create_glacier_mec_landunit

       ! Advance clm time step
       
       call t_startf ('lc_clm2_adv_timestep')
       call advance_timestep()
       call t_stopf ('lc_clm2_adv_timestep')

    end do

    ! Finish accumulation of attribute vector and average and zero out partial sum and counter
    
    if (avg_count /= 0) then
       recip = 1.0_r8/(real(avg_count,r8))
       l2x_l_SUM%rAttr(:,:) = l2x_l_SUM%rAttr(:,:) * recip
    endif
    call mct_aVect_copy( l2x_l_SUM, l2x_l )
    call mct_aVect_zero( l2x_l_SUM) 
    avg_count = 0                   

    if (create_glacier_mec_landunit) then
       call seq_infodata_GetData(infodata, glcrun_alarm = glcrun_alarm )
       if (glcrun_alarm) then
          if (avg_count_sno /= 0) then
             recip = 1.0_r8/(real(avg_count_sno,r8))
             s2x_s_SUM%rAttr(:,:) = s2x_s_SUM%rAttr(:,:) * recip
          endif
          call mct_aVect_copy( s2x_s_SUM, s2x_s )
          call mct_aVect_zero( s2x_s_SUM)
          avg_count_sno = 0
       endif
    endif

    ! Check that internal clock is in sync with master clock

    call get_curr_date( yr, mon, day, tod, offset=-dtime )
    ymd = yr*10000 + mon*100 + day
    tod = tod
    if ( .not. seq_timemgr_EClockDateInSync( EClock, ymd, tod ) )then
       call seq_timemgr_EclockGetData( EClock, curr_ymd=ymd_sync, curr_tod=tod_sync )
       write(iulog,*)' clm ymd=',ymd     ,'  clm tod= ',tod
       write(iulog,*)'sync ymd=',ymd_sync,' sync tod= ',tod_sync
       call endrun( sub//":: CLM clock not in sync with Master Sync clock" )
    end if
    
    ! Reset shr logging to my original values

    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)
  
#if (defined _MEMTRACE)
    if(masterproc) then
       lbnum=1
       call memmon_dump_fort('memmon.out','lnd_run_mct:end::',lbnum)
       call memmon_reset_addr()
    endif
#endif

    first_call  = .false.

  end subroutine lnd_run_mct

!---------------------------------------------------------------------------
!BOP
!
! !IROUTINE: lnd_final_mct
!
! !INTERFACE:
  subroutine lnd_final_mct( EClock, cdata_l, x2l_l, l2x_l, &
                                    cdata_s, x2s_s, s2x_s )
!
! !DESCRIPTION:
! Finalize land surface model
!
!------------------------------------------------------------------------------
!
    use seq_cdata_mod   ,only : seq_cdata, seq_cdata_setptrs
    use seq_timemgr_mod ,only : seq_timemgr_EClockGetData, seq_timemgr_StopAlarmIsOn, &
                                seq_timemgr_RestartAlarmIsOn, seq_timemgr_EClockDateInSync
    use mct_mod
    use esmf
   implicit none
! !ARGUMENTS:
    type(ESMF_Clock) , intent(in)    :: EClock    ! Input synchronization clock from driver
    type(seq_cdata)  , intent(in)    :: cdata_l   ! Input driver data for land model
    type(mct_aVect)  , intent(inout) :: x2l_l     ! Import state to land model
    type(mct_aVect)  , intent(inout) :: l2x_l     ! Export state from land model
    type(seq_cdata)  , intent(in)    :: cdata_s   ! Input driver data for snow model (land-ice)
    type(mct_aVect)  , intent(inout) :: x2s_s     ! Import state for snow model
    type(mct_aVect)  , intent(inout) :: s2x_s     ! Export state for snow model
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!---------------------------------------------------------------------------

   ! fill this in
  end subroutine lnd_final_mct

!=================================================================================

!---------------------------------------------------------------------------
!BOP
!
! !IROUTINE: lnd_SetgsMap_mct
!
! !INTERFACE:
  subroutine lnd_SetgsMap_mct( mpicom_lnd, LNDID, gsMap_lnd )
!-------------------------------------------------------------------
!
! !DESCRIPTION:
!
! Set the MCT GS map for the land model
!
!-------------------------------------------------------------------
! !USES:
    use shr_kind_mod , only : r8 => shr_kind_r8
    use decompMod    , only : get_proc_bounds, ldecomp
    use domainMod    , only : ldomain
    use mct_mod      , only : mct_gsMap, mct_gsMap_init
    implicit none
! !ARGUMENTS:
    integer        , intent(in)  :: mpicom_lnd    ! MPI communicator for the clm land model
    integer        , intent(in)  :: LNDID         ! Land model identifyer number
    type(mct_gsMap), intent(out) :: gsMap_lnd     ! Resulting MCT GS map for the land model
!
! !LOCAL VARIABLES:
    integer,allocatable :: gindex(:)  ! Number the local grid points
    integer :: i, j, n, gi            ! Indices
    integer :: lsize,gsize            ! GS Map size
    integer :: ier                    ! Error code
    integer :: begg, endg             ! Beginning/Ending grid cell index
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!---------------------------------------------------------------------------

    ! Build the land grid numbering for MCT
    ! NOTE:  Numbering scheme is: West to East and South to North
    ! starting at south pole.  Should be the same as what's used in SCRIP
    
    call get_proc_bounds(begg, endg)

    allocate(gindex(begg:endg),stat=ier)

    ! number the local grid

    do n = begg, endg
       gindex(n) = ldecomp%gdc2glo(n)
    end do
    lsize = endg-begg+1
    gsize = ldomain%ni * ldomain%nj

    call mct_gsMap_init( gsMap_lnd, gindex, mpicom_lnd, LNDID, lsize, gsize )

    deallocate(gindex)

  end subroutine lnd_SetgsMap_mct

!=================================================================================

!---------------------------------------------------------------------------
!BOP
!
! !IROUTINE: lnd_export_mct
!
! !INTERFACE:
  subroutine lnd_export_mct( clm_l2a, l2x_l, begg, endg )   
!
! !DESCRIPTION:
!
! Convert the data to be sent from the clm model to the coupler from clm data types
! to MCT data types.
! 
!---------------------------------------------------------------------------
! !USES:
    use shr_kind_mod       , only : r8 => shr_kind_r8
    use clm_varctl         , only : iulog
    use clm_time_manager   , only : get_nstep, get_step_size  
    use clm_atmlnd         , only : lnd2atm_type
    use seq_drydep_mod     , only : n_drydep
    use shr_megan_mod      , only : shr_megan_mechcomps_n
    use clm_cpl_indices
    use clmtype
    implicit none
! !ARGUMENTS:
    type(lnd2atm_type), intent(inout) :: clm_l2a    ! clm land to atmosphere exchange data type
    type(mct_aVect)   , intent(inout) :: l2x_l      ! Land to coupler export state on land grid
    integer           , intent(in)    :: begg       ! beginning grid cell index
    integer           , intent(in)    :: endg       ! ending grid cell index
!
! !LOCAL VARIABLES:
    integer  :: g,i                           ! indices
    integer  :: ier                           ! error status
    integer  :: nstep                         ! time step index
    integer  :: dtime                         ! time step   
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!---------------------------------------------------------------------------
    
    ! cesm sign convention is that fluxes are positive downward

    l2x_l%rAttr(:,:) = 0.0_r8

    do g = begg,endg
       i = 1 + (g-begg)
       l2x_l%rAttr(index_l2x_Sl_t,i)        =  clm_l2a%t_rad(g)
       l2x_l%rAttr(index_l2x_Sl_snowh,i)    =  clm_l2a%h2osno(g)
       l2x_l%rAttr(index_l2x_Sl_avsdr,i)    =  clm_l2a%albd(g,1)
       l2x_l%rAttr(index_l2x_Sl_anidr,i)    =  clm_l2a%albd(g,2)
       l2x_l%rAttr(index_l2x_Sl_avsdf,i)    =  clm_l2a%albi(g,1)
       l2x_l%rAttr(index_l2x_Sl_anidf,i)    =  clm_l2a%albi(g,2)
       l2x_l%rAttr(index_l2x_Sl_tref,i)     =  clm_l2a%t_ref2m(g)
       l2x_l%rAttr(index_l2x_Sl_qref,i)     =  clm_l2a%q_ref2m(g)
       l2x_l%rAttr(index_l2x_Sl_u10,i)      =  clm_l2a%u_ref10m(g)
       l2x_l%rAttr(index_l2x_Fall_taux,i)   = -clm_l2a%taux(g)
       l2x_l%rAttr(index_l2x_Fall_tauy,i)   = -clm_l2a%tauy(g)
       l2x_l%rAttr(index_l2x_Fall_lat,i)    = -clm_l2a%eflx_lh_tot(g)
       l2x_l%rAttr(index_l2x_Fall_sen,i)    = -clm_l2a%eflx_sh_tot(g)
       l2x_l%rAttr(index_l2x_Fall_lwup,i)   = -clm_l2a%eflx_lwrad_out(g)
       l2x_l%rAttr(index_l2x_Fall_evap,i)   = -clm_l2a%qflx_evap_tot(g)
       l2x_l%rAttr(index_l2x_Fall_swnet,i)  =  clm_l2a%fsa(g)
       if (index_l2x_Fall_fco2_lnd /= 0) then
          l2x_l%rAttr(index_l2x_Fall_fco2_lnd,i) = -clm_l2a%nee(g)  
       end if

       ! Additional fields for DUST, PROGSSLT, dry-deposition and VOC
       ! These are now standard fields, but the check on the index makes sure the driver handles them
       if (index_l2x_Sl_ram1      /= 0 )  l2x_l%rAttr(index_l2x_Sl_ram1,i) = clm_l2a%ram1(g)
       if (index_l2x_Sl_fv        /= 0 )  l2x_l%rAttr(index_l2x_Sl_fv,i)   = clm_l2a%fv(g)
       if (index_l2x_Sl_soilw     /= 0 )  l2x_l%rAttr(index_l2x_Sl_soilw,i)   = clm_l2a%h2osoi_vol(g,1)
       if (index_l2x_Fall_flxdst1 /= 0 )  l2x_l%rAttr(index_l2x_Fall_flxdst1,i)= -clm_l2a%flxdst(g,1)
       if (index_l2x_Fall_flxdst2 /= 0 )  l2x_l%rAttr(index_l2x_Fall_flxdst2,i)= -clm_l2a%flxdst(g,2)
       if (index_l2x_Fall_flxdst3 /= 0 )  l2x_l%rAttr(index_l2x_Fall_flxdst3,i)= -clm_l2a%flxdst(g,3)
       if (index_l2x_Fall_flxdst4 /= 0 )  l2x_l%rAttr(index_l2x_Fall_flxdst4,i)= -clm_l2a%flxdst(g,4)


       ! for dry dep velocities
       if (index_l2x_Sl_ddvel     /= 0 )  then
          l2x_l%rAttr(index_l2x_Sl_ddvel:index_l2x_Sl_ddvel+n_drydep-1,i) = &
               clm_l2a%ddvel(g,:n_drydep)
       end if

       ! for MEGAN VOC emis fluxes
       if (index_l2x_Fall_flxvoc  /= 0 ) then
          l2x_l%rAttr(index_l2x_Fall_flxvoc:index_l2x_Fall_flxvoc+shr_megan_mechcomps_n-1,i) = &
               -clm_l2a%flxvoc(g,:shr_megan_mechcomps_n)
       end if

#ifdef LCH4
       if (index_l2x_Fall_methane /= 0) then
          l2x_l%rAttr(index_l2x_Fall_methane,i) = -clm_l2a%flux_ch4(g) 
       endif
#endif

       ! sign convention is positive downward with 
       ! hierarchy of atm/glc/lnd/rof/ice/ocn.  so water sent from land to rof is positive

       l2x_l%rattr(index_l2x_Flrl_rofliq,i) = clm_l2a%rofliq(g)
       l2x_l%rattr(index_l2x_Flrl_rofice,i) = clm_l2a%rofice(g)

    end do

  end subroutine lnd_export_mct

!====================================================================================

!---------------------------------------------------------------------------
!BOP
!
! !IROUTINE: lnd_import_mct
!
! !INTERFACE:
  subroutine lnd_import_mct( x2l_l, a2l, begg, endg )
!
! !DESCRIPTION:
!
! Convert the input data from the coupler to the land model from MCT import state
! into internal clm data types.
!
!---------------------------------------------------------------------------
! !USES:
    use shr_kind_mod    , only: r8 => shr_kind_r8
    use clm_atmlnd      , only: atm2lnd_type
    use clm_varctl      , only: co2_type, co2_ppmv, iulog, use_c13
    use clm_varctl      , only: startyear_experiment, endyear_experiment, add_temperature
    use clm_varctl      , only: add_co2
    use clm_varcon      , only: rair, o2_molar_const
    use clm_varcon      , only: c13ratio
    use shr_const_mod   , only: SHR_CONST_TKFRZ
    use abortutils      , only: endrun
    use mct_mod         , only: mct_aVect
    use clm_cpl_indices
    use clmtype
    use domainMod        , only : ldomain
    use clm_time_manager,only : get_curr_date, get_nstep, get_curr_calday, get_step_size
    use netcdf

#ifdef CPL_BYPASS
    use shr_const_mod    , only : SHR_CONST_STEBOL
    use shr_kind_mod     , only : r8 => shr_kind_r8, CL => shr_kind_CL
    use fileutils        , only : getavu, relavu
    use spmdmod          , only : masterproc, mpicom, iam, npes, MPI_REAL8, MPI_INTEGER, MPI_STATUS_SIZE
    use ncdio_pio        , only : pio_subsystem
    use shr_pio_mod      , only : shr_pio_getiotype
    use clm_nlUtilsMod   , only : find_nlgroup_name

    use clm_varctl       , only : metdata_type, metdata_bypass, metdata_biases, co2_file, aero_file
    use clm_varctl       , only : const_climate_hist
    use clm_varctl       , only : startdate_add_temperature, startdate_add_co2

    use controlMod       , only : NLFilename

#endif

    implicit none
! !ARGUMENTS:
    type(mct_aVect)   , intent(inout) :: x2l_l   ! Driver MCT import state to land model
    type(atm2lnd_type), intent(inout) :: a2l     ! clm internal input data type
    integer           , intent(in)    :: begg	
    integer           , intent(in)    :: endg
!
! !LOCAL VARIABLES:
    integer  :: g,i,nstep,ier        ! indices, number of steps, and error code
    real(r8) :: forc_rainc           ! rainxy Atm flux mm/s
    real(r8) :: e                    ! vapor pressure (Pa)
    real(r8) :: qsat                 ! saturation specific humidity (kg/kg)
    real(r8) :: forc_rainl           ! rainxy Atm flux mm/s
    real(r8) :: forc_snowc           ! snowfxy Atm flux  mm/s
    real(r8) :: forc_snowl           ! snowfxl Atm flux  mm/s
    real(r8) :: co2_ppmv_diag        ! temporary
    real(r8) :: co2_ppmv_prog        ! temporary
    real(r8) :: co2_ppmv_val         ! temporary
    integer  :: co2_type_idx         ! integer flag for co2_type options
    real(r8) :: esatw                ! saturation vapor pressure over water (Pa)
    real(r8) :: esati                ! saturation vapor pressure over ice (Pa)
    real(r8) :: a0,a1,a2,a3,a4,a5,a6 ! coefficients for esat over water
    real(r8) :: b0,b1,b2,b3,b4,b5,b6 ! coefficients for esat over ice
    real(r8) :: tdc, t               ! Kelvins to Celcius function and its input

    character(len=32), parameter :: sub = 'lnd_import_mct'

    integer ::  yr, mon, day, tod

#ifdef CPL_BYPASS

    integer  :: m,thism              ! indices
    real(r8) :: vp                   ! water vapor pressure (Pa)
    integer  :: thisng, np, num, nu_nml, nml_error
    integer  :: ng_all(100000)
    real(r8) :: swndf, swndr, swvdf, swvdr, ratio_rvrf, frac, q
    real(r8) :: thiscosz, avgcosz, szenith
    integer  :: swrad_period_len, swrad_period_start, thishr, thismin
    real(r8) :: timetemp(2)
    real(r8) :: latixy(500000), longxy(500000)
    integer ::  ierr, varid, dimid, nindex(2), caldaym(13)
    integer ::  ncid, met_ncids(14), mask_ncid, thisncid, ng, tm
    integer ::  aindex(2), tindex(14,2), starti(3), counti(3)
    integer ::  grid_map(500000), zone_map(500000)
    integer ::  met_nvars, nyears_spinup, nyears_trans, starti_site, endi_site
    real(r8) :: smap05_lat(360), smap05_lon(720)
    real(r8) :: smapt62_lat(94), smapt62_lon(192)
    real(r8) :: smap2_lat(96), smap2_lon(144)
    real(r8) :: thisdist, mindist, thislon
    real(r8) :: tbot, tempndep(1,1,158), thiscalday, wt1(14), wt2(14), thisdoy
    real(r8) :: ea
    real(r8) :: site_metdata(14,12)
    real(r8) :: var_month_mean(12)
    !real(r8) :: hdm1(720,360,1), hdm2(720,360,1)
    !real(r8) :: lnfm1(192,94,2920)
    !real(r8) :: ndep1(144,96,1), ndep2(144,96,1)
    !real(r8) :: aerodata(14,144,96,14)
    integer  :: lnfmind(2)
    integer  :: var_month_count(12)
    integer*2 :: temp(1,500000)
    integer :: xtoget, ytoget, thisx, thisy, calday_start
    integer :: sdate_addt, sy_addt, sm_addt, sd_addt
    integer :: sdate_addco2, sy_addco2, sm_addco2, sd_addco2
    character(len=200) metsource_str, thisline

    integer :: av, v, n, nummetdims, g3, gtoget, ztoget, line, mystart, tod_start, thistimelen
    character(len=20) aerovars(14), metvars(14)
    character(len=3) zst
    integer :: stream_year_first_lightng, stream_year_last_lightng, model_year_align_lightng
    integer :: stream_year_first_popdens, stream_year_last_popdens, model_year_align_popdens
    integer :: stream_year_first_ndep,    stream_year_last_ndep,    model_year_align_ndep
    character(len=CL)  :: metdata_fname
    character(len=CL)  :: lightngmapalgo = 'bilinear'! Mapping alogrithm
    character(len=CL)  :: popdensmapalgo = 'bilinear'
    character(len=CL)  :: ndepmapalgo    = 'bilinear'
    character(len=CL)  :: stream_fldFileName_lightng ! lightning stream filename to read
    character(len=CL)  :: stream_fldFileName_popdens ! poplulation density stream filename
    character(len=CL)  :: stream_fldFileName_ndep    ! nitrogen deposition stream filename
    logical :: use_sitedata, has_zonefile, use_daymet, use_livneh
    integer status(MPI_STATUS_SIZE)

    data caldaym / 1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366 /
#endif

    ! Constants to compute vapor pressure
    parameter (a0=6.107799961_r8    , a1=4.436518521e-01_r8, &
               a2=1.428945805e-02_r8, a3=2.650648471e-04_r8, &
               a4=3.031240396e-06_r8, a5=2.034080948e-08_r8, &
               a6=6.136820929e-11_r8)

    parameter (b0=6.109177956_r8    , b1=5.034698970e-01_r8, &
               b2=1.886013408e-02_r8, b3=4.176223716e-04_r8, &
               b4=5.824720280e-06_r8, b5=4.838803174e-08_r8, &
               b6=1.838826904e-10_r8)
    !
    ! function declarations
    !
    tdc(t) = min( 50._r8, max(-50._r8,(t-SHR_CONST_TKFRZ)) )
    esatw(t) = 100._r8*(a0+t*(a1+t*(a2+t*(a3+t*(a4+t*(a5+t*a6))))))
    esati(t) = 100._r8*(b0+t*(b1+t*(b2+t*(b3+t*(b4+t*(b5+t*b6))))))

#ifdef CPL_BYPASS
    namelist /light_streams/         &
        stream_year_first_lightng,  &
        stream_year_last_lightng,   &
        model_year_align_lightng,   &
        lightngmapalgo,             &
        stream_fldFileName_lightng

    namelist /popd_streams/          &
        stream_year_first_popdens,  &
        stream_year_last_popdens,   &
        model_year_align_popdens,   &
        popdensmapalgo,             &
        stream_fldFileName_popdens

    namelist /ndepdyn_nml/        &
        stream_year_first_ndep,  &
        stream_year_last_ndep,   &
        model_year_align_ndep,   &
        ndepmapalgo,             &
        stream_fldFileName_ndep

    stream_fldFileName_lightng = ' '
    stream_fldFileName_popdens = ' '
#endif

!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
! 27 February 2008: Keith Oleson; Forcing height change
!
!EOP
!---------------------------------------------------------------------------


    nstep = get_nstep()	
   
    co2_type_idx = 0
    if (co2_type == 'prognostic') then
       co2_type_idx = 1
    else if (co2_type == 'diagnostic') then
       co2_type_idx = 2
    end if
    if (co2_type == 'prognostic' .and. index_x2l_Sa_co2prog == 0) then
       call endrun( sub//' ERROR: must have nonzero index_x2l_Sa_co2prog for co2_type equal to prognostic' )
    else if (co2_type == 'diagnostic' .and. index_x2l_Sa_co2diag == 0) then
       call endrun( sub//' ERROR: must have nonzero index_x2l_Sa_co2diag for co2_type equal to diagnostic' )
    end if
    
    ! Note that the precipitation fluxes received  from the coupler
    ! are in units of kg/s/m^2. To convert these precipitation rates
    ! in units of mm/sec, one must divide by 1000 kg/m^3 and multiply
    ! by 1000 mm/m resulting in an overall factor of unity.
    ! Below the units are therefore given in mm/s.
    
    call get_curr_date( yr, mon, day, tod )

    do g = begg,endg
        i = 1 + (g - begg)
       
        ! Determine flooding input, sign convention is positive downward and
        ! hierarchy is atm/glc/lnd/rof/ice/ocn.  so water sent from rof to land is negative,
        ! change the sign to indicate addition of water to system.

 	    a2l%forc_flood(g)   = -x2l_l%rattr(index_x2l_Flrr_flood,i)

        a2l%volr(g)   =  x2l_l%rattr(index_x2l_Slrr_volr,i) &
                      * (ldomain%area(g) * 1.e6_r8)

#ifdef CPL_BYPASS
  !-----------------------------------------------------------------------------------------------------

        !read forcing data directly, bypass coupler
        a2l%forc_flood(g) = 0._r8
        a2l%volr(g)       = 0._r8

        !Get meteorological data, concatenated to include whole record
        !Note we only do this at the first timestep and keep the whole forcing dataset in the memory

        !-----------------------------------Meteorological forcing  -----------------------------------

        thiscalday = get_curr_calday()

        !on first timestep, read all the met data for relevant gridcell(s) and store in array.
        !   Met data are held in short integer format to save memory.
        !   Each node must have enough memory to hold these data.
        met_nvars=7
        if (metdata_type(1:3) == 'cpl') met_nvars=14

        if (a2l%loaded_bypassdata == 0) then
          !meteorological forcing
          if (index(metdata_type, 'qian') .gt. 0) then
            a2l%metsource = 0
          else if (index(metdata_type,'cru') .gt. 0) then
            a2l%metsource = 1
          else if (index(metdata_type,'site') .gt. 0) then
            a2l%metsource = 2
          else if (index(metdata_type,'princeton') .gt. 0) then
            a2l%metsource = 3
          else if (index(metdata_type,'gswp3') .gt. 0) then
            a2l%metsource = 4
          else if (index(metdata_type,'cpl') .gt. 0) then
            a2l%metsource = 5
          else
            call endrun( sub//' ERROR: Invalid met data source for cpl_bypass' )
          end if

          use_livneh = .false.
          use_daymet = .false.
          if(index(metdata_type, 'livneh') .gt. 0) then
              use_livneh = .true.
          else if (index(metdata_type, 'daymet') .gt. 0) then
              use_daymet = .true.
          end if

          metvars(1) = 'TBOT'
          metvars(2) = 'PSRF'
          metvars(3) = 'QBOT'
          if (a2l%metsource .eq. 2) metvars(3) = 'RH'
          if (a2l%metsource .ne. 5) metvars(4) = 'FSDS'
          if (a2l%metsource .ne. 5) metvars(5) = 'PRECTmms'
          if (a2l%metsource .ne. 5) metvars(6) = 'WIND'
          metvars(4) = 'FSDS'
          metvars(5) = 'PRECTmms'
          metvars(6) = 'WIND'
          metvars(7) = 'FLDS'
          if (a2l%metsource .eq. 5) then
              metvars(4) = 'SWNDF'
              metvars(5) = 'RAINC'
              metvars(6) = 'U'
              metvars(8) = 'SWNDR'
              metvars(9) = 'SWVDF'
              metvars(10) = 'SWVDR'
              metvars(11) = 'RAINL'
              metvars(12) = 'SNOWC'
              metvars(13) = 'SNOWL'
              metvars(14) = 'V'
          else
              metvars(4) = 'FSDS'
              metvars(5) = 'PRECTmms'
              metvars(6) = 'WIND'
          end if

          !set defaults
          a2l%startyear_met       = 1901
          a2l%endyear_met_spinup  = 1920
          if (a2l%metsource == 0) then
            metsource_str = 'qian'
            a2l%startyear_met       = 1948
            a2l%endyear_met_spinup  = 1972
            a2l%endyear_met_trans   = 2004
          else if (a2l%metsource == 1) then
            metsource_str = 'cruncep'
            a2l%endyear_met_trans  = 2016
          else if (a2l%metsource == 2) then
            metsource_str = 'site'
            !get year information from file
            ierr = nf90_open(trim(metdata_bypass) // '/all_hourly.nc', nf90_nowrite, ncid)
            ierr = nf90_inq_varid(ncid, 'start_year', varid)
            ierr = nf90_get_var(ncid, varid, a2l%startyear_met)
            ierr = nf90_inq_varid(ncid, 'end_year', varid)
            ierr = nf90_get_var(ncid, varid, a2l%endyear_met_spinup)
            ierr = nf90_close(ncid)
            a2l%endyear_met_trans = a2l%endyear_met_spinup
          else if (a2l%metsource == 3) then
            metsource_str = 'princeton'
            a2l%endyear_met_trans = 2012
          else if (a2l%metsource == 4) then
            a2l%endyear_met_trans  = 2014
            if(index(metdata_type, 'v1') .gt. 0) a2l%endyear_met_trans  = 2010
          else if (a2l%metsource == 5) then
            a2l%startyear_met      = 566 !76
            a2l%endyear_met_spinup = 590 !100
            a2l%endyear_met_trans  = 590 !100
          end if

          if (use_livneh) then
              a2l%startyear_met      = 1950
              a2l%endyear_met_spinup = 1969
          else if (use_daymet) then
              a2l%startyear_met      = 1980
              a2l%endyear_met_spinup = a2l%endyear_met_trans
          end if

          nyears_spinup = a2l%endyear_met_spinup - &
                             a2l%startyear_met + 1
          nyears_trans  = a2l%endyear_met_trans - &
                             a2l%startyear_met  + 1

          !check for site data in run directory (monthly mean T, precip)
          inquire(file=trim(metdata_biases), exist=use_sitedata)

          !get grid lat/lon information, zone mappings
          inquire(file=trim(metdata_bypass) // '/zone_mappings.txt', exist=has_zonefile)
          if (has_zonefile) then
            open(unit=13, file=trim(metdata_bypass) // '/zone_mappings.txt')
          else if (a2l%metsource .ne. 2) then
            call endrun( sub//' ERROR: Zone mapping file does not exist for cpl_bypass' )
          end if

          if (a2l%metsource .ne. 2) then
            ng = 0     !number of points
            do v=1,500000
              read(13,*, end=10), longxy(v), latixy(v), zone_map(v), grid_map(v)
              ng = ng + 1
            end do
10          continue
            close(unit=13)

            !Figure out the closest point and which zone file to open
            mindist=99999
            do g3 = 1,ng
              ! in CPL_BYPASS met dataset, longitude is in format of 0-360, but 'ldomain%lonc(g)' may or may not.
              if (ldomain%lonc(g) .lt. 0) then
                if (longxy(g3) >= 180) longxy(g3) = longxy(g3)-360._r8
              else if (ldomain%lonc(g) .ge. 180) then
                if (longxy(g3) < 0) longxy(g3) = longxy(g3) + 360._r8
              end if
              thisdist = 100*((latixy(g3) - ldomain%latc(g))**2 + &
                              (longxy(g3) - ldomain%lonc(g))**2)**0.5
              if (thisdist .lt. mindist) then
                mindist = thisdist
                ztoget = zone_map(g3)
                gtoget = grid_map(g3)
              end if
            end do
          else
            gtoget = 1
          end if

          !get the site metdata for bias correction if they exist (lat/lons must match domain file)
          if (use_sitedata) then
            open(unit=9, file=trim(metdata_biases),status='old')
            read(9,*) thisline
            site_metdata(:,:)=-999._r8
            do while ((site_metdata(1,1) .lt. ldomain%lonc(g) - 0.01 .or. &
                    site_metdata(1,1) .gt. ldomain%lonc(g) + 0.01) .and. &
                      (site_metdata(2,1) .lt. ldomain%latc(g) - 0.01 .or. &
                       site_metdata(2,1) .gt. ldomain%latc(g) + 0.01))
              read(9,*) site_metdata(1:7,1)
              if (site_metdata(1,1) .lt. 0) site_metdata(1,1) = site_metdata(1,1)+360._r8
            end do
            do line=2,12
              read(9,*) site_metdata(1:7,line)
            end do
            close(unit=9)
          end if

          do v=1,met_nvars
            write(zst, '(I3)') 100+ztoget
            if (a2l%metsource == 0) then
                metdata_fname =  trim(metsource_str) // '_' // trim(metvars(v)) // '_z' // zst(2:3) // '.nc'
            else if (a2l%metsource == 1) then
                metdata_fname = 'CRUNCEP.v5_' // trim(metvars(v)) // '_1901-2013_z' // zst(2:3) // '.nc'
                if (use_livneh .and. ztoget .ge. 16 .and. ztoget .le. 20) then
                    metdata_fname = 'CRUNCEP5_Livneh_' // trim(metvars(v)) // '_1950-2013_z' // zst(2:3) // '.nc'
                else if (use_daymet .and. ztoget .ge. 16 .and. ztoget .le. 20) then
                    metdata_fname = 'CRUNCEP5_Daymet3_' // trim(metvars(v)) // '_1980-2013_z' // zst(2:3) // '.nc'
                end if
            else if (a2l%metsource == 2) then
                metdata_fname = 'all_hourly.nc'
            else if (a2l%metsource == 3) then
               metdata_fname = 'Princeton_' // trim(metvars(v)) // '_1901-2012_z' // zst(2:3) // '.nc'
                if (use_livneh .and. ztoget .ge. 16 .and. ztoget .le. 20) then
                    metdata_fname = 'Princeton_Livneh_' // trim(metvars(v)) // '_1950-2012_z' // zst(2:3) // '.nc'
                else if (use_daymet .and. ztoget .ge. 16 .and. ztoget .le. 20) then
                    metdata_fname = 'Princeton_Daymet3_' // trim(metvars(v)) // '_1980-2012_z' // zst(2:3) // '.nc'
                end if
            else if (a2l%metsource == 4) then
                metdata_fname = 'GSWP3_' // trim(metvars(v)) // '_1901-2014_z' // zst(2:3) // '.nc'
                if(index(metdata_type, 'v1') .gt. 0) &
                    metdata_fname = 'GSWP3_' // trim(metvars(v)) // '_1901-2010_z' // zst(2:3) // '.nc'

                if (use_livneh .and. ztoget .ge. 16 .and. ztoget .le. 20) then
                    metdata_fname = 'GSWP3_Livneh_' // trim(metvars(v)) // '_1950-2010_z' // zst(2:3) // '.nc'
                else if (use_daymet .and. ztoget .ge. 16 .and. ztoget .le. 20) then
                    metdata_fname = 'GSWP3_Daymet3_' // trim(metvars(v)) // '_1980-2010_z' // zst(2:3) // '.nc'
                end if
            else if (a2l%metsource == 5) then
                    !metdata_fname = 'WCYCL1850S.ne30_' // trim(metvars(v)) // '_0076-0100_z' // zst(2:3) // '.nc'
                    metdata_fname = 'CBGC1850S.ne30_' // trim(metvars(v)) // '_0566-0590_z' // zst(2:3) // '.nc'
            end if

            ierr = nf90_open(trim(metdata_bypass) // '/' // trim(metdata_fname), NF90_NOWRITE, met_ncids(v))
            if (ierr .ne. 0) call endrun(msg=' ERROR: Failed to open cpl_bypass input meteorology file' )

            !get timestep information
            ierr = nf90_inq_dimid(met_ncids(v), 'DTIME', dimid)
            ierr = nf90_Inquire_Dimension(met_ncids(v), dimid, len = a2l%timelen(v))

            starti(1) = 1
            counti(1) = 2
            ierr = nf90_inq_varid(met_ncids(v), 'DTIME', varid)
            ierr = nf90_get_var(met_ncids(v), varid, timetemp, starti(1:1), counti(1:1))
            a2l%timeres(v)        = (timetemp(2)-timetemp(1))*24._r8
            a2l%npf(v)            = 86400d0*(timetemp(2)-timetemp(1))/get_step_size()
            a2l%timelen_spinup(v) = nyears_spinup*(365*nint(24./a2l%timeres(v)))

            ierr = nf90_inq_varid(met_ncids(v), trim(metvars(v)), varid)

            !get the conversion factors
            ierr = nf90_get_att(met_ncids(v), varid, 'scale_factor', a2l%scale_factors(v))
            if (ierr .ne. 0) a2l%scale_factors(v) = 1.0d0

            ierr = nf90_get_att(met_ncids(v), varid, 'add_offset', a2l%add_offsets(v))
            if (ierr .ne. 0) a2l%add_offsets(v) = 0.0d0

            !get the met data
            starti(1) = 1
            starti(2) = gtoget
            counti(1) = a2l%timelen_spinup(v)
            counti(2) = 1
            if (.not. const_climate_hist .and. (yr .ge. 1850 .or. use_sitedata)) counti(1) = a2l%timelen(v)

            if (i == 1 .and. v == 1)  then
              allocate(a2l%atm_input       (met_nvars,begg:endg,1,1:counti(1)))
            end if

            ierr = nf90_get_var(met_ncids(v), varid, a2l%atm_input(v,g:g,1,1:counti(1)), starti(1:2), counti(1:2))
            ierr = nf90_close(met_ncids(v))

            if (use_sitedata .and. v == 1) then
                starti_site = max((nint(site_metdata(4,1))-a2l%startyear_met) * &
                                     365*nint(24./a2l%timeres(v))+1,1)
                endi_site   = (min(a2l%endyear_met_trans,nint(site_metdata(5,1))) - &
                                     a2l%startyear_met+1)*(365*nint(24./a2l%timeres(v)))
            end if

            a2l%var_offset(v,g,:) = 0._r8
            a2l%var_mult(v,g,:)   = 1._r8

            if (use_sitedata) then
              !Compute monthly biases for site vs. reanalysis
              var_month_mean(:)  = 0._r8
              var_month_count(:) = 0
              do i=starti_site, endi_site
                thisdoy = mod(i,365*nint(24./a2l%timeres(v)))/(nint(24./a2l%timeres(v)))+1
                do m=1,12
                  if (thisdoy .ge. caldaym(m) .and. thisdoy .lt. caldaym(m+1)) thism = m
                enddo
                var_month_mean(thism) = var_month_mean(thism) + (a2l%atm_input(v,g,1,i)* &
                                          a2l%scale_factors(v) + a2l%add_offsets(v))
                var_month_count(thism) = var_month_count(thism)+1
              end do

              do m = 1,12
                var_month_mean(m) = var_month_mean(m)/var_month_count(m)
                !calculate offset and linear bias factors for temperature and precipitation
                if (v .eq. 1) a2l%var_offset(v,g,m) = (site_metdata(6,m)+SHR_CONST_TKFRZ) - var_month_mean(m)
                if (v .eq. 5 .and. var_month_mean(m) .gt. 0) &
                      a2l%var_mult(v,g,m) = (site_metdata(7,m))/(caldaym(m+1)-caldaym(m))/24._r8/ &
                                                      3600._r8 / var_month_mean(m)
              end do
            end if

            !Align spinups and transient simulations
            !figure out which year to start with (assuming spinups always use integer multiple of met cycles)
            mystart = a2l%startyear_met
            do while (mystart > 1850)
              mystart = mystart - nyears_spinup
            end do
            if (a2l%metsource == 5) mystart=1850

            if (yr .lt. 1850) then
              !a2l%tindex(g,v,1) = (mod(yr-1,nyears_spinup) + (1850-mystart)) * 365 * nint(24./a2l%timeres(v))
              a2l%tindex(g,v,1) = mod(yr+1849-mystart,nyears_spinup) * 365 * nint(24./a2l%timeres(v))
            else if (yr .le. a2l%endyear_met_spinup) then
              !a2l%tindex(g,v,1) = (mod(yr-1850,nyears_spinup) + (1850-mystart)) * 365 * nint(24./a2l%timeres(v))
              a2l%tindex(g,v,1) = mod(yr-mystart,nyears_spinup) * 365 * nint(24./a2l%timeres(v))
            else
              a2l%tindex(g,v,1) = (yr - a2l%startyear_met) * 365 * nint(24./a2l%timeres(v))
            end if
            !adjust for starts not at beginning of year (but currently MUST begin at hour 0)
            a2l%tindex(g,v,1) = a2l%tindex(g,v,1) + (caldaym(mon)+day-2)* &
                                         nint(24./a2l%timeres(v))

            a2l%tindex(g,v,2) = a2l%tindex(g,v,1) + 1
            if (a2l%tindex(g,v,1) == 0) then
              a2l%tindex(g,v,1) = a2l%timelen(v)
              if (yr .le. a2l%endyear_met_spinup) a2l%tindex(g,v,1) = a2l%timelen_spinup(v)
             end if
          end do    !end variable loop
        else
          do v=1,met_nvars
            if (a2l%npf(v) - 1._r8 .gt. 1e-3) then
              if (v .eq. 4 .or. v .eq. 5 .or. (v .ge. 8 .and. v .le. 13)) then    !rad/Precipitation
                if (mod(tod/get_step_size(),nint(a2l%npf(v))) == 1 .and. nstep .gt. 3) then
                  a2l%tindex(g,v,1) = a2l%tindex(g,v,1)+1
                  a2l%tindex(g,v,2) = a2l%tindex(g,v,2)+1
                end if
              else
                if (mod(tod/get_step_size()-1,nint(a2l%npf(v))) <= a2l%npf(v)/2._r8 .and. &
                    mod(tod/get_step_size(),nint(a2l%npf(v))) > a2l%npf(v)/2._r8) then
                  a2l%tindex(g,v,1) = a2l%tindex(g,v,1)+1
                  a2l%tindex(g,v,2) = a2l%tindex(g,v,2)+1
                end if
              end if
            else
              a2l%tindex(g,v,1) = a2l%tindex(g,v,1)+nint(1/a2l%npf(v))
              a2l%tindex(g,v,2) = a2l%tindex(g,v,2)+nint(1/a2l%npf(v))
            end if

            if (const_climate_hist .or. yr .le. a2l%startyear_met) then
              if (a2l%tindex(g,v,1) .gt. a2l%timelen_spinup(v)) a2l%tindex(g,v,1) = 1
              if (a2l%tindex(g,v,2) .gt. a2l%timelen_spinup(v)) a2l%tindex(g,v,2) = 1
            else if (yr .gt. a2l%endyear_met_trans) then
              if (a2l%tindex(g,v,1) .gt. a2l%timelen(v)) then
                 a2l%tindex(g,v,1) = a2l%timelen(v)-a2l%timelen_spinup(v)+1
              end if
              if (a2l%tindex(g,v,2) .gt. a2l%timelen(v)) then
                 a2l%tindex(g,v,2) = a2l%timelen(v)-a2l%timelen_spinup(v)+1
              end if
            end if

            !if (yr .gt. a2l%startyear_met) then
            !  if (a2l%tindex(g,v,1) .gt. a2l%timelen(v)) a2l%tindex(g,v,1) = 1
            !  if (a2l%tindex(g,v,2) .gt. a2l%timelen(v)) a2l%tindex(g,v,2) = 1
            !else
            !  if (a2l%tindex(g,v,1) .gt. a2l%timelen_spinup(v)) a2l%tindex(g,v,1) = 1
            !  if (a2l%tindex(g,v,2) .gt. a2l%timelen_spinup(v)) a2l%tindex(g,v,2) = 1
            !end if
          end do
        end if

        tindex = a2l%tindex(g,:,:)

        !get weights for linear interpolation
        do v=1,met_nvars
          if (a2l%npf(v) - 1._r8 .gt. 1e-3) then
               wt1(v) = 1._r8 - (mod((tod+86400)/get_step_size()-a2l%npf(v)/2._r8, &
                   a2l%npf(v))*1._r8)/a2l%npf(v)
               wt2(v) = 1._r8 - wt1(v)
          else
             wt1(v) = 0._r8
             wt2(v) = 1._r8
          end if
        end do

        !Air temperature
        a2l%forc_t(g)  = min(((a2l%atm_input(1,g,1,tindex(1,1))*a2l%scale_factors(1)+ &
                               a2l%add_offsets(1))*wt1(1) + (a2l%atm_input(1,g,1,tindex(1,2))* &
                               a2l%scale_factors(1)+a2l%add_offsets(1))*wt2(1)) * &
                                              a2l%var_mult(1,g,mon) + a2l%var_offset(1,g,mon), 323._r8)
        a2l%forc_th(g) = min(((a2l%atm_input(1,g,1,tindex(1,1))*a2l%scale_factors(1)+ &
                               a2l%add_offsets(1))*wt1(1) + (a2l%atm_input(1,g,1,tindex(1,2))* &
                               a2l%scale_factors(1)+a2l%add_offsets(1))*wt2(1)) * &
                               a2l%var_mult(1,g,mon) + a2l%var_offset(1,g,mon), 323._r8)

        tbot = a2l%forc_t(g)

        !Air pressure
        a2l%forc_pbot(g) = max(((a2l%atm_input(2,g,1,tindex(2,1))*a2l%scale_factors(2)+ &
                                 a2l%add_offsets(2))*wt1(2) + (a2l%atm_input(2,g,1,tindex(2,2)) &
                                 *a2l%scale_factors(2)+a2l%add_offsets(2))*wt2(2)) * &
                                 a2l%var_mult(2,g,mon) + a2l%var_offset(2,g,mon), 4e4_r8)
        !Specific humidity
        a2l%forc_q(g) = max(((a2l%atm_input(3,g,1,tindex(3,1))*a2l%scale_factors(3)+ &
                                  a2l%add_offsets(3))*wt1(3) + (a2l%atm_input(3,g,1,tindex(3,2)) &
                                  *a2l%scale_factors(3)+a2l%add_offsets(3))*wt2(3)) * &
                                  a2l%var_mult(3,g,mon) + a2l%var_offset(3,g,mon), 1e-9_r8)
        !
        if (tbot > SHR_CONST_TKFRZ) then
          e = esatw(tdc(tbot))
        else
          e = esati(tdc(tbot))
        end if
        qsat = 0.622_r8*e / (a2l%forc_pbot(g) - 0.378_r8*e)
        if (a2l%metsource == 2) then                            !convert RH to qbot, when input is actually RH
          a2l%forc_q(g) = qsat * a2l%forc_q(g) / 100.0_r8
        else if(a2l%forc_q(g)>qsat) then     ! data checking for specific humidity
          a2l%forc_q(g) = qsat
        end if

        !use longwave from file if provided
        a2l%forc_lwrad(g) = ((a2l%atm_input(7,g,1,tindex(7,1))*a2l%scale_factors(7)+ &
                                                        a2l%add_offsets(7))*wt1(7) + (a2l%atm_input(7,g,1,tindex(7,2)) &
                                                        *a2l%scale_factors(7)+a2l%add_offsets(7))*wt2(7)) * &
                                                        a2l%var_mult(7,g,mon) + a2l%var_offset(7,g,mon)
        if (a2l%forc_lwrad(g) .le. 50 .or. a2l%forc_lwrad(g) .ge. 600) then
        !Longwave radiation (calculated from air temperature, humidity)
            e =  a2l%forc_pbot(g) * a2l%forc_q(g) / &
                 (0.622_R8 + 0.378_R8 * a2l%forc_q(g) )
            ea = 0.70_R8 + 5.95e-05_R8 * 0.01_R8 * e * exp(1500.0_R8/tbot)
            a2l%forc_lwrad(g) = ea * SHR_CONST_STEBOL * tbot**4
        end if

        !Shortwave radiation (cosine zenith angle interpolation)
        thishr = (tod-get_step_size()/2)/3600
        if (thishr < 0) thishr=thishr+24
        thismin = mod((tod-get_step_size()/2)/60, 60)
        thiscosz = max(cos(szenith(ldomain%lonc(g),ldomain%latc(g),0,int(thiscalday),thishr,thismin,0)* &
                        3.14159265358979/180.0d0), 0.001d0)
        avgcosz = 0d0
        if (a2l%npf(4) - 1._r8 .gt. 1e-3) then
          swrad_period_len   = get_step_size()*nint(a2l%npf(4))
          swrad_period_start = ((tod-get_step_size()/2)/swrad_period_len) * swrad_period_len
          !set to last period if first model timestep of the day
          if (tod-get_step_size()/2 < 0) swrad_period_start = ((86400-get_step_size()/2)/swrad_period_len) * swrad_period_len

          do tm=1,nint(a2l%npf(4))
            !Get the average cosine zenith angle over the time resolution of the input data
            thishr  = (swrad_period_start+(tm-1)*get_step_size()+get_step_size()/2)/3600
            if (thishr > 23) thishr=thishr-24
            thismin = mod((swrad_period_start+(tm-1)*get_step_size()+get_step_size()/2)/60, 60)
            avgcosz  = avgcosz + max(cos(szenith(ldomain%lonc(g),ldomain%latc(g),0,int(thiscalday),thishr, thismin, 0) &
                       *3.14159265358979/180.0d0), 0.001d0)/a2l%npf(4)
          end do
        else
          avgcosz = thiscosz
        end if
        if (thiscosz > 0.001d0) then
          wt2(4) = min(thiscosz/avgcosz, 10.0_r8)
        else
          wt2(4) = 0d0
        end if

        if (a2l%metsource == 5) then
            wt2(4)=1.0   !cosz interp not working
            wt2(8:10)=1.0
            swndf = max(((a2l%atm_input(4,g,1,tindex(4,2))*a2l%scale_factors(4)+ &
                                     a2l%add_offsets(4))*wt2(4)), 0.0_r8)
            swndr = max(((a2l%atm_input(8,g,1,tindex(8,2))*a2l%scale_factors(8)+ &
                                     a2l%add_offsets(8))*wt2(8)), 0.0_r8)
            swvdf = max(((a2l%atm_input(9,g,1,tindex(9,2))*a2l%scale_factors(9)+ &
                                     a2l%add_offsets(9))*wt2(9)), 0.0_r8)
            swvdr = max(((a2l%atm_input(10,g,1,tindex(10,2))*a2l%scale_factors(10)+ &
                                     a2l%add_offsets(10))*wt2(10)), 0.0_r8)
            a2l%forc_solad(g,2) = swndr
            a2l%forc_solad(g,1) = swvdr
            a2l%forc_solai(g,2) = swndf
            a2l%forc_solai(g,1) = swvdf
        else
            swndr = max(((a2l%atm_input(4,g,1,tindex(4,2))*a2l%scale_factors(4)+ &
                                     a2l%add_offsets(4))*wt2(4)) * 0.50_R8, 0.0_r8)
            swndf = max(((a2l%atm_input(4,g,1,tindex(4,2))*a2l%scale_factors(4)+ &
                                    a2l%add_offsets(4))*wt2(4))*0.50_R8, 0.0_r8)
            swvdr = max(((a2l%atm_input(4,g,1,tindex(4,2))*a2l%scale_factors(4)+ &
                                    a2l%add_offsets(4))*wt2(4))*0.50_R8, 0.0_r8)
            swvdf = max(((a2l%atm_input(4,g,1,tindex(4,2))*a2l%scale_factors(4)+ &
                                    a2l%add_offsets(4))*wt2(4))*0.50_R8, 0.0_r8)
            ratio_rvrf =   min(0.99_R8,max(0.29548_R8 + 0.00504_R8*swndr &
                           -1.4957e-05_R8*swndr**2 + 1.4881e-08_R8*swndr**3,0.01_R8))
            a2l%forc_solad(g,2) = ratio_rvrf*swndr
            a2l%forc_solai(g,2) = (1._R8 - ratio_rvrf)*swndf
            ratio_rvrf =   min(0.99_R8,max(0.17639_R8 + 0.00380_R8*swvdr  &
                               -9.0039e-06_R8*swvdr**2 +8.1351e-09_R8*swvdr**3,0.01_R8))
            a2l%forc_solad(g,1) = ratio_rvrf*swvdr
            a2l%forc_solai(g,1) = (1._R8 - ratio_rvrf)*swvdf
        end if
        !Rain and snow
        if (a2l%metsource == 5) then
          forc_rainc = max((((a2l%atm_input(5,g,1,tindex(5,2))*a2l%scale_factors(5)+ &
                                        a2l%add_offsets(5)))*a2l%var_mult(5,g,mon) + &
                                        a2l%var_offset(5,g,mon)), 0.0_r8)
          forc_rainl = max((((a2l%atm_input(11,g,1,tindex(11,2))*a2l%scale_factors(11)+ &
                                        a2l%add_offsets(11)))*a2l%var_mult(11,g,mon) + &
                                        a2l%var_offset(11,g,mon)), 0.0_r8)
          forc_snowc = max((((a2l%atm_input(12,g,1,tindex(12,2))*a2l%scale_factors(12)+ &
                                        a2l%add_offsets(12)))*a2l%var_mult(12,g,mon) + &
                                        a2l%var_offset(12,g,mon)), 0.0_r8)
          forc_snowl = max((((a2l%atm_input(13,g,1,tindex(13,2))*a2l%scale_factors(13)+ &
                                        a2l%add_offsets(13)))*a2l%var_mult(13,g,mon) + &
                                          a2l%var_offset(13,g,mon)), 0.0_r8)
        else
          frac = (a2l%forc_t(g) - SHR_CONST_TKFRZ)*0.5_R8       ! ramp near freezing
          frac = min(1.0_R8,max(0.0_R8,frac))           ! bound in [0,1]
          !Don't interpolate rainfall data
          forc_rainc = 0.1_R8 * frac * max((((a2l%atm_input(5,g,1,tindex(5,2))*a2l%scale_factors(5)+ &
                                        a2l%add_offsets(5)))*a2l%var_mult(5,g,mon) + &
                                        a2l%var_offset(5,g,mon)), 0.0_r8)
          forc_rainl = 0.9_R8 * frac * max((((a2l%atm_input(5,g,1,tindex(5,2))*a2l%scale_factors(5)+ &
                                         a2l%add_offsets(5)))*a2l%var_mult(5,g,mon) + &
                                         a2l%var_offset(5,g,mon)), 0.0_r8)
          forc_snowc = 0.1_R8 * (1.0_R8 - frac) * max((((a2l%atm_input(5,g,1,tindex(5,2))*a2l%scale_factors(5)+ &
                  a2l%add_offsets(5)))*a2l%var_mult(5,g,mon) + a2l%var_offset(5,g,mon)), 0.0_r8)
          forc_snowl = 0.9_R8 * (1.0_R8 - frac) * max((((a2l%atm_input(5,g,1,tindex(5,2))*a2l%scale_factors(5)+ &
                  a2l%add_offsets(5))) * a2l%var_mult(5,g,mon) + a2l%var_offset(5,g,mon)), 0.0_r8)
        end if
        !Wind
        a2l%forc_u(g) = (a2l%atm_input(6,g,1,tindex(6,1))*a2l%scale_factors(6)+ &
                                     a2l%add_offsets(6))*wt1(6) + (a2l%atm_input(6,g,1,tindex(6,2))* &
                                     a2l%scale_factors(6)+a2l%add_offsets(6))*wt2(6)
        if (a2l%metsource == 5) then
          a2l%forc_v(g) = (a2l%atm_input(14,g,1,tindex(14,1))*a2l%scale_factors(14)+ &
                                     a2l%add_offsets(14))*wt1(14) + (a2l%atm_input(14,g,1,tindex(14,2))* &
                                     a2l%scale_factors(14)+a2l%add_offsets(14))*wt2(14)
        else
            a2l%forc_v(g) = 0.0_R8
        end if
        a2l%forc_hgt(g) = 30.0_R8 !(a2l%atm_input(8,g,1,tindex(1))*wt1 + &
                                             !a2l%atm_input(8,g,1,tindex(2))*wt2)    ! zgcmxy  Atm state, default=30m

  !------------------------------------Fire data -------------------------------------------------------

        nindex(1) = yr-1848
        nindex(2) = nindex(1)+1
        if (yr .lt. 1850 .or. const_climate_hist) nindex(1:2) = 2
        if (yr .ge. 2010 .and. .not. const_climate_hist) nindex(1:2) = 161

#ifndef NOFIRE
        !if (use_cn .and. .not.use_nofire) then
          if (a2l%loaded_bypassdata == 0 .or. (mon .eq. 1 .and. day .eq. 1 .and. tod .eq. 0)) then
            if (masterproc .and. i .eq. 1) then
              ! Read pop_dens streams namelist to get filename
              nu_nml = getavu()
              open(nu_nml, file=trim(NLFilename), status='old', iostat=nml_error )
              call find_nlgroup_name(nu_nml, 'popd_streams', status=nml_error)
              if (nml_error == 0) then
                  read(nu_nml, nml=popd_streams,iostat=nml_error)
                  if (nml_error /= 0) then
                      call endrun(msg='ERROR reading popdens namelist')
                  end if
              end if
              close(nu_nml)
              call relavu( nu_nml )

              ierr = nf90_open(trim(stream_fldFileName_popdens), NF90_NOWRITE, ncid)
              ierr = nf90_inq_varid(ncid, 'lat', varid)
              ierr = nf90_get_var(ncid, varid, smap05_lat)
              ierr = nf90_inq_varid(ncid, 'lon', varid)
              ierr = nf90_get_var(ncid, varid, smap05_lon)
              ierr = nf90_inq_varid(ncid, 'hdm', varid)
              starti(1:2) = 1
              starti(3)   = nindex(1)
              counti(1) = 720
              counti(2) = 360
              counti(3) = 1
              ierr = nf90_get_var(ncid, varid, a2l%hdm1, starti, counti)
              starti(3) = nindex(2)
              if (nindex(1) .ne. nindex(2)) then
                  ierr = nf90_get_var(ncid, varid, a2l%hdm2, starti, counti)
              else
                  a2l%hdm2 = a2l%hdm1
              end if
              ierr = nf90_close(ncid)
            end if

            if (i .eq. 1) then
              call mpi_bcast (a2l%hdm1, 360*720, MPI_REAL8, 0, mpicom, ier)
              call mpi_bcast (a2l%hdm2, 360*720, MPI_REAL8, 0, mpicom, ier)
              call mpi_bcast (smap05_lon, 720, MPI_REAL8, 0, mpicom, ier)
              call mpi_bcast (smap05_lat, 360, MPI_REAL8, 0, mpicom, ier)
            end if
          end if

          !figure out which point to get
          if (a2l%loaded_bypassdata == 0) then
            mindist=99999
            do thisx = 1,720
              do thisy = 1,360
                  if (ldomain%lonc(g) .lt. 0) then
                      if (smap05_lon(thisx) >= 180) smap05_lon(thisx) = smap05_lon(thisx)-360._r8
                  else if (ldomain%lonc(g) .ge. 180) then
                      if (smap05_lon(thisx) < 0) smap05_lon(thisx) = smap05_lon(thisx) + 360._r8
                  end if
                  thisdist = 100*((smap05_lat(thisy) - ldomain%latc(g))**2 + &
                          (smap05_lon(thisx) - ldomain%lonc(g))**2)**0.5
                  if (thisdist .lt. mindist) then
                      mindist = thisdist
                      a2l%hdmind(g,1) = thisx
                      a2l%hdmind(g,2) = thisy
                  end if
              end do
            end do
          end if
          !get weights for interpolation
          wt1(1) = 1._r8 - (thiscalday -1._r8)/365._r8
          wt2(1) = 1._r8 - wt1(1)
          a2l%forc_hdm(g) = a2l%hdm1(a2l%hdmind(g,1),a2l%hdmind(g,2),1)*wt1(1) + &
                                     a2l%hdm2(a2l%hdmind(g,1),a2l%hdmind(g,2),1)*wt2(1)

          if (a2l%loaded_bypassdata .eq. 0 .and. masterproc .and. i .eq. 1) then
            ! Read light_streams namelist to get filename
            nu_nml = getavu()
            open( nu_nml, file=trim(NLFilename), status='old', iostat=nml_error )
            call find_nlgroup_name(nu_nml, 'light_streams', status=nml_error)
            if (nml_error == 0) then
              read(nu_nml, nml=light_streams,iostat=nml_error)
              if (nml_error /= 0) then
                call endrun(msg='ERROR reading light_streams namelist')
              end if
            end if
            close(nu_nml)
            call relavu( nu_nml )

            !Get all of the data (master processor only)
            allocate(a2l%lnfm_all       (192,94,2920))
            ierr = nf90_open(trim(stream_fldFileName_lightng), NF90_NOWRITE, ncid)
            ierr = nf90_inq_varid(ncid, 'lat', varid)
            ierr = nf90_get_var(ncid, varid, smapt62_lat)
            ierr = nf90_inq_varid(ncid, 'lon', varid)
            ierr = nf90_get_var(ncid, varid, smapt62_lon)
            ierr = nf90_inq_varid(ncid, 'lnfm', varid)
            ierr = nf90_get_var(ncid, varid, a2l%lnfm_all)
            ierr = nf90_close(ncid)
          end if
          if (a2l%loaded_bypassdata .eq. 0 .and. i .eq. 1) then
            call mpi_bcast (smapt62_lon, 192, MPI_REAL8, 0, mpicom, ier)
            call mpi_bcast (smapt62_lat, 94, MPI_REAL8, 0, mpicom, ier)
          end if
          if (a2l%loaded_bypassdata .eq. 0) then
            mindist=99999
            do thisx = 1,192
              do thisy = 1,94
                if (ldomain%lonc(g) .lt. 0) then
                  if (smapt62_lon(thisx) >= 180) smapt62_lon(thisx) = smapt62_lon(thisx)-360._r8
                else if (ldomain%lonc(g) .ge. 180) then
                  if (smapt62_lon(thisx) < 0) smapt62_lon(thisx) = smapt62_lon(thisx) + 360._r8
                end if
                thisdist = 100*((smapt62_lat(thisy) - ldomain%latc(g))**2 + &
                            (smapt62_lon(thisx) - ldomain%lonc(g))**2)**0.5
                if (thisdist .lt. mindist) then
                  mindist = thisdist
                  lnfmind(1) = thisx
                  lnfmind(2) = thisy
                end if
              end do
            end do
            if (masterproc) then
              a2l%lnfm(g,:) = a2l%lnfm_all(lnfmind(1),lnfmind(2),:)
              do np = 1,npes-1
                if (i == 1) then
                  call mpi_recv(thisng,  1, MPI_INTEGER, np, 100000+np, mpicom, status, ier)
                  ng_all(np) = thisng
                end if
                if (i <= ng_all(np)) then
                  call mpi_recv(lnfmind, 2, MPI_INTEGER, np, 200000+np, mpicom, status, ier)
                  call mpi_send(a2l%lnfm_all(lnfmind(1),lnfmind(2),:), 2920, &
                            MPI_REAL8, np, 300000+np, mpicom, ier)
                end if
              end do
            else
              if (i == 1)  call mpi_send(thisng,  1, MPI_INTEGER, 0, 100000+iam, mpicom, ier)
              call mpi_send(lnfmind, 2, MPI_INTEGER, 0, 200000+iam, mpicom, ier)
              call mpi_recv(a2l%lnfm(g,:), 2920, MPI_REAL8, 0, 300000+iam, mpicom, status, ier)
            end if

            !Lightning data is 3-hourly.  Does not currently interpolate.
            a2l%forc_lnfm(g) = a2l%lnfm(g, ((int(thiscalday)-1)*8+tod/(3600*3))+1)

          !end if ! end of 'if use_cn .and. .not.use_nofire
#endif
          !------------------------------------Nitrogen deposition----------------------------------------------

          !DMR note - ndep will NOT be correct if more than 1850 years of model
          !spinup (model year > 1850)
          nindex(1) = min(max(yr-1848,2), 168)
          nindex(2) = min(nindex(1)+1, 168)

          if (a2l%loaded_bypassdata .eq. 0 .or. (mon .eq. 1 .and. day .eq. 1 .and. tod .eq. 0)) then
            if (masterproc .and. i .eq. 1) then
              nu_nml = getavu()
              open( nu_nml, file=trim(NLFilename), status='old', iostat=nml_error )
              call find_nlgroup_name(nu_nml, 'ndepdyn_nml', status=nml_error)
              if (nml_error == 0) then
                read(nu_nml, nml=ndepdyn_nml,iostat=nml_error)
                if (nml_error /= 0) then
                  call endrun(msg='ERROR reading ndep namelist')
                end if
              end if
              close(nu_nml)
              call relavu( nu_nml )

              ierr = nf90_open(trim(stream_fldFileName_ndep), nf90_nowrite, ncid)
              ierr = nf90_inq_varid(ncid, 'lat', varid)
              ierr = nf90_get_var(ncid, varid, smap2_lat)
              ierr = nf90_inq_varid(ncid, 'lon', varid)
              ierr = nf90_get_var(ncid, varid, smap2_lon)
              ierr = nf90_inq_varid(ncid, 'NDEP_year', varid)
              starti(1:2) = 1
              starti(3)   = nindex(1)
              counti(1)   = 144
              counti(2)   = 96
              counti(3)   = 1
              ierr = nf90_get_var(ncid, varid, a2l%ndep1, starti, counti)
              if (nindex(1) .ne. nindex(2)) then
                starti(3) = nindex(2)
                ierr = nf90_get_var(ncid, varid, a2l%ndep2, starti, counti)
              else
                a2l%ndep2 = a2l%ndep1
              end if
              ierr = nf90_close(ncid)
             end if
             if (i .eq. 1) then
               call mpi_bcast (a2l%ndep1, 144*96, MPI_REAL8, 0, mpicom, ier)
               call mpi_bcast (a2l%ndep2, 144*96, MPI_REAL8, 0, mpicom, ier)
               call mpi_bcast (smap2_lon, 144, MPI_REAL8, 0, mpicom, ier)
               call mpi_bcast (smap2_lat, 96, MPI_REAL8, 0, mpicom, ier)
             end if
          end if

          if (a2l%loaded_bypassdata .eq. 0) then
            mindist=99999
            do thisx = 1,144
              do thisy = 1,96
                if (ldomain%lonc(g) .lt. 0) then
                  if (smap2_lon(thisx) >= 180) smap2_lon(thisx) = smap2_lon(thisx)-360._r8
                else if (ldomain%lonc(g) .ge. 180) then
                  if (smap2_lon(thisx) < 0) smap2_lon(thisx) = smap2_lon(thisx) + 360._r8
                end if
                thislon = smap2_lon(thisx)
                thisdist = 100*((smap2_lat(thisy) - ldomain%latc(g))**2 + &
                              (thislon - ldomain%lonc(g))**2)**0.5
                if (thisdist .lt. mindist) then
                  mindist = thisdist
                  a2l%ndepind(g,1) = thisx
                  a2l%ndepind(g,2) = thisy
                end if
              end do
            end do
          end if

          !get weights for interpolation
          wt1(1) = 1._r8 - (thiscalday -1._r8)/365._r8
          wt2(1) = 1._r8 - wt1(1)

          a2l%forc_ndep(g)    = (a2l%ndep1(a2l%ndepind(g,1),a2l%ndepind(g,2),1)*wt1(1) + &
                                              a2l%ndep2(a2l%ndepind(g,1),a2l%ndepind(g,2),1)*wt2(1)) / (365._r8 * 86400._r8)
        end if

   !------------------------------------Aerosol forcing--------------------------------------------------
        if (a2l%loaded_bypassdata .eq. 0 .or. (mon .eq. 1 .and. day .eq. 1 .and. tod .eq. 0)) then
          if (masterproc .and. i .eq. 1) then
            aerovars(1) = 'BCDEPWET'
            aerovars(2) = 'BCPHODRY'
            aerovars(3) = 'BCPHIDRY'
            aerovars(4) = 'OCDEPWET'
            aerovars(5) = 'OCPHODRY'
            aerovars(6) = 'OCPHIDRY'
            aerovars(7) = 'DSTX01DD'
            aerovars(8) = 'DSTX02DD'
            aerovars(9) = 'DSTX03DD'
            aerovars(10) = 'DSTX04DD'
            aerovars(11) = 'DSTX01WD'
            aerovars(12) = 'DSTX02WD'
            aerovars(13) = 'DSTX03WD'
            aerovars(14) = 'DSTX04WD'
            ierr = nf90_open(trim(aero_file), nf90_nowrite, ncid)
            ierr = nf90_inq_varid(ncid, 'lat', varid)
            ierr = nf90_get_var(ncid, varid, smap2_lat)
            ierr = nf90_inq_varid(ncid, 'lon', varid)
            ierr = nf90_get_var(ncid, varid, smap2_lon)
            starti(1:2) = 1
            starti(3)   = max((min(yr,2100)-1849)*12+1, 13)-1
            counti(1)   = 144
            counti(2)   = 96
            counti(3)   = 14
            do av=1,14
              ierr = nf90_inq_varid(ncid, trim(aerovars(av)), varid)
              ierr = nf90_get_var(ncid, varid, a2l%aerodata(av,:,:,:), starti, counti)
            end do
            ierr = nf90_close(ncid)
          end if
          if (i .eq. 1) then
             call mpi_bcast (a2l%aerodata, 14*144*96*14, MPI_REAL8, 0, mpicom, ier)
          end if
        end if

        !Use ndep grid indices since they're on the same grid
        if (a2l%loaded_bypassdata .eq. 0) then
            mindist=99999
            do thisx = 1,144
              do thisy = 1,96
                if (ldomain%lonc(g) .lt. 0) then
                  if (smap2_lon(thisx) >= 180) smap2_lon(thisx) = smap2_lon(thisx)-360._r8
                else if (ldomain%lonc(g) .ge. 180) then
                  if (smap2_lon(thisx) < 0) smap2_lon(thisx) = smap2_lon(thisx) + 360._r8
                end if
                thislon = smap2_lon(thisx)
                thisdist = 100*((smap2_lat(thisy) - ldomain%latc(g))**2 + &
                              (thislon - ldomain%lonc(g))**2)**0.5
                if (thisdist .lt. mindist) then
                  mindist = thisdist
                  a2l%ndepind(g,1) = thisx
                  a2l%ndepind(g,2) = thisy
                end if
              end do
            end do
        end if

        !get weights for interpolation (note this method doesn't get the month boundaries quite right..)
        aindex(1) = mon+1
        if (thiscalday .le. (caldaym(mon+1)+caldaym(mon))/2._r8) then
           wt1(1) = 0.5_r8 + (thiscalday-caldaym(mon))/(caldaym(mon+1)-caldaym(mon))
           aindex(2) = aindex(1)-1
        else
           wt1(1) = 1.0_r8 - (thiscalday-(caldaym(mon+1)+caldaym(mon))/2._r8)/   &
                          (caldaym(mon+1)-caldaym(mon))
           aindex(2) = aindex(1)+1
        end if
        wt2(1) = 1._r8 - wt1(1)

        do av = 1,14
          a2l%forc_aer(g,av)  =  a2l%aerodata(av,a2l%ndepind(g,1), &
            a2l%ndepind(g,2),aindex(1))*wt1(1)+a2l%aerodata(av,a2l%ndepind(g,1), &
            a2l%ndepind(g,2),aindex(2))*wt2(1)
        end do

  !---------------------------------- ADD T & CO2 --------------------------------------------------------
       !Parse startdate for adding temperature
       if (startdate_add_temperature .ne. '') then
         call get_curr_date( yr, mon, day, tod )
         read(startdate_add_temperature,*) sdate_addt
         sy_addt     = sdate_addt/10000
         sm_addt     = (sdate_addt-sy_addt*10000)/100
         sd_addt     = sdate_addt-sy_addt*10000-sm_addt*100
         read(startdate_add_co2,*) sdate_addco2
         sy_addco2     = sdate_addco2/10000
         sm_addco2     = (sdate_addco2-sy_addco2*10000)/100
         sd_addco2     = sdate_addco2-sy_addco2*10000-sm_addt*100
       end if
       if (startdate_add_temperature .ne. '') then
         if ((yr == sy_addt .and. mon == sm_addt .and. day >= sd_addt) .or. &
             (yr == sy_addt .and. mon > sm_addt) .or. (yr > sy_addt)) then
           a2l%forc_t(g) = a2l%forc_t(g) + add_temperature
           a2l%forc_th(g) = a2l%forc_th(g) + add_temperature
         end if
       end if

  !---------------------------------- CO2 -------------------------------------------------------------------

       if (co2_type_idx /= 0) then
          !atmospheric CO2 (to be used for transient simulations only)
          if (a2l%loaded_bypassdata .eq. 0) then
            ierr = nf90_open(trim(co2_file), nf90_nowrite, ncid)
            if (ierr .ne. 0) call endrun(msg=' ERROR: Failed to open cpl_bypass input CO2 file' )
            ierr = nf90_inq_dimid(ncid, 'time', dimid)
            ierr = nf90_Inquire_Dimension(ncid, dimid, len = thistimelen)
            ierr = nf90_inq_varid(ncid, 'CO2', varid)
            ierr = nf90_get_var(ncid, varid, a2l%co2_input(:,:,1:thistimelen))
            ierr = nf90_inq_varid(ncid, 'C13O2', varid)
            ierr = nf90_get_var(ncid, varid, a2l%c13o2_input(:,:,1:thistimelen))
            ierr = nf90_close(ncid)
          end if

          !get weights/indices for interpolation (assume values represent annual averages)
          nindex(1) = min(max(yr,1850),2100)-1764
          if (thiscalday .le. 182.5) then
            nindex(2) = nindex(1)-1
          else
            nindex(2) = nindex(1)+1
          end if
          wt1(1) = 1._r8 - abs((182.5 - (thiscalday -1._r8))/365._r8)
          wt2(1) = 1._r8 - wt1(1)

          co2_ppmv_val = a2l%co2_input(1,1,nindex(1))*wt1(1) + a2l%co2_input(1,1,nindex(2))*wt2(1)
          if (startdate_add_co2 .ne. '' .and. add_co2 .ne. 0._r8) then
            if ((yr == sy_addco2 .and. mon == sm_addco2 .and. day >= sd_addco2) .or. &
                (yr == sy_addco2 .and. mon > sm_addco2) .or. (yr > sy_addco2)) then
              co2_ppmv_val = co2_ppmv_val + add_co2
            end if
          end if
          if (use_c13) then
            a2l%forc_pc13o2(g) = (a2l%c13o2_input(1,1,nindex(1))*wt1(1) + &
               a2l%c13o2_input(1,1,nindex(2))*wt2(1)) * 1.e-6_r8 * a2l%forc_pbot(g)
          end if
          !TEST (FACE-like experiment begins in 2010)
          !if (yr .ge. 2010) a2l%co2_input = 550.

       else if (co2_type_idx == 0) then

          ! CO2 constant, value from namelist
          co2_ppmv_val = co2_ppmv
          if (use_c13) then
            a2l%forc_pc13o2(g) = co2_ppmv_val * c13ratio * 1.e-6_r8 &
                                             * a2l%forc_pbot(g)
          end if

       else

         call endrun( sub//' ERROR: Invalid co2_type_idx, must be 0 or not (constant or diagnostic) for CPL_BYPASS' )

       end if
     
  !-----------------------------------------------------------------------------------------------------
#else
        a2l%forc_hgt(g)     = x2l_l%rAttr(index_x2l_Sa_z,i)         ! zgcmxy  Atm state m
        a2l%forc_u(g)       = x2l_l%rAttr(index_x2l_Sa_u,i)         ! forc_uxy  Atm state m/s
        a2l%forc_v(g)       = x2l_l%rAttr(index_x2l_Sa_v,i)         ! forc_vxy  Atm state m/s
        a2l%forc_th(g)      = x2l_l%rAttr(index_x2l_Sa_ptem,i)      ! forc_thxy Atm state K
        a2l%forc_q(g)       = x2l_l%rAttr(index_x2l_Sa_shum,i)      ! forc_qxy  Atm state kg/kg
        a2l%forc_pbot(g)    = x2l_l%rAttr(index_x2l_Sa_pbot,i)      ! ptcmxy  Atm state Pa
        a2l%forc_t(g)       = x2l_l%rAttr(index_x2l_Sa_tbot,i)      ! forc_txy  Atm state K
        a2l%forc_lwrad(g)   = x2l_l%rAttr(index_x2l_Faxa_lwdn,i)    ! flwdsxy Atm flux  W/m^2
        forc_rainc          = x2l_l%rAttr(index_x2l_Faxa_rainc,i)   ! mm/s
        forc_rainl          = x2l_l%rAttr(index_x2l_Faxa_rainl,i)   ! mm/s
        forc_snowc          = x2l_l%rAttr(index_x2l_Faxa_snowc,i)   ! mm/s
        forc_snowl          = x2l_l%rAttr(index_x2l_Faxa_snowl,i)   ! mm/s
        a2l%forc_solad(g,2) = x2l_l%rAttr(index_x2l_Faxa_swndr,i)   ! forc_sollxy  Atm flux  W/m^2
        a2l%forc_solad(g,1) = x2l_l%rAttr(index_x2l_Faxa_swvdr,i)   ! forc_solsxy  Atm flux  W/m^2
        a2l%forc_solai(g,2) = x2l_l%rAttr(index_x2l_Faxa_swndf,i)   ! forc_solldxy Atm flux  W/m^2
        a2l%forc_solai(g,1) = x2l_l%rAttr(index_x2l_Faxa_swvdf,i)   ! forc_solsdxy Atm flux  W/m^2

        ! atmosphere coupling, for prognostic/prescribed aerosols
        a2l%forc_aer(g,1)  =  x2l_l%rAttr(index_x2l_Faxa_bcphidry,i)
        a2l%forc_aer(g,2)  =  x2l_l%rAttr(index_x2l_Faxa_bcphodry,i)
        a2l%forc_aer(g,3)  =  x2l_l%rAttr(index_x2l_Faxa_bcphiwet,i)
        a2l%forc_aer(g,4)  =  x2l_l%rAttr(index_x2l_Faxa_ocphidry,i)
        a2l%forc_aer(g,5)  =  x2l_l%rAttr(index_x2l_Faxa_ocphodry,i)
        a2l%forc_aer(g,6)  =  x2l_l%rAttr(index_x2l_Faxa_ocphiwet,i)
        a2l%forc_aer(g,7)  =  x2l_l%rAttr(index_x2l_Faxa_dstwet1,i)
        a2l%forc_aer(g,8)  =  x2l_l%rAttr(index_x2l_Faxa_dstdry1,i)
        a2l%forc_aer(g,9)  =  x2l_l%rAttr(index_x2l_Faxa_dstwet2,i)
        a2l%forc_aer(g,10) =  x2l_l%rAttr(index_x2l_Faxa_dstdry2,i)
        a2l%forc_aer(g,11) =  x2l_l%rAttr(index_x2l_Faxa_dstwet3,i)
        a2l%forc_aer(g,12) =  x2l_l%rAttr(index_x2l_Faxa_dstdry3,i)
        a2l%forc_aer(g,13) =  x2l_l%rAttr(index_x2l_Faxa_dstwet4,i)
        a2l%forc_aer(g,14) =  x2l_l%rAttr(index_x2l_Faxa_dstdry4,i)

        ! Determine optional receive fields

        if (index_x2l_Sa_co2prog /= 0) then
           co2_ppmv_prog = x2l_l%rAttr(index_x2l_Sa_co2prog,i)   ! co2 atm state prognostic
        else
           co2_ppmv_prog = co2_ppmv
        end if
 
        if (index_x2l_Sa_co2diag /= 0) then
           co2_ppmv_diag = x2l_l%rAttr(index_x2l_Sa_co2diag,i)   ! co2 atm state diagnostic
        else
           co2_ppmv_diag = co2_ppmv
        end if

        if (co2_type_idx == 1) then
           co2_ppmv_val = co2_ppmv_prog
        else if (co2_type_idx == 2) then

           co2_ppmv_val = co2_ppmv_diag
           if (use_c13) then
             a2l%forc_pc13o2(g) = co2_ppmv_val * c13ratio * 1.e-6_r8 * a2l%forc_pbot(g)
           end if

           if (yr .ge. startyear_experiment .and. yr .le. endyear_experiment) then
              co2_ppmv_val = co2_ppmv_val + add_co2
           end if

        else
           co2_ppmv_val = co2_ppmv
           if (use_c13) then
              a2l%forc_pc13o2(g) = co2_ppmv_val * c13ratio * 1.e-6_r8 * a2l%forc_pbot(g)
           end if
        end if

#endif


#ifdef LCH4
        if (index_x2l_Sa_methane /= 0) then
           a2l%forc_pch4(g) = x2l_l%rAttr(index_x2l_Sa_methane,i)
        endif
#endif

        ! Determine derived quantities for required fields
        a2l%forc_hgt_u(g) = a2l%forc_hgt(g)    !observational height of wind [m]
        a2l%forc_hgt_t(g) = a2l%forc_hgt(g)    !observational height of temperature [m]
        a2l%forc_hgt_q(g) = a2l%forc_hgt(g)    !observational height of humidity [m]
        a2l%forc_vp(g)    = a2l%forc_q(g) * a2l%forc_pbot(g) &
                            / (0.622_r8 + 0.378_r8 * a2l%forc_q(g))
        a2l%forc_rho(g)   = (a2l%forc_pbot(g) - 0.378_r8 * a2l%forc_vp(g)) &
                            / (rair * a2l%forc_t(g))
        a2l%forc_po2(g)   = o2_molar_const * a2l%forc_pbot(g)
        a2l%forc_wind(g)  = sqrt(a2l%forc_u(g)**2 + a2l%forc_v(g)**2)
        a2l%forc_solar(g) = a2l%forc_solad(g,1) + a2l%forc_solai(g,1) + &
                            a2l%forc_solad(g,2) + a2l%forc_solai(g,2)
        a2l%forc_rain(g)  = forc_rainc + forc_rainl
        a2l%forc_snow(g)  = forc_snowc + forc_snowl
        a2l%rainf    (g)  = a2l%forc_rain(g) + a2l%forc_snow(g)

        if (a2l%forc_t(g) > SHR_CONST_TKFRZ) then
           e = esatw(tdc(a2l%forc_t(g)))
        else
           e = esati(tdc(a2l%forc_t(g)))
        end if
        qsat           = 0.622_r8*e / (a2l%forc_pbot(g) - 0.378_r8*e)
        a2l%forc_rh(g) = 100.0_r8*(a2l%forc_q(g) / qsat)
        ! Make sure relative humidity is properly bounded
        a2l%forc_rh(g) = min( 100.0_r8, a2l%forc_rh(g) )
        a2l%forc_rh(g) = max(   0.0_r8, a2l%forc_rh(g) )
        
        ! Determine derived quantities for optional fields
        ! Note that the following does unit conversions from ppmv to partial pressures (Pa)
        ! Note that forc_pbot is in Pa

        a2l%forc_pco2(g)   = co2_ppmv_val * 1.e-6_r8 * a2l%forc_pbot(g)

     end do

#ifdef CPL_BYPASS
    a2l%loaded_bypassdata = 1
#endif


   end subroutine lnd_import_mct

!===============================================================================

!---------------------------------------------------------------------------
!BOP
!
! !IROUTINE: lnd_domain_mct
!
! !INTERFACE:
  subroutine lnd_domain_mct( lsize, gsMap_l, dom_l )
!
! !DESCRIPTION:
!
! Send the land model domain information to the coupler
!
!---------------------------------------------------------------------------
! !USES:
    use shr_kind_mod, only : r8 => shr_kind_r8
    use clm_varcon  , only : re
    use domainMod   , only : ldomain
    use decompMod   , only : get_proc_bounds
    use spmdMod     , only : iam
    use mct_mod     , only : mct_gsMap, mct_gGrid, mct_gGrid_importIAttr, &
                             mct_gGrid_importRAttr, mct_gGrid_init,       &
                             mct_gsMap_orderedPoints
    use seq_flds_mod
    implicit none
! !ARGUMENTS:
    integer        , intent(in)    :: lsize     ! land model domain data size
    type(mct_gsMap), intent(inout) :: gsMap_l   ! Output land model MCT GS map
    type(mct_ggrid), intent(out)   :: dom_l     ! Output domain information for land model
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!---------------------------------------------------------------------------
    !
    ! Local Variables
    !
    integer :: g,i,j              ! index
    integer :: begg, endg         ! beginning and ending gridcell indices
    real(r8), pointer :: data(:)  ! temporary
    integer , pointer :: idata(:) ! temporary
    !-------------------------------------------------------------------
    !
    ! Initialize mct domain type
    ! lat/lon in degrees,  area in radians^2, mask is 1 (land), 0 (non-land)
    ! Note that in addition land carries around landfrac for the purposes of domain checking
    ! 
    call mct_gGrid_init( GGrid=dom_l, CoordChars=trim(seq_flds_dom_coord), &
       OtherChars=trim(seq_flds_dom_other), lsize=lsize )
    !
    ! Allocate memory
    !
    allocate(data(lsize))
    !
    ! Determine global gridpoint number attribute, GlobGridNum, which is set automatically by MCT
    !
    call mct_gsMap_orderedPoints(gsMap_l, iam, idata)
    call mct_gGrid_importIAttr(dom_l,'GlobGridNum',idata,lsize)
    !
    ! Determine domain (numbering scheme is: West to East and South to North to South pole)
    ! Initialize attribute vector with special value
    !
    data(:) = -9999.0_R8 
    call mct_gGrid_importRAttr(dom_l,"lat"  ,data,lsize) 
    call mct_gGrid_importRAttr(dom_l,"lon"  ,data,lsize) 
    call mct_gGrid_importRAttr(dom_l,"area" ,data,lsize) 
    call mct_gGrid_importRAttr(dom_l,"aream",data,lsize) 
    data(:) = 0.0_R8     
    call mct_gGrid_importRAttr(dom_l,"mask" ,data,lsize) 
    !
    ! Determine bounds
    !
    call get_proc_bounds(begg, endg)
    !
    ! Fill in correct values for domain components
    ! Note aream will be filled in in the atm-lnd mapper
    !
    do g = begg,endg
       i = 1 + (g - begg)
       data(i) = ldomain%lonc(g)
    end do
    call mct_gGrid_importRattr(dom_l,"lon",data,lsize) 

    do g = begg,endg
       i = 1 + (g - begg)
       data(i) = ldomain%latc(g)
    end do
    call mct_gGrid_importRattr(dom_l,"lat",data,lsize) 

    do g = begg,endg
       i = 1 + (g - begg)
       data(i) = ldomain%area(g)/(re*re)
    end do
    call mct_gGrid_importRattr(dom_l,"area",data,lsize) 

    do g = begg,endg
       i = 1 + (g - begg)
       data(i) = real(ldomain%mask(g), r8)
    end do
    call mct_gGrid_importRattr(dom_l,"mask",data,lsize) 

    do g = begg,endg
       i = 1 + (g - begg)
       data(i) = real(ldomain%frac(g), r8)
    end do
    call mct_gGrid_importRattr(dom_l,"frac",data,lsize) 

    deallocate(data)
    deallocate(idata)

  end subroutine lnd_domain_mct
    
!===============================================================================

  subroutine sno_export_mct( s2x, s2x_s )   

    use clm_glclnd      , only : lnd2glc_type
    use decompMod       , only : get_proc_bounds
    use clm_cpl_indices 
    use clm_varctl       , only : iulog

    type(lnd2glc_type), intent(inout) :: s2x
    type(mct_aVect)   , intent(inout) :: s2x_s

    integer :: g,i,num
    integer :: begg, endg    ! beginning and ending gridcell indices
    !-----------------------------------------------------

    call get_proc_bounds(begg, endg)

    ! qice is positive if ice is growing, negative if melting

    s2x_s%rAttr(:,:) = 0.0_r8
    do g = begg,endg
       i = 1 + (g-begg)
       do num = 1,glc_nec
          s2x_s%rAttr(index_s2x_Ss_tsrf(num),i)   = s2x%tsrf(g,num)
          s2x_s%rAttr(index_s2x_Ss_topo(num),i)   = s2x%topo(g,num)
          s2x_s%rAttr(index_s2x_Fgss_qice(num),i) = s2x%qice(g,num)
       end do
    end do  

  end subroutine sno_export_mct

!====================================================================================

  subroutine sno_import_mct( x2s_s, x2s )

    use clm_glclnd      , only: glc2lnd_type
    use decompMod       , only: get_proc_bounds
    use mct_mod         , only: mct_aVect
    use clm_cpl_indices
    !
    ! Arguments
    !
    type(mct_aVect)   , intent(inout) :: x2s_s
    type(glc2lnd_type), intent(inout) :: x2s
    !
    ! Local Variables
    !
    integer  :: g,i,num
    integer  :: begg, endg   ! beginning and ending gridcell indices
    !-----------------------------------------------------

    call get_proc_bounds(begg, endg)

    do g = begg,endg
       i = 1 + (g - begg)
       do num = 1,glc_nec
          x2s%frac(g,num)  = x2s_s%rAttr(index_x2s_Sg_frac(num),i)
          x2s%topo(g,num)  = x2s_s%rAttr(index_x2s_Sg_topo(num),i)
          x2s%hflx(g,num)  = x2s_s%rAttr(index_x2s_Fsgg_hflx(num),i)
          x2s%rofi(g,num)  = x2s_s%rAttr(index_x2s_Fsgg_rofi(num),i)
          x2s%rofl(g,num)  = x2s_s%rAttr(index_x2s_Fsgg_rofl(num),i)
       end do
     end do  

   end subroutine sno_import_mct

!====================================================================================

end module lnd_comp_mct



double precision function szenith(xcoor, ycoor, ltm, jday, hr, min, offset)
  !Function to calcualte solar zenith angle
  !Used in coupler bypass mode to compute inerpolation for incoming solar

  use shr_kind_mod , only: r8 => shr_kind_r8, cl=>shr_kind_cl
  implicit none
  !inputs
  real(r8) xcoor, ycoor, offset_min
  integer jday, hr, min, ltm, offset
  !working variables
  real(r8) d2r, r2d, lsn, latrad, decrad, decdeg, ha
  real(r8) hangle, harad, saltrad, saltdeg, sazirad, sazideg
  real(r8) szendeg,szenrad

  real pi
  parameter(pi = 3.14159265358979)
  offset_min = offset/60d0   !note assumes 1hr or smaller timestep
  min = min - offset_min

  !adjust time for offsets
  if (min < 0) then
    hr = hr - 1
    min = min+60
  end if
  if (min >= 60) then
    hr = hr+1
    min = min-60
  end if
  if (hr < 0) then
    hr = hr+24
    jday = jday-1
  end if
  if (hr >= 24) then
    hr = hr-24
    jday = jday+1
  end if

  if (jday < 1) jday = 1
  if (xcoor > 180d0) xcoor = xcoor-360d0

  d2r     = pi/180d0
  r2d     = 1/d2r
  lsn     = 12.0d0+((ltm-xcoor)/15.0d0)
  latrad  = ycoor*d2r
  decrad  = 23.45*d2r*sin(d2r*360d0*(284d0+jday)/365d0)
  decdeg  = decrad*r2d
  ha      = hr+min/60.0d0
  hangle  = (lsn-ha)*60.0d0
  harad   = hangle*0.0043633d0

  saltrad = asin((sin(latrad)*sin(decrad))+(cos(latrad)*cos(decrad) &
       *cos(harad)))
  saltdeg = saltrad * r2d
  sazirad = asin(cos(decrad)*sin(harad)/cos(saltrad))
  sazideg = sazirad * r2d

  IF (saltdeg.LT.0.0d0 .OR. saltrad.GT.180.0d0) THEN  ! sun is below horizon
     saltdeg = 0.0d0
     saltrad = 0.0d0
     szendeg = 90.0d0
     szenrad = 90.0d0*d2r
     !mass    = 1229d0**.5d0             ! if solaralt=0 -> sin(0)=0
  ELSE
     szendeg = 90d0-saltdeg
     szenrad = szendeg*d2r
  ENDIF
  szenith = szendeg

end function szenith
