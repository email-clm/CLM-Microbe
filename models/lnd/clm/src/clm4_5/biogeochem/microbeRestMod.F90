module microbeRestMod
#ifdef MICROBE

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: microbeRestMod
!
! !DESCRIPTION:
! Reads from or writes restart data
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils,   only : endrun
!
! !PUBLIC TYPES:
  implicit none
  private
! save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: microbeRest
!
! !REVISION HISTORY:
! Created by Xiaofeng Xu in 2013
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: micrbeRest
!
! !INTERFACE:
  subroutine microbeRest( ncid, flag )
!
! !DESCRIPTION:
! Read/Write biogeophysics information to/from restart file.
!
! !USES:
    use clmtype
    use ncdio_pio
    use decompMod     , only : get_proc_bounds
    use clm_varctl    , only : nsrest
    use clm_time_manager , only : is_restart
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(inout) :: ncid ! netcdf id
    character(len=*), intent(in) :: flag     ! 'read' or 'write'
!
! !CALLED FROM:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: c,l,g,j      ! indices
    integer :: begp, endp   ! per-proc beginning and ending pft indices
    integer :: begc, endc   ! per-proc beginning and ending column indices
    integer :: begl, endl   ! per-proc beginning and ending landunit indices
    integer :: begg, endg   ! per-proc gridcell ending gridcell indices
    logical :: readvar      ! determine if variable is on initial file
    character(len=128) :: varname         ! temporary
    !type(gridcell_type), pointer :: gptr  ! pointer to gridcell derived subtype
    type(landunit_type), pointer :: lptr  ! pointer to landunit derived subtype
    type(column_type)  , pointer :: cptr  ! pointer to column derived subtype
    type(pft_type)     , pointer :: pptr  ! pointer to pft derived subtype
!-----------------------------------------------------------------------

    ! Set pointers into derived type

    !gptr => clm3%g
    lptr => lun !clm3%g%l
    cptr => col !clm3%g%l%c
    pptr => pft !clm3%g%l%c%p
!    finundated					=> cws%finundated

    ! column microbial state variable - cdocs
  if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='pH', xtype=ncd_double, &
            dim1name='column', long_name='pH', units='unitless')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='pH', data=cps%pH, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='soilpH_unsat', xtype=ncd_double, &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='soilpH in unsat', units='unitless')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='soilpH_unsat', data=cps%soilpH_unsat, &
            dim1name='column', switchdim=.true., &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='soilpH_sat', xtype=ncd_double, &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='soilpH in saturation fraction', units='unitless')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='soilpH_sat', data=cps%soilpH_sat, &
            dim1name='column', switchdim=.true., &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

   if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='FINUNDATED', xtype=ncd_double, &
            dim1name='column', long_name='fraction of inundated area', units='m^2/m^2')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='FINUNDATED', data=cws%finundated, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

   if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='fsat_pre', xtype=ncd_double, &
            dim1name='column', long_name='previous time step fraction of inundated area', units='m^2/m^2')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='fsat_pre', data=cmic%fsat_pre, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

   if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='fsat', xtype=ncd_double, &
            dim1name='column', long_name='fraction of inundated area', units='m^2/m^2')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='fsat', data=cws%fsat, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

   if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='waterhead_unsat', xtype=ncd_double, &
            dim1name='column', long_name='fraction of inundated area', units='m^2/m^2')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='waterhead_unsat', data=cmic%waterhead_unsat, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if
!20150831
   !~ if (flag == 'define') then
       !~ call ncd_defvar(ncid=ncid, varname='rr', xtype=ncd_double, &
            !~ dim1name='column', long_name='respiration', units='mol/s')
    !~ else if (flag == 'read' .or. flag == 'write') then
       !~ call ncd_io(varname='rr', data=pcf%rr, &
            !~ dim1name='column', ncid=ncid, flag=flag, readvar=readvar)
       !~ if (flag=='read' .and. .not. readvar) then
          !~ if (is_restart()) call endrun()
       !~ end if
    !~ end if
    
   !~ if (flag == 'define') then
       !~ call ncd_defvar(ncid=ncid, varname='froot_mr', xtype=ncd_double, &
            !~ dim1name='column', long_name='fine root maintenance respiration', units='mol/s')
    !~ else if (flag == 'read' .or. flag == 'write') then
       !~ call ncd_io(varname='froot_mr', data=pcf%froot_mr, &
            !~ dim1name='column', ncid=ncid, flag=flag, readvar=readvar)
       !~ if (flag=='read' .and. .not. readvar) then
          !~ if (is_restart()) call endrun()
       !~ end if
    !~ end if

    !~ if (flag == 'define') then
       !~ call ncd_defvar(ncid=ncid, varname='rootfr', xtype=ncd_double, &
            !~ dim1name='column', dim2name='levgrnd', switchdim=.true., &
            !~ long_name='root rerspiration fraction', units='%')
    !~ else if (flag == 'read' .or. flag == 'write') then
       !~ call ncd_io(varname='rootfr', data=pps%rootfr, &
            !~ dim1name='column', switchdim=.true., &
            !~ ncid=ncid, flag=flag, readvar=readvar)
       !~ if (flag=='read' .and. .not. readvar) then
          !~ if (is_restart()) call endrun()
       !~ end if
    !~ end if
    
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='hr_vr', xtype=ncd_double, &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='hetrotrophic respiration fraction', units='mol/m^2')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='hr_vr', data=ccf%hr_vr, &
            dim1name='column', switchdim=.true., &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if
 !20150831   
    !~ if (flag == 'define') then
       !~ call ncd_defvar(ncid=ncid, varname='sminn_vr', xtype=ncd_double, &
            !~ dim1name='column', dim2name='levgrnd', switchdim=.true., &
            !~ long_name='soil mineral nitrogen', units='mol/m^2')
    !~ else if (flag == 'read' .or. flag == 'write') then
       !~ call ncd_io(varname='sminn_vr', data=cns%sminn_vr, &
            !~ dim1name='column', switchdim=.true., &
            !~ ncid=ncid, flag=flag, readvar=readvar)
       !~ if (flag=='read' .and. .not. readvar) then
          !~ if (is_restart()) call endrun()
       !~ end if
    !~ end if

!20150830
   if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='froot_r', xtype=ncd_double, &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='fine root respiration', units='mol/s')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='froot_r', data=cmic%froot_r, &
            dim1name='column', switchdim=.true., &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if
    
!~ #if (defined HUM_HOL)	
   !~ if (flag == 'define') then
       !~ call ncd_defvar(ncid=ncid, varname='qflx_lat_aqu_layer', xtype=ncd_double, &
            !~ dim1name='column', dim2name='levgrnd', switchdim=.true., &
            !~ long_name='lateral water flux', units='mol/s')
    !~ else if (flag == 'read' .or. flag == 'write') then
       !~ call ncd_io(varname='qflx_lat_aqu_layer', data=cwf%qflx_lat_aqu_layer, &
            !~ dim1name='column', switchdim=.true., &
            !~ ncid=ncid, flag=flag, readvar=readvar)
       !~ if (flag=='read' .and. .not. readvar) then
          !~ if (is_restart()) call endrun()
       !~ end if
    !~ end if
!~ #endif
!20150830
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='CACES_PROD', xtype=ncd_double, &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='acetate production ', units='mol/m^3')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='CACES_PROD', data=cmic%caces_prod, &
            dim1name='column', switchdim=.true., &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if
   
   if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='CACES_UNSAT_PROD', xtype=ncd_double, &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='acetate production in unsaturated fraction', units='mol/m^3')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='CACES_UNSAT_PROD', data=cmic%caces_unsat_prod, &
            dim1name='column', switchdim=.true., &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if
    
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='CACES_SAT_PROD', xtype=ncd_double, &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='acetate production in saturated fraction', units='mol/m^3')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='CACES_SAT_PROD', data=cmic%caces_sat_prod, &
            dim1name='column', switchdim=.true., &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='CACES_PROD_H2', xtype=ncd_double, &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='acetogenesis ', units='mol/m^3')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='CACES_PROD_H2', data=cmic%caces_prod_h2, &
            dim1name='column', switchdim=.true., &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if
   
   if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='CACES_UNSAT_PROD_H2', xtype=ncd_double, &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='acetogenesis in unsaturated fraction', units='mol/m^3')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='CACES_UNSAT_PROD_H2', data=cmic%caces_unsat_prod_h2, &
            dim1name='column', switchdim=.true., &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if
    
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='CACES_SAT_PROD_H2', xtype=ncd_double, &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='acetogenesis in saturated fraction', units='mol/m^3')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='CACES_SAT_PROD_H2', data=cmic%caces_sat_prod_h2, &
            dim1name='column', switchdim=.true., &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='CDOCS_PRE', xtype=ncd_double, &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='previous time step DOC concentration', units='mol/m^3')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='CDOCS_PRE', data=cmic%cdocs_pre, &
            dim1name='column', switchdim=.true., &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if
    
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='CDONS_MIN', xtype=ncd_double, &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='dissolved organic nitrogen minterlization ', units='mol/m^3')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='CDONS_MIN', data=cmic%cdons_min, &
            dim1name='column', switchdim=.true., &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if
    
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='CDOCS', xtype=ncd_double, &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='dissolved organic carbon', units='mol/m^3')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='CDOCS', data=cmic%cdocs, &
            dim1name='column', switchdim=.true., &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='CDOCS_UNSAT', xtype=ncd_double, &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='dissolved organic carbon in unsaturated fraction', units='mol/m^3')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='CDOCS_UNSAT', data=cmic%cdocs_unsat, &
            dim1name='column', switchdim=.true., &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if
    
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='CDOCS_SAT', xtype=ncd_double, &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='dissolved organic carbon in saturated fraction', units='mol/m^3')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='CDOCS_SAT', data=cmic%cdocs_sat, &
            dim1name='column', switchdim=.true., &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if
    
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='CDONS', xtype=ncd_double, &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='dissolved organic nitrogen', units='mol/m^3')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='CDONS', data=cmic%cdons, &
            dim1name='column', switchdim=.true., &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='CDONS_UNSAT', xtype=ncd_double, &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='dissolved organic nitrogen in unsaturated fraction', units='mol/m^3')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='CDONS_UNSAT', data=cmic%cdons_unsat, &
            dim1name='column', switchdim=.true., &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if
    
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='CDONS_SAT', xtype=ncd_double, &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='dissolved organic nitrogen in saturated fraction', units='mol/m^3')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='CDONS_SAT', data=cmic%cdons_sat, &
            dim1name='column', switchdim=.true., &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if
    
    ! column microbial state variable - cmicbiocs
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='CMICBIOCS', xtype=ncd_double, &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='microbial biomass carbon', units='molC/m^3')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='CMICBIOCS', data=cmic%cmicbiocs, &
            dim1name='column', switchdim=.true., &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if
 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='CMICBIONS', xtype=ncd_double, &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='microbial biomass nitrogen', units='molN/m^3')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='CMICBIONS', data=cmic%cmicbions, &
            dim1name='column', switchdim=.true., &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if 
    
    ! column microbial state variable - caces
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='CACES', xtype=ncd_double, &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='acetic acid ', units='mol/m^3')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='CACES', data=cmic%caces, &
            dim1name='column', switchdim=.true., &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if
    
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='CACES_UNSAT', xtype=ncd_double, &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='acetic acid in unsaturated fraction ', units='mol/m^3')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='CACES_UNSAT', data=cmic%caces_unsat, &
            dim1name='column', switchdim=.true., &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if
    
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='CACES_SAT', xtype=ncd_double, &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='acetic acid in saturated fraction', units='mol/m^3')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='CACES_SAT', data=cmic%caces_sat, &
            dim1name='column', switchdim=.true., &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if
    
    ! column microbial state variable - cacebios
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='CACEBIOS', xtype=ncd_double, &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='methanogenesis based on acetic acid', units='mol/m^3')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='CACEBIOS', data=cmic%cacebios, &
            dim1name='column', switchdim=.true., &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if
    
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='CACEBIOS_UNSAT', xtype=ncd_double, &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='methanogenesis based on acetic acid in unsaturated fraction', units='mol/m^3')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='CACEBIOS_UNSAT', data=cmic%cacebios_unsat, &
            dim1name='column', switchdim=.true., &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if
    
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='CACEBIOS_SAT', xtype=ncd_double, &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='methanogenesis based on acetic acid in saturated fraction soil column', units='mol/m^3')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='CACEBIOS_SAT', data=cmic%cacebios_sat, &
            dim1name='column', switchdim=.true., &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if
    
    ! column microbial state variable - cco2bios
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='CCO2BIOS', xtype=ncd_double, &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='methanogensis based on co2', units='mol/m^3')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='CCO2BIOS', data=cmic%cco2bios, &
            dim1name='column', switchdim=.true., &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='CCO2BIOS_UNSAT', xtype=ncd_double, &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='methanogensis based on co2 in unsaturated fraction', units='mol/m^3')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='CCO2BIOS_UNSAT', data=cmic%cco2bios_unsat, &
            dim1name='column', switchdim=.true., &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='CCO2BIOS_SAT', xtype=ncd_double, &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='methanogensis based on co2 in saturated fraction', units='mol/m^3')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='CCO2BIOS_SAT', data=cmic%cco2bios_sat, &
            dim1name='column', switchdim=.true., &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    !~ ! column microbial state variable - caerch4bios
        if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='CAERCH4BIOS', xtype=ncd_double, &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='aerobic oxidationo of ch4', units='mol/m^3')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='CAERCH4BIOS', data=cmic%caerch4bios, &
            dim1name='column', switchdim=.true., &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='CAERCH4BIOS_UNSAT', xtype=ncd_double, &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='aerobic oxidationo of ch4 in unsaturated fraction', units='mol/m^3')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='CAERCH4BIOS_UNSAT', data=cmic%caerch4bios_unsat, &
            dim1name='column', switchdim=.true., &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='CAERCH4BIOS_SAT', xtype=ncd_double, &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='aerobic oxidationo of ch4 in saturated fraction', units='mol/m^3')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='CAERCH4BIOS_SAT', data=cmic%caerch4bios_sat, &
            dim1name='column', switchdim=.true., &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    !~ ! column microbial state variable - canaerch4bios    
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='CANAERCH4BIOS', xtype=ncd_double, &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='biomass for microbes for anaerobic methane oxidation', units='mol/m^3')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='CANAERCH4BIOS', data=cmic%canaerch4bios, &
            dim1name='column', switchdim=.true., &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='CANAERCH4BIOS_UNSAT', xtype=ncd_double, &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='biomass for microbes for anaerobic methane oxidation in unsaturated fraction', units='mol/m^3')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='CANAERCH4BIOS_UNSAT', data=cmic%canaerch4bios_unsat, &
            dim1name='column', switchdim=.true., &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='CANAERCH4BIOS_SAT', xtype=ncd_double, &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='biomass for microbes for anaerobic methane oxidation in saturated fraction', units='mol/m^3')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='CANAERCH4BIOS_SAT', data=cmic%canaerch4bios_sat, &
            dim1name='column', switchdim=.true., &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if
   
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='CCON_CH4S', xtype=ncd_double, &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='column ch4 concentration', units='mol/m^3')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='CCON_CH4S', data=cmic%ccon_ch4s, &
            dim1name='column', switchdim=.true., &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if
    
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='CCON_CH4S_UNSAT', xtype=ncd_double, &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='column ch4 concentration in unsaturated fraction', units='mol/m^3')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='CCON_CH4S_UNSAT', data=cmic%ccon_ch4s_unsat, &
            dim1name='column', switchdim=.true., &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='CCON_CH4S_SAT', xtype=ncd_double, &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='column ch4 concentration in saturated fraction', units='mol/m^3')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='CCON_CH4S_SAT', data=cmic%ccon_ch4s_sat, &
            dim1name='column', switchdim=.true., &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='CCON_CO2S', xtype=ncd_double, &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='column co2 concentration', units='mol/m^3')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='CCON_CO2S', data=cmic%ccon_co2s, &
            dim1name='column', switchdim=.true., &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if
    
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='CCON_CO2S_UNSAT', xtype=ncd_double, &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='column co2 concentration in unsaturated fraction', units='mol/m^3')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='CCON_CO2S_UNSAT', data=cmic%ccon_co2s_unsat, &
            dim1name='column', switchdim=.true., &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='CCON_CO2S_SAT', xtype=ncd_double, &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='column co2 concentration in saturated fraction', units='mol/m^3')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='CCON_CO2S_SAT', data=cmic%ccon_co2s_sat, &
            dim1name='column', switchdim=.true., &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if


    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='CCON_H2S', xtype=ncd_double, &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='column h2 concentration', units='mol/m^3')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='CCON_H2S', data=cmic%ccon_h2s, &
            dim1name='column', switchdim=.true., &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if
    
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='CCON_H2S_UNSAT', xtype=ncd_double, &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='column h2 concentration in unsaturated fraction', units='mol/m^3')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='CCON_H2S_UNSAT', data=cmic%ccon_h2s_unsat, &
            dim1name='column', switchdim=.true., &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if
    
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='CCON_H2S_SAT', xtype=ncd_double, &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='column h2 concentration in saturated fraction', units='mol/m^3')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='CCON_H2S_SAT', data=cmic%ccon_h2s_sat, &
            dim1name='column', switchdim=.true., &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if
    
    
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='CCON_O2S', xtype=ncd_double, &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='column o2 concentration', units='mol/m^3')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='CCON_O2S', data=cmic%ccon_o2s, &
            dim1name='column', switchdim=.true., &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if
    
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='CCON_O2S_UNSAT', xtype=ncd_double, &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='column o2 concentration in unsaturated fraction', units='mol/m^3')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='CCON_O2S_UNSAT', data=cmic%ccon_o2s_unsat, &
            dim1name='column', switchdim=.true., &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if
    
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='CCON_O2S_SAT', xtype=ncd_double, &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='column o2 concentration in saturated fraction', units='mol/m^3')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='CCON_O2S_SAT', data=cmic%ccon_o2s_sat, &
            dim1name='column', switchdim=.true., &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if
    
	
    ! pft ch4 state variable - tempavg_agnpp
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='tempavg_agnpp', xtype=ncd_double,  &
            dim1name='pft',long_name='Temp. Average AGNPP',units='gC/m^2/s')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='tempavg_agnpp', data=pcf%tempavg_agnpp, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun
       end if
    end if

    ! pft ch4 state variable - tempavg_bgnpp
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='tempavg_bgnpp', xtype=ncd_double,  &
            dim1name='pft',long_name='Temp. Average BGNPP',units='gC/m^2/s')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='tempavg_bgnpp', data=pcf%tempavg_bgnpp, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun
       end if
    end if

    ! pft ch4 state variable - annavg_agnpp
   if (flag == 'define') then
      call ncd_defvar(ncid=ncid, varname='annavg_agnpp', xtype=ncd_double,  &
            dim1name='pft',long_name='Ann. Average AGNPP',units='gC/m^2/s')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='annavg_agnpp', data=pcf%annavg_agnpp, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun
       end if
    end if

    ! pft ch4 state variable - annavg_bgnpp
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='annavg_bgnpp', xtype=ncd_double,  &
            dim1name='pft',long_name='Ann. Average BGNPP',units='gC/m^2/s')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='annavg_bgnpp', data=pcf%annavg_bgnpp, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun
       end if
    end if
    
  end subroutine microbeRest

#endif

end module microbeRestMod
