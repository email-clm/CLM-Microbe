module seq_flds_mod

   use shr_sys_mod,    only : shr_sys_abort
   use seq_drydep_mod, only : seq_drydep_init, seq_drydep_read, lnd_drydep

   implicit none
   public
   save

   interface seq_flds_lookup; module procedure &
     seq_flds_esmf_metadata_get
   end interface

   integer, parameter :: CX = 800           ! use local extra-long char
   character(len=CX)  :: seq_drydep_fields  ! List of dry-deposition fields

   !----------------------------------------------------------------------------
   ! for the domain
   !----------------------------------------------------------------------------

   character(CX) :: seq_flds_dom_coord = 'lat:lon'
   character(CX) :: seq_flds_dom_other = 'area:aream:mask:frac:ascale'

   !----------------------------------------------------------------------------
   ! state + flux fields
   !----------------------------------------------------------------------------

   character(CX) :: seq_flds_a2x_states 
   character(CX) :: seq_flds_a2x_fluxes 
   character(CX) :: seq_flds_x2a_states 
   character(CX) :: seq_flds_x2a_fluxes

   character(CX) :: seq_flds_i2x_states 
   character(CX) :: seq_flds_i2x_fluxes 
   character(CX) :: seq_flds_x2i_states 
   character(CX) :: seq_flds_x2i_fluxes

   character(CX) :: seq_flds_l2x_states 
   character(CX) :: seq_flds_l2x_fluxes 
   character(CX) :: seq_flds_x2l_states 
   character(CX) :: seq_flds_x2l_fluxes

   character(CX) :: seq_flds_o2x_states 
   character(CX) :: seq_flds_o2x_fluxes 
   character(CX) :: seq_flds_x2o_states 
   character(CX) :: seq_flds_x2o_fluxes

   character(CX) :: seq_flds_g2x_states 
   character(CX) :: seq_flds_g2x_fluxes 
   character(CX) :: seq_flds_x2g_states 
   character(CX) :: seq_flds_x2g_fluxes

   character(CX) :: seq_flds_s2x_states 
   character(CX) :: seq_flds_s2x_fluxes 
   character(CX) :: seq_flds_x2s_states 
   character(CX) :: seq_flds_x2s_fluxes

   character(CX) :: seq_flds_xao_albedo
   character(CX) :: seq_flds_xao_states 
   character(CX) :: seq_flds_xao_fluxes

   character(CX) :: seq_flds_r2x_states 
   character(CX) :: seq_flds_r2x_fluxes

   !----------------------------------------------------------------------------
   ! combined state/flux fields
   !----------------------------------------------------------------------------

   character(CX) :: seq_flds_dom_fields 
   character(CX) :: seq_flds_a2x_fields 
   character(CX) :: seq_flds_x2a_fields 
   character(CX) :: seq_flds_i2x_fields 
   character(CX) :: seq_flds_x2i_fields 
   character(CX) :: seq_flds_l2x_fields 
   character(CX) :: seq_flds_x2l_fields 
   character(CX) :: seq_flds_o2x_fields 
   character(CX) :: seq_flds_x2o_fields 
   character(CX) :: seq_flds_xao_fields 
   character(CX) :: seq_flds_r2x_fields
   character(CX) :: seq_flds_g2x_fields 
   character(CX) :: seq_flds_x2g_fields 
   character(CX) :: seq_flds_s2x_fields 
   character(CX) :: seq_flds_x2s_fields 
   character(CX) :: stringtmp

   !----------------------------------------------------------------------------
   ! component names
   !----------------------------------------------------------------------------

   character(32) :: seq_flds_atmname='atm'
   character(32) :: seq_flds_ocnname='ocn'
   character(32) :: seq_flds_icename='ice'
   character(32) :: seq_flds_lndname='lnd'
   character(32) :: seq_flds_glcname='glc'
   character(32) :: seq_flds_snoname='sno'
   character(32) :: seq_flds_rtmname='roff'

!----------------------------------------------------------------------------
 contains
!----------------------------------------------------------------------------

subroutine seq_flds_set()

   !----------------------------------------------------------------------------
   ! atm fields: atm->drv
   !----------------------------------------------------------------------------
   ! DEF = default fields.  CCSM expects at least these fields to be present.
   ! Remove at your own risk.
   ! 
   ! If you add a field, add a name to the string below.  You may also
   ! document the additional field here. The relevant component must
   ! also be modified to set or access the new field.
   !
   !----------------------------------------------------------------------------
   ! States
   seq_flds_a2x_states = &
         'Sa_z'        &    ! bottom atm level height         DEF
      //':Sa_u'        &    ! bottom atm level zon wind       DEF
      //':Sa_v'        &    ! bottom atm level mer wind       DEF
      //':Sa_tbot'     &    ! bottom atm level temp           DEF
      //':Sa_ptem'     &    ! bottom atm level pot temp       DEF
      //':Sa_shum'     &    ! bottom atm level spec hum       DEF
      //':Sa_dens'     &    ! bottom atm level air den        DEF
      //':Sa_pbot'     &    ! bottom atm level pressurea      DEF
      //':Sa_pslv'          ! sea level atm pressure          DEF

   ! Fluxes
   seq_flds_a2x_fluxes = &
         'Faxa_lwdn'   &    ! downward lw heat flux           DEF
      //':Faxa_rainc'  &    ! prec: liquid "convective"       DEF
      //':Faxa_rainl'  &    ! prec: liquid "large scale"      DEF
      //':Faxa_snowc'  &    ! prec: frozen "convective"       DEF
      //':Faxa_snowl'  &    ! prec: frozen "large scale"      DEF
      //':Faxa_swndr'  &    ! sw: nir direct  downward        DEF
      //':Faxa_swvdr'  &    ! sw: vis direct  downward        DEF
      //':Faxa_swndf'  &    ! nir diffuse downward            DEF
      //':Faxa_swvdf'  &    ! sw: vis diffuse downward        DEF
      //':Faxa_swnet'  &    ! sw: net                         DEF
      //':Faxa_bcphidry' &  ! flux: Black Carbon hydrophilic dry deposition   DEF
      //':Faxa_bcphodry' &  ! flux: Black Carbon hydrophobic dry deposition   DEF
      //':Faxa_bcphiwet' &  ! flux: Black Carbon hydrophilic wet deposition   DEF
      //':Faxa_ocphidry' &  ! flux: Organic Carbon hydrophilic dry deposition DEF
      //':Faxa_ocphodry' &  ! flux: Organic Carbon hydrophobic dry deposition DEF
      //':Faxa_ocphiwet' &  ! flux: Organic Carbon hydrophilic dry deposition DEF
      //':Faxa_dstwet1'  &  ! flux: Size 1 dust -- wet deposition DEF
      //':Faxa_dstwet2'  &  ! flux: Size 2 dust -- wet deposition DEF
      //':Faxa_dstwet3'  &  ! flux: Size 3 dust -- wet deposition DEF
      //':Faxa_dstwet4'  &  ! flux: Size 4 dust -- wet deposition DEF
      //':Faxa_dstdry1'  &  ! flux: Size 1 dust -- dry deposition DEF
      //':Faxa_dstdry2'  &  ! flux: Size 2 dust -- dry deposition DEF
      //':Faxa_dstdry3'  &  ! flux: Size 3 dust -- dry deposition DEF
      //':Faxa_dstdry4'     ! flux: Size 4 dust -- dry deposition DEF

   !----------------------------------------------------------------------------
   ! atm fields: drv->atm
   !----------------------------------------------------------------------------
   ! DEF = default fields.  CCSM expects at least these fields to be present.
   ! Remove at your own risk.
   ! 
   ! If you add a field, add a name to the string below.  You may also
   ! document the additional field here. The relevant component must
   ! also be modified to set or access the new field.
   ! 
   !----------------------------------------------------------------------------
   ! States
   seq_flds_x2a_states = &
         'Sx_tref'     &    ! 2m reference temperature        DEF
      //':Sx_qref'     &    ! 2m reference specific humidity  DEF
      //':Sx_avsdr'    &    ! albedo, visible, direct         DEF
      //':Sx_anidr'    &    ! albedo, near-ir, direct         DEF
      //':Sx_avsdf'    &    ! albedo, visible, diffuse        DEF
      //':Sx_anidf'    &    ! albedo, near-ir, diffuse        DEF
      //':Sx_t'        &    ! surface temperature             DEF
      //':So_t'        &    ! sea surface temperature         DEF
      //':Sl_snowh'    &    ! surface snow depth              DEF
      //':Si_snowh'    &    ! surface snow depth              DEF
      //':Sx_lfrac'    &    ! surface land fraction           DEF
      //':Sx_ifrac'    &    ! surface ice fraction            DEF
      //':Sx_ofrac'    &    ! surface ocn fraction            DEF
      //':So_ustar'    &    ! needed for isoptope calc        DEF
      //':So_re'       &    ! needed for isoptope calc        DEF
      //':So_ssq'      &    ! needed for isoptope calc        DEF
      //':Sl_fv'       &    !                                 DEF
      //':Sl_ram1'     &    !                                 DEF
      //':Sx_u10'           ! 10m wind                        DEF

   ! Fluxes
   seq_flds_x2a_fluxes = &
         'Faxx_taux'   &    ! wind stress, zonal              DEF
      //':Faxx_tauy'   &    ! wind stress, meridional         DEF
      //':Faxx_lat'    &    ! latent          heat flux       DEF
      //':Faxx_sen'    &    ! sensible        heat flux       DEF
      //':Faxx_lwup'   &    ! upward longwave heat flux       DEF
      //':Faxx_evap'   &    ! evaporation    water flux       DEF
      //':Fall_flxdst1' &   ! dust flux bin 1                 DEF
      //':Fall_flxdst2' &   ! dust flux bin 2                 DEF
      //':Fall_flxdst3' &   ! dust flux bin 3                 DEF
      //':Fall_flxdst4'     ! dust flux bin 4                 DEF



   !----------------------------------------------------------------------------
   ! ice fields: ice->drv
   !----------------------------------------------------------------------------
   ! DEF = default fields.  CCSM expects at least these fields to be present.
   ! Remove at your own risk.
   ! 
   ! If you add a field, add a name to the string below.  You may also
   ! document the additional field here. The relevant component must
   ! also be modified to set or access the new field.
   ! 
   !----------------------------------------------------------------------------
   ! States
   seq_flds_i2x_states = &
         'Si_t'        &    ! temperature                     DEF
      //':Si_tref'     &    ! 2m reference temperature        DEF
      //':Si_qref'     &    ! 2m reference specific humidity  DEF
      //':Si_ifrac'    &    ! fractional ice cov wrt ocean    DEF
      //':Si_avsdr'    &    ! albedo: visible, direct         DEF
      //':Si_anidr'    &    ! albedo: near ir, direct         DEF
      //':Si_avsdf'    &    ! albedo: visible, diffuse        DEF
      //':Si_anidf'    &    ! albedo: near ir, diffuse        DEF
      //':Si_snowh'    &    ! surface snow depth (m)          DEF 
      //':Si_u10'           ! 10m wind                        DEF   

   ! Fluxes
   seq_flds_i2x_fluxes = &
         'Faii_taux'   &    ! wind stress, zonal              DEF
      //':Faii_tauy'   &    ! wind stress, meridional         DEF
      //':Faii_lat'    &    ! latent          heat flux       DEF
      //':Faii_sen'    &    ! sensible        heat flux       DEF
      //':Faii_lwup'   &    ! upward longwave heat flux       DEF
      //':Faii_evap'   &    ! evaporation    water flux       DEF
      //':Faii_swnet'  &    ! shortwave: net absorbed         DEF
      //':Fioi_swpen'  &    ! net SW penetrating ice          DEF
      //':Fioi_melth'  &    ! heat  flux from melting ice     DEF
      //':Fioi_meltw'  &    ! water flux from melting ice     DEF
      //':Fioi_salt'   &    ! salt  flux from melting ice     DEF
      //':Fioi_taux'   &    ! ice/ocn stress, zonal           DEF
      //':Fioi_tauy'        ! ice/ocn stress, meridional      DEF

   !----------------------------------------------------------------------------
   ! ice fields: drv->ice
   !----------------------------------------------------------------------------
   ! DEF = default fields.  CCSM expects at least these fields to be present.
   ! Remove at your own risk.
   ! 
   ! If you add a field, add a name to the string below.  You may also
   ! document the additional field here. The relevant component must
   ! also be modified to set or access the new field.
   !
   !----------------------------------------------------------------------------
   ! States
   seq_flds_x2i_states = &
         'So_t'        &    ! ocean layer temperature         DEF
      //':So_s'        &    ! ocn salinity                    DEF
      //':So_u'        &    ! ocn u velocity                  DEF
      //':So_v'        &    ! ocn v velocity                  DEF
      //':So_dhdx'     &    ! ocn surface slope, zonal        DEF
      //':So_dhdy'     &    ! ocn surface slope, merid        DEF
      //':Sa_u'        &    ! atm u velocity                  DEF
      //':Sa_v'        &    ! atm v velocity                  DEF
      //':Sa_z'        &    ! atm bottom layer height         DEF
      //':Sa_ptem'     &    ! atm potential temp              DEF
      //':Sa_tbot'     &    ! atm bottom temp                 DEF
      //':Sa_pbot'     &    ! atm bottom pressure             DEF
      //':Sa_shum'     &    ! atm specfic humidity            DEF
      //':Sa_dens'          ! atm air density                 DEF

   ! Fluxes
   seq_flds_x2i_fluxes = &
         'Fioo_q'      &    ! ocn freeze or melt heat         DEF
      //':Faxa_swndr'  &    ! atm sw near-ir, direct          DEF
      //':Faxa_swvdr'  &    ! atm sw visable, direct          DEF
      //':Faxa_swndf'  &    ! atm sw near-ir, diffuse         DEF
      //':Faxa_swvdf'  &    ! atm sw visable, diffuse         DEF
      //':Faxa_lwdn'   &    ! long-wave down                  DEF
      //':Faxa_rain'   &    ! prec: liquid                    DEF
      //':Faxa_snow'   &    ! prec: frozen                    DEF
      //':Faxa_bcphidry' &  ! flux: Black Carbon hydrophilic dry deposition   DEF
      //':Faxa_bcphodry' &  ! flux: Black Carbon hydrophobic dry deposition   DEF
      //':Faxa_bcphiwet' &  ! flux: Black Carbon hydrophilic wet deposition   DEF
      //':Faxa_ocphidry' &  ! flux: Organic Carbon hydrophilic dry deposition DEF
      //':Faxa_ocphodry' &  ! flux: Organic Carbon hydrophobic dry deposition DEF
      //':Faxa_ocphiwet' &  ! flux: Organic Carbon hydrophilic dry deposition DEF
      //':Faxa_dstwet1'  &  ! flux: Size 1 dust -- wet deposition DEF
      //':Faxa_dstwet2'  &  ! flux: Size 2 dust -- wet deposition DEF
      //':Faxa_dstwet3'  &  ! flux: Size 3 dust -- wet deposition DEF
      //':Faxa_dstwet4'  &  ! flux: Size 4 dust -- wet deposition DEF
      //':Faxa_dstdry1'  &  ! flux: Size 1 dust -- dry deposition DEF
      //':Faxa_dstdry2'  &  ! flux: Size 2 dust -- dry deposition DEF
      //':Faxa_dstdry3'  &  ! flux: Size 3 dust -- dry deposition DEF
      //':Faxa_dstdry4'     ! flux: Size 4 dust -- dry deposition DEF

   !----------------------------------------------------------------------------
   ! lnd fields: lnd->drv 
   !----------------------------------------------------------------------------
   ! DEF = default fields.  CCSM expects at least these fields to be present.
   ! Remove at your own risk.
   ! 
   ! If you add a field, add a name to the string below.  You may also
   ! document the additional field here. The relevant component must
   ! also be modified to set or access the new field.
   !
   !----------------------------------------------------------------------------
   ! States
    seq_flds_l2x_states = &
         'Sl_t'        &    ! temperature                     DEF
      //':Sl_tref'     &    ! 2m reference temperature        DEF
      //':Sl_qref'     &    ! 2m reference specific humidity  DEF
      //':Sl_avsdr'    &    ! albedo: direct , visible        DEF
      //':Sl_anidr'    &    ! albedo: direct , near-ir        DEF
      //':Sl_avsdf'    &    ! albedo: diffuse, visible        DEF
      //':Sl_anidf'    &    ! albedo: diffuse, near-ir        DEF
      //':Sl_snowh'    &    ! snow height                     DEF
      //':Sl_landfrac' &    ! fractional land                 DEF
      //':Sl_fv'       &    ! friction velocity               DEF
      //':Sl_ram1'     &    ! aerodynamical resistance        DEF
      //':Sl_u10'           ! 10m wind                        DEF

   ! Fluxes
    seq_flds_l2x_fluxes = &
         'Fall_taux'   &    ! wind stress, zonal              DEF
      //':Fall_tauy'   &    ! wind stress, meridional         DEF
      //':Fall_lat'    &    ! latent          heat flux       DEF
      //':Fall_sen'    &    ! sensible        heat flux       DEF
      //':Fall_lwup'   &    ! upward longwave heat flux       DEF
      //':Fall_evap'   &    ! evaporation    water flux       DEF
      //':Fall_swnet'  &    ! shortwave: net absorbed         DEF
      //':Fall_flxdst1' &   ! dust flux bin 1                 DEF
      //':Fall_flxdst2' &   ! dust flux bin 2                 DEF
      //':Fall_flxdst3' &   ! dust flux bin 3                 DEF
      //':Fall_flxdst4'     ! dust flux bin 4                 DEF

   !----------------------------------------------------------------------------
   ! lnd fields: drv->lnd
   !----------------------------------------------------------------------------
   ! DEF = default fields.  CCSM expects at least these fields to be present.
   ! Remove at your own risk.
   ! 
   ! If you add a field, add a name to the string below.  You may also
   ! document the additional field here. The relevant component must
   ! also be modified to set or access the new field.
   !
   !----------------------------------------------------------------------------
   ! States
    seq_flds_x2l_states = &
         'Sa_z'        &    ! bottom atm level height         DEF
      //':Sa_u'        &    ! bottom atm level zon wind       DEF
      //':Sa_v'        &    ! bottom atm level mer wind       DEF
      //':Sa_tbot'     &    ! bottom atm level temp           DEF
      //':Sa_ptem'     &    ! bottom atm level pot temp       DEF
      //':Sa_shum'     &    ! bottom atm level spec hum       DEF
      //':Sa_pbot'          ! bottom atm level pressure       DEF

   ! Fluxes
    seq_flds_x2l_fluxes = &
         'Faxa_lwdn'   &    ! downward longwave heat flux     DEF
      //':Faxa_rainc'  &    ! precip: liquid, convective      DEF
      //':Faxa_rainl'  &    ! precip: liquid, large-scale     DEF
      //':Faxa_snowc'  &    ! precip: frozen, convective      DEF
      //':Faxa_snowl'  &    ! precip: frozen, large-scale     DEF
      //':Faxa_swndr'  &    ! shortwave: nir direct  down     DEF
      //':Faxa_swvdr'  &    ! shortwave: vis direct  down     DEF
      //':Faxa_swndf'  &    ! shortwave: nir diffuse down     DEF
      //':Faxa_swvdf'  &    ! shortwave: vis diffuse down     DEF
      //':Faxa_bcphidry' &  ! flux: Black Carbon hydrophilic dry deposition   DEF
      //':Faxa_bcphodry' &  ! flux: Black Carbon hydrophobic dry deposition   DEF
      //':Faxa_bcphiwet' &  ! flux: Black Carbon hydrophilic wet deposition   DEF
      //':Faxa_ocphidry' &  ! flux: Organic Carbon hydrophilic dry deposition DEF
      //':Faxa_ocphodry' &  ! flux: Organic Carbon hydrophobic dry deposition DEF
      //':Faxa_ocphiwet' &  ! flux: Organic Carbon hydrophilic dry deposition DEF
      //':Faxa_dstwet1'  &  ! flux: Size 1 dust -- wet deposition DEF
      //':Faxa_dstwet2'  &  ! flux: Size 2 dust -- wet deposition DEF
      //':Faxa_dstwet3'  &  ! flux: Size 3 dust -- wet deposition DEF
      //':Faxa_dstwet4'  &  ! flux: Size 4 dust -- wet deposition DEF
      //':Faxa_dstdry1'  &  ! flux: Size 1 dust -- dry deposition DEF
      //':Faxa_dstdry2'  &  ! flux: Size 2 dust -- dry deposition DEF
      //':Faxa_dstdry3'  &  ! flux: Size 3 dust -- dry deposition DEF
      //':Faxa_dstdry4'     ! flux: Size 4 dust -- dry deposition DEF

   !----------------------------------------------------------------------------
   ! ocn fields: ocn->drv
   !----------------------------------------------------------------------------
   ! DEF = default fields.  CCSM expects at least these fields to be present.
   ! Remove at your own risk.
   ! 
   ! If you add a field, add a name to the string below.  You may also
   ! document the additional field here. The relevant component must
   ! also be modified to set or access the new field.
   !
   !----------------------------------------------------------------------------
   ! States
    seq_flds_o2x_states = &
         'So_t'        &    ! temperature                     DEF
      //':So_u'        &    ! velocity, zonal                 DEF
      //':So_v'        &    ! velocity, meridional            DEF
      //':So_s'        &    ! salinity                        DEF
      //':So_dhdx'     &    ! surface slope, zonal            DEF
      //':So_dhdy'          ! surface slope, meridional       DEF

   ! Fluxes
    seq_flds_o2x_fluxes = &
         'Fioo_q'           ! heat of fusion (q>0) melt pot (q<0)  DEF 

   !----------------------------------------------------------------------------
   ! ocn fields: drv->ocn
   !----------------------------------------------------------------------------
   ! DEF = default fields.  CCSM expects at least these fields to be present.
   ! Remove at your own risk.
   ! 
   ! If you add a field, add a name to the string below.  You may also
   ! document the additional field here. The relevant component must
   ! also be modified to set or access the new field.
   !
   !----------------------------------------------------------------------------
   ! States
    seq_flds_x2o_states = &
         'Si_ifrac'    &    ! state: ice fraction wrt ocean   DEF
      //':Sa_pslv'     &    ! state: sea level pressure       DEF
      //':Sx_duu10n'        ! state: 10m wind speed squared   DEF 

   ! Fluxes
    seq_flds_x2o_fluxes = &
         'Foxx_taux'   &    ! wind stress: zonal              DEF
      //':Foxx_tauy'   &    ! wind stress: meridional         DEF
      //':Foxx_swnet'  &    ! heat flux: shortwave net        DEF
      //':Foxx_lat'    &    ! heat flux: latent               DEF
      //':Foxx_sen'    &    ! heat flux: sensible             DEF
      //':Foxx_lwdn'   &    ! heat flux: long-wave down       DEF
      //':Foxx_lwup'   &    ! heat flux: long-wave up         DEF
      //':Foxx_melth'  &    ! heat flux: melt                 DEF
      //':Foxx_salt'   &    ! salt flux                       DEF
      //':Foxx_prec'   &    ! water flux: rain+snow           DEF
      //':Foxx_snow'   &    ! water flux: snow                DEF
      //':Foxx_rain'   &    ! water flux: rain                DEF
      //':Foxx_evap'   &    ! water flux: evap                DEF
      //':Foxx_meltw'  &    ! water flux: melt                DEF
      //':Forr_roff'   &    ! water flux: runoff              DEF
      //':Forr_ioff'   &    ! water flux: frozen runoff       DEF
      //':Foxx_bcphidry' &  ! flux: Black Carbon hydrophilic dry deposition   DEF
      //':Foxx_bcphodry' &  ! flux: Black Carbon hydrophobic dry deposition   DEF
      //':Foxx_bcphiwet' &  ! flux: Black Carbon hydrophilic wet deposition   DEF
      //':Foxx_ocphidry' &  ! flux: Organic Carbon hydrophilic dry deposition DEF
      //':Foxx_ocphodry' &  ! flux: Organic Carbon hydrophobic dry deposition DEF
      //':Foxx_ocphiwet' &  ! flux: Organic Carbon hydrophilic dry deposition DEF
      //':Foxx_dstwet1'  &  ! flux: Size 1 dust -- wet deposition DEF
      //':Foxx_dstwet2'  &  ! flux: Size 2 dust -- wet deposition DEF
      //':Foxx_dstwet3'  &  ! flux: Size 3 dust -- wet deposition DEF
      //':Foxx_dstwet4'  &  ! flux: Size 4 dust -- wet deposition DEF
      //':Foxx_dstdry1'  &  ! flux: Size 1 dust -- dry deposition DEF
      //':Foxx_dstdry2'  &  ! flux: Size 2 dust -- dry deposition DEF
      //':Foxx_dstdry3'  &  ! flux: Size 3 dust -- dry deposition DEF
      //':Foxx_dstdry4'     ! flux: Size 4 dust -- dry deposition DEF

   !----------------------------------------------------------------------------
   ! hub computed fields: atm/ocn states/fluxes
   !----------------------------------------------------------------------------
   ! States
    seq_flds_xao_albedo = &
         'So_avsdr'    &    ! albedo: visible, direct         DEF 
      //':So_anidr'    &    ! albedo: near ir, direct         DEF 
      //':So_avsdf'    &    ! albedo: visible, diffuse        DEF 
      //':So_anidf'         ! albedo: near ir, diffuse        DEF 

    seq_flds_xao_states = &
         'So_tref'     &    ! 2m reference temperature        DEF 
      //':So_qref'     &    ! 2m reference specific humidity  DEF 
      //':Sx_duu10n'   &    ! diag 10m wind speed squared     DEF
      //':So_ustar'    &    ! ustar                           DEF
      //':So_ssq'      &    ! surface saturation spec. hum.   DEF
      //':So_re'       &    ! sqrt of exch. coeff (tracers)   DEF
      //':So_u10'           ! state: 10m wind speed           DEF 
   
   ! Fluxes
    seq_flds_xao_fluxes = &
         'Faox_taux'   &    ! wind stress, zonal              DEF 
      //':Faox_tauy'   &    ! wind stress, meridional         DEF 
      //':Faox_lat'    &    ! latent          heat flux       DEF 
      //':Faox_sen'    &    ! sensible        heat flux       DEF 
      //':Faox_evap'   &    ! evaporation    water flux       DEF 
      //':Faox_lwup'        ! upward longwave heat flux       DEF 

   !----------------------------------------------------------------------------
   ! run-off field 
   !----------------------------------------------------------------------------
   ! DEF = default fields.  CCSM expects at least these fields to be present.
   ! Remove at your own risk.
   ! 
   ! If you add a field, add a name to the string below.  You may also
   ! document the additional field here. The relevant component must
   ! also be modified to set or access the new field.
   !
   !----------------------------------------------------------------------------
   ! States
    seq_flds_r2x_states = &
         ' '
   ! Fluxes
    seq_flds_r2x_fluxes = &
         'Forr_roff'       &  ! runoff to ocean              DEF
      //':Forr_ioff'          ! frozen runoff to ocean       DEF

   !----------------------------------------------------------------------------
   ! glc fields: glc->drv
   !----------------------------------------------------------------------------
   ! DEF = default fields.  CCSM expects at least these fields to be present.
   ! Remove at your own risk.
   ! 
   ! If you add a field, add a name to the string below.  You may also
   ! document the additional field here. The relevant component must
   ! also be modified to set or access the new field.
   !
   !----------------------------------------------------------------------------
   ! States
    seq_flds_g2x_states = '' &
#if (defined GLC_NEC_10 || defined GLC_NEC_5 || defined GLC_NEC_3 || defined GLC_NEC_1)
      // 'Sg_frac01'       &  ! fraction   
      //':Sg_topo01'       &  ! topography 
#endif
#if (defined GLC_NEC_10 || defined GLC_NEC_5 || defined GLC_NEC_3 )
      //':Sg_frac02'       &  !
      //':Sg_topo02'       &  !
      //':Sg_frac03'       &  !
      //':Sg_topo03'       &  !
#endif
#if (defined GLC_NEC_10 || defined GLC_NEC_5 )
      //':Sg_frac04'       &  !
      //':Sg_topo04'       &  !
      //':Sg_frac05'       &  !
      //':Sg_topo05'       &  !
#endif
#if (defined GLC_NEC_10 )
      //':Sg_frac06'       &  !
      //':Sg_topo06'       &  !
      //':Sg_frac07'       &  !
      //':Sg_topo07'       &  !
      //':Sg_frac08'       &  !
      //':Sg_topo08'       &  !
      //':Sg_frac09'       &  !
      //':Sg_topo09'       &  !
      //':Sg_frac10'       &  !
      //':Sg_topo10'       &  !
#endif
      //''                    ! empty

   ! Fluxes
    seq_flds_g2x_fluxes = '' &
#if (defined GLC_NEC_10 || defined GLC_NEC_5 || defined GLC_NEC_3 || defined GLC_NEC_1)
      // 'Fsgg_rofi01'     &   ! roff to land/sno     
      //':Fsgg_rofl01'     &   ! roff to land/sno     
      //':Fsgg_hflx01'     &   ! heat flux to land/sno
#endif
#if (defined GLC_NEC_10 || defined GLC_NEC_5 || defined GLC_NEC_3 )
      //':Fsgg_rofi02'     &   !
      //':Fsgg_rofl02'     &   !
      //':Fsgg_hflx02'     &   !
      //':Fsgg_rofi03'     &   !
      //':Fsgg_rofl03'     &   !
      //':Fsgg_hflx03'     &   !
#endif
#if (defined GLC_NEC_10 || defined GLC_NEC_5 )
      //':Fsgg_rofi04'     &   !
      //':Fsgg_rofl04'     &   !
      //':Fsgg_hflx04'     &   !
      //':Fsgg_rofi05'     &   !
      //':Fsgg_rofl05'     &   !
      //':Fsgg_hflx05'     &   !
#endif
#if (defined GLC_NEC_10 )
      //':Fsgg_rofi06'     &   !
      //':Fsgg_rofl06'     &   !
      //':Fsgg_hflx06'     &   !
      //':Fsgg_rofi07'     &   !
      //':Fsgg_rofl07'     &   !
      //':Fsgg_hflx07'     &   !
      //':Fsgg_rofi08'     &   !
      //':Fsgg_rofl08'     &   !
      //':Fsgg_hflx08'     &   !
      //':Fsgg_rofi09'     &   !
      //':Fsgg_rofl09'     &   !
      //':Fsgg_hflx09'     &   !
      //':Fsgg_rofi10'     &   !
      //':Fsgg_rofl10'     &   !
      //':Fsgg_hflx10'     &   !
#endif
      //''                     ! empty

   !----------------------------------------------------------------------------
   ! glc fields: drv->glc
   !----------------------------------------------------------------------------
   ! DEF = default fields.  CCSM expects at least these fields to be present.
   ! Remove at your own risk.
   ! 
   ! If you add a field, add a name to the string below.  You may also
   ! document the additional field here. The relevant component must
   ! also be modified to set or access the new field.
   !
   !----------------------------------------------------------------------------
   ! States
    seq_flds_x2g_states = '' &
#if (defined GLC_NEC_10 || defined GLC_NEC_5 || defined GLC_NEC_3 || defined GLC_NEC_1)
      // 'Ss_tsrf01'       &  ! surface temp 
      //':Ss_topo01'       &  ! topography
#endif
#if (defined GLC_NEC_10 || defined GLC_NEC_5 || defined GLC_NEC_3 )
      //':Ss_tsrf02'       &  ! 
      //':Ss_topo02'       &  ! 
      //':Ss_tsrf03'       &  ! 
      //':Ss_topo03'       &  ! 
#endif
#if (defined GLC_NEC_10 || defined GLC_NEC_5 )
      //':Ss_tsrf04'       &  ! 
      //':Ss_topo04'       &  ! 
      //':Ss_tsrf05'       &  ! 
      //':Ss_topo05'       &  ! 
#endif
#if (defined GLC_NEC_10 )
      //':Ss_tsrf06'       &  ! 
      //':Ss_topo06'       &  ! 
      //':Ss_tsrf07'       &  ! 
      //':Ss_topo07'       &  ! 
      //':Ss_tsrf08'       &  ! 
      //':Ss_topo08'       &  ! 
      //':Ss_tsrf09'       &  ! 
      //':Ss_topo09'       &  ! 
      //':Ss_tsrf10'       &  ! 
      //':Ss_topo10'       &  ! 
#endif
      //''                     ! empty

   ! Fluxes
    seq_flds_x2g_fluxes = '' &
#if (defined GLC_NEC_10 || defined GLC_NEC_5 || defined GLC_NEC_3 || defined GLC_NEC_1)
      // 'Fgss_qice01'     &   ! land/sno qice to glc
#endif
#if (defined GLC_NEC_10 || defined GLC_NEC_5 || defined GLC_NEC_3 )
      //':Fgss_qice02'     &   !
      //':Fgss_qice03'     &   !
#endif
#if (defined GLC_NEC_10 || defined GLC_NEC_5 )
      //':Fgss_qice04'     &   !
      //':Fgss_qice05'     &   !
#endif
#if (defined GLC_NEC_10 )
      //':Fgss_qice06'     &   !
      //':Fgss_qice07'     &   !
      //':Fgss_qice08'     &   !
      //':Fgss_qice09'     &   !
      //':Fgss_qice10'     &   !
#endif
      //''                     ! empty

   !----------------------------------------------------------------------------
   ! sno fields: sno->drv
   !----------------------------------------------------------------------------
    seq_flds_s2x_states = seq_flds_x2g_states
    seq_flds_s2x_fluxes = seq_flds_x2g_fluxes

   !----------------------------------------------------------------------------
   ! sno fields: drv->sno
   !----------------------------------------------------------------------------
    seq_flds_x2s_states = seq_flds_g2x_states
    seq_flds_x2s_fluxes = seq_flds_g2x_fluxes

    !----------------------------------------------------------------------------
    ! optional fields
    !----------------------------------------------------------------------------

#if (defined CO2A) 
    seq_flds_a2x_states   = trim(seq_flds_a2x_states) &
        //':Sa_co2prog'    &  ! bottom atm level prognostic co2
        //':Sa_co2diag'       ! bottom atm level diagnostic co2

    seq_flds_x2l_states   = trim(seq_flds_x2l_states) &
        //':Sa_co2prog'    &  ! atm prognostic co2
        //':Sa_co2diag'       ! atm diagnostic co2

#elif (defined CO2B) 
   seq_flds_a2x_states   = trim(seq_flds_a2x_states) &
        //':Sa_co2prog'    &  ! bottom atm level prognostic co2
        //':Sa_co2diag'       ! bottom atm level diagnostic co2

   seq_flds_x2a_fluxes   = trim(seq_flds_x2a_fluxes) &
        //':Faxx_fco2_lnd'    ! co2 flux from lnd

   seq_flds_l2x_fluxes   = trim(seq_flds_l2x_fluxes) &
        //':Fall_fco2_lnd'    ! co2 flux from lnd

   seq_flds_x2l_states   = trim(seq_flds_x2l_states) &
        //':Sa_co2prog'    &  ! atm prognostic co2
        //':Sa_co2diag'       ! atm diagnostic co2

#elif (defined CO2C) 
   seq_flds_a2x_states   = trim(seq_flds_a2x_states) &
        //':Sa_co2prog'    &  ! bottom atm level prognostic co2
        //':Sa_co2diag'       ! bottom atm level diagnostic co2

   seq_flds_x2a_fluxes   = trim(seq_flds_x2a_fluxes) &
        //':Faxx_fco2_lnd' &  ! co2 flux from lnd
        //':Faxx_fco2_ocn'    ! co2 flux from ocn

   seq_flds_l2x_fluxes   = trim(seq_flds_l2x_fluxes) &
        //':Fall_fco2_lnd'    ! co2 flux from lnd

   seq_flds_x2l_states   = trim(seq_flds_x2l_states) &
        //':Sa_co2prog'    &  ! atm prognostic co2
        //':Sa_co2diag'       ! atm diagnostic co2

   seq_flds_o2x_fluxes   = trim(seq_flds_o2x_fluxes) &
        //':Faoo_fco2'        ! co2 flux                       

   seq_flds_x2o_states   = trim(seq_flds_x2o_states) &
        //':Sa_co2prog'    &  ! bottom atm level prognostic co2
        //':Sa_co2diag'       ! atm diagnostic co2

#elif (defined CO2_DMSA)
   seq_flds_a2x_states   = trim(seq_flds_a2x_states) &
        //':Sa_co2prog'    &  ! bottom atm level prognostic co2
        //':Sa_co2diag'       ! bottom atm level diagnostic co2

   seq_flds_x2a_fluxes   = trim(seq_flds_x2a_fluxes) &
        //':Faxx_fco2_lnd' &  ! co2 flux from lnd
        //':Faxx_fco2_ocn' &  ! co2 flux from ocn
        //':Faxx_fdms'        ! dms flux

   seq_flds_l2x_fluxes   = trim(seq_flds_l2x_fluxes) &
        //':Fall_fco2_lnd'    ! co2 flux from lnd

   seq_flds_x2l_states   = trim(seq_flds_x2l_states) &
        //':Sa_co2prog'    &  ! atm prognostic co2
        //':Sa_co2diag'       ! atm diagnostic co2

   seq_flds_o2x_fluxes   = trim(seq_flds_o2x_fluxes) &
        //':Faoo_fco2'     &  ! co2 flux                       
        //':Faoo_fdms'        ! dms flux

   seq_flds_x2o_states   = trim(seq_flds_x2o_states) &
        //':Sa_co2prog'    &  ! bottom atm level prognostic co2
        //':Sa_co2diag'       ! atm diagnostic co2
#endif

#if (defined VOC )
    seq_flds_x2a_fluxes = trim(seq_flds_x2a_fluxes) &
         //':Fall_flxvoc1' &   ! VOC flux bin 1
         //':Fall_flxvoc2' &   ! VOC flux bin 2
         //':Fall_flxvoc3' &   ! VOC flux bin 3
         //':Fall_flxvoc4' &   ! VOC flux bin 4
         //':Fall_flxvoc5'     ! VOC flux bin 5

    seq_flds_l2x_fluxes = trim(seq_flds_l2x_fluxes) &
         //':Fall_flxvoc1' &   ! VOC flux bin 1
         //':Fall_flxvoc2' &   ! VOC flux bin 2
         //':Fall_flxvoc3' &   ! VOC flux bin 3
         //':Fall_flxvoc4' &   ! VOC flux bin 4
         //':Fall_flxvoc5'     ! VOC flux bin 5
#endif


    !-----------------------------------------------------------------------------
    ! Dry Deposition initialization
    ! First read namelist and figure out the drydep field list to pass
    ! Then check if file exists and if not, n_drydep will be zero
    ! Then add dry deposition fields to land export and atmosphere import states
    ! Then initialize dry deposition fields
    ! Note: CAM and CLM will then call seq_drydep_setHCoeff
    !-----------------------------------------------------------------------------
 
    call seq_drydep_read(nlfilename='drv_flds_in', seq_drydep_fields=seq_drydep_fields)

    if ( lnd_drydep ) then
       seq_flds_x2a_states = trim(seq_flds_x2a_states) // trim(seq_drydep_fields)
       seq_flds_l2x_states = trim(seq_flds_l2x_states) // trim(seq_drydep_fields)
    endif

    call seq_drydep_init( )

    !----------------------------------------------------------------------------
    ! state + flux fields
    !----------------------------------------------------------------------------

    !validate lengths of states strings
    if (len_trim(seq_flds_a2x_states ) .ge. len(seq_flds_a2x_states )) &
        call shr_sys_abort('maximum length of string has been exceeded for: seq_flds_a2x_states ')
    if (len_trim(seq_flds_x2a_states ) .ge. len(seq_flds_x2a_states )) &
        call shr_sys_abort('maximum length of string has been exceeded for: seq_flds_x2a_states ')
    if (len_trim(seq_flds_i2x_states ) .ge. len(seq_flds_i2x_states )) &
        call shr_sys_abort('maximum length of string has been exceeded for: seq_flds_i2x_states ')
    if (len_trim(seq_flds_x2i_states ) .ge. len(seq_flds_x2i_states )) &
        call shr_sys_abort('maximum length of string has been exceeded for: seq_flds_x2i_states ')
    if (len_trim(seq_flds_l2x_states ) .ge. len(seq_flds_l2x_states )) &
        call shr_sys_abort('maximum length of string has been exceeded for: seq_flds_l2x_states ')
    if (len_trim(seq_flds_x2l_states ) .ge. len(seq_flds_x2l_states )) &
        call shr_sys_abort('maximum length of string has been exceeded for: seq_flds_x2l_states ')
    if (len_trim(seq_flds_o2x_states ) .ge. len(seq_flds_o2x_states )) &
        call shr_sys_abort('maximum length of string has been exceeded for: seq_flds_o2x_states ')
    if (len_trim(seq_flds_x2o_states ) .ge. len(seq_flds_x2o_states )) &
        call shr_sys_abort('maximum length of string has been exceeded for: seq_flds_x2o_states ')
    if (len_trim(seq_flds_g2x_states ) .ge. len(seq_flds_g2x_states )) &
        call shr_sys_abort('maximum length of string has been exceeded for: seq_flds_g2x_states ')
    if (len_trim(seq_flds_x2g_states ) .ge. len(seq_flds_x2g_states )) &
        call shr_sys_abort('maximum length of string has been exceeded for: seq_flds_x2g_states ')
    if (len_trim(seq_flds_s2x_states ) .ge. len(seq_flds_s2x_states )) &
        call shr_sys_abort('maximum length of string has been exceeded for: seq_flds_s2x_states ')
    if (len_trim(seq_flds_x2s_states ) .ge. len(seq_flds_x2s_states )) &
        call shr_sys_abort('maximum length of string has been exceeded for: seq_flds_x2s_states ')
    if (len_trim(seq_flds_xao_albedo ) .ge. len(seq_flds_xao_albedo )) &
        call shr_sys_abort('maximum length of string has been exceeded for: seq_flds_xao_albedo ')
    if (len_trim(seq_flds_xao_states ) .ge. len(seq_flds_xao_states )) &
        call shr_sys_abort('maximum length of string has been exceeded for: seq_flds_xao_states ')
    if (len_trim(seq_flds_r2x_states ) .ge. len(seq_flds_r2x_states )) &
        call shr_sys_abort('maximum length of string has been exceeded for: seq_flds_r2x_states ')


    !validate lengths of fluxes strings
    if (len_trim(seq_flds_a2x_fluxes ) .ge. len(seq_flds_a2x_fluxes )) &
        call shr_sys_abort('maximum length of string has been exceeded for: seq_flds_a2x_fluxes ')
    if (len_trim(seq_flds_x2a_fluxes ) .ge. len(seq_flds_x2a_fluxes )) &
        call shr_sys_abort('maximum length of string has been exceeded for: seq_flds_x2a_fluxes ')
    if (len_trim(seq_flds_i2x_fluxes ) .ge. len(seq_flds_i2x_fluxes )) &
        call shr_sys_abort('maximum length of string has been exceeded for: seq_flds_i2x_fluxes ')
    if (len_trim(seq_flds_x2i_fluxes ) .ge. len(seq_flds_x2i_fluxes )) &
        call shr_sys_abort('maximum length of string has been exceeded for: seq_flds_x2i_fluxes ')
    if (len_trim(seq_flds_l2x_fluxes ) .ge. len(seq_flds_l2x_fluxes )) &
        call shr_sys_abort('maximum length of string has been exceeded for: seq_flds_l2x_fluxes ')
    if (len_trim(seq_flds_x2l_fluxes ) .ge. len(seq_flds_x2l_fluxes )) &
        call shr_sys_abort('maximum length of string has been exceeded for: seq_flds_x2l_fluxes ')
    if (len_trim(seq_flds_o2x_fluxes ) .ge. len(seq_flds_o2x_fluxes )) &
        call shr_sys_abort('maximum length of string has been exceeded for: seq_flds_o2x_fluxes ')
    if (len_trim(seq_flds_x2o_fluxes ) .ge. len(seq_flds_x2o_fluxes )) &
        call shr_sys_abort('maximum length of string has been exceeded for: seq_flds_x2o_fluxes ')
    if (len_trim(seq_flds_g2x_fluxes ) .ge. len(seq_flds_g2x_fluxes )) &
        call shr_sys_abort('maximum length of string has been exceeded for: seq_flds_g2x_fluxes ')
    if (len_trim(seq_flds_x2g_fluxes ) .ge. len(seq_flds_x2g_fluxes )) &
        call shr_sys_abort('maximum length of string has been exceeded for: seq_flds_x2g_fluxes ')
    if (len_trim(seq_flds_s2x_fluxes ) .ge. len(seq_flds_s2x_fluxes )) &
        call shr_sys_abort('maximum length of string has been exceeded for: seq_flds_s2x_fluxes ')
    if (len_trim(seq_flds_x2s_fluxes ) .ge. len(seq_flds_x2s_fluxes )) &
        call shr_sys_abort('maximum length of string has been exceeded for: seq_flds_x2s_fluxes ')
    if (len_trim(seq_flds_xao_fluxes ) .ge. len(seq_flds_xao_fluxes )) &
        call shr_sys_abort('maximum length of string has been exceeded for: seq_flds_xao_fluxes ')
    if (len_trim(seq_flds_r2x_fluxes ) .ge. len(seq_flds_r2x_fluxes )) &
        call shr_sys_abort('maximum length of string has been exceeded for: seq_flds_r2x_fluxes ')


    !validate lengths of fields strings
    if ((len_trim(seq_flds_dom_coord)  + len_trim(seq_flds_dom_other))  .ge. len(seq_flds_dom_fields)) &
        call shr_sys_abort('maximum length of string has been exceeded for: seq_flds_dom_fields')
    if ((len_trim(seq_flds_a2x_states) + len_trim(seq_flds_a2x_fluxes)) .ge. len(seq_flds_a2x_fields)) &
        call shr_sys_abort('maximum length of string has been exceeded for: seq_flds_a2x_fields')
    if ((len_trim(seq_flds_x2a_states) + len_trim(seq_flds_x2a_fluxes)) .ge. len(seq_flds_x2a_fields)) &
        call shr_sys_abort('maximum length of string has been exceeded for: seq_flds_x2a_fields')
    if ((len_trim(seq_flds_i2x_states) + len_trim(seq_flds_i2x_fluxes)) .ge. len(seq_flds_i2x_fields)) &
        call shr_sys_abort('maximum length of string has been exceeded for: seq_flds_i2x_fields')
    if ((len_trim(seq_flds_x2i_states) + len_trim(seq_flds_x2i_fluxes)) .ge. len(seq_flds_x2i_fields)) &
        call shr_sys_abort('maximum length of string has been exceeded for: seq_flds_x2i_fields')
    if ((len_trim(seq_flds_l2x_states) + len_trim(seq_flds_l2x_fluxes)) .ge. len(seq_flds_l2x_fields)) &
        call shr_sys_abort('maximum length of string has been exceeded for: seq_flds_l2x_fields')
    if ((len_trim(seq_flds_x2l_states) + len_trim(seq_flds_x2l_fluxes)) .ge. len(seq_flds_x2l_fields)) &
        call shr_sys_abort('maximum length of string has been exceeded for: seq_flds_x2l_fields')
    if ((len_trim(seq_flds_o2x_states) + len_trim(seq_flds_o2x_fluxes)) .ge. len(seq_flds_o2x_fields)) &
        call shr_sys_abort('maximum length of string has been exceeded for: seq_flds_o2x_fields')
    if ((len_trim(seq_flds_x2o_states) + len_trim(seq_flds_x2o_fluxes)) .ge. len(seq_flds_x2o_fields)) &
        call shr_sys_abort('maximum length of string has been exceeded for: seq_flds_x2o_fields')
    if ((len_trim(seq_flds_g2x_states) + len_trim(seq_flds_g2x_fluxes)) .ge. len(seq_flds_g2x_fields)) &
        call shr_sys_abort('maximum length of string has been exceeded for: seq_flds_g2x_fields')
    if ((len_trim(seq_flds_x2g_states) + len_trim(seq_flds_x2g_fluxes)) .ge. len(seq_flds_x2g_fields)) &
        call shr_sys_abort('maximum length of string has been exceeded for: seq_flds_x2g_fields')
    if ((len_trim(seq_flds_s2x_states) + len_trim(seq_flds_s2x_fluxes)) .ge. len(seq_flds_s2x_fields)) &
        call shr_sys_abort('maximum length of string has been exceeded for: seq_flds_s2x_fields')
    if ((len_trim(seq_flds_x2s_states) + len_trim(seq_flds_x2s_fluxes)) .ge. len(seq_flds_x2s_fields)) &
        call shr_sys_abort('maximum length of string has been exceeded for: seq_flds_x2s_fields')
    if ((len_trim(seq_flds_xao_albedo) + len_trim(seq_flds_xao_states) + len_trim(seq_flds_xao_fluxes)) .ge. len(seq_flds_xao_fields)) &
        call shr_sys_abort('maximum length of string has been exceeded for: seq_flds_xao_fields')
    if ((len_trim(seq_flds_r2x_states) + len_trim(seq_flds_r2x_fluxes)) .ge. len(seq_flds_r2x_fields)) &
        call shr_sys_abort('maximum length of string has been exceeded for: seq_flds_r2x_fields')


!    seq_flds_dom_fields = trim(seq_flds_dom_coord)//":"//trim(seq_flds_dom_other)
!    seq_flds_a2x_fields = trim(seq_flds_a2x_states)//":"//trim(seq_flds_a2x_fluxes)
!    seq_flds_x2a_fields = trim(seq_flds_x2a_states)//":"//trim(seq_flds_x2a_fluxes)
!    seq_flds_i2x_fields = trim(seq_flds_i2x_states)//":"//trim(seq_flds_i2x_fluxes)
!    seq_flds_x2i_fields = trim(seq_flds_x2i_states)//":"//trim(seq_flds_x2i_fluxes)
!    seq_flds_l2x_fields = trim(seq_flds_l2x_states)//":"//trim(seq_flds_l2x_fluxes)
!    seq_flds_x2l_fields = trim(seq_flds_x2l_states)//":"//trim(seq_flds_x2l_fluxes)
!    seq_flds_o2x_fields = trim(seq_flds_o2x_states)//":"//trim(seq_flds_o2x_fluxes)
!    seq_flds_x2o_fields = trim(seq_flds_x2o_states)//":"//trim(seq_flds_x2o_fluxes)
!    seq_flds_g2x_fields = trim(seq_flds_g2x_states)//":"//trim(seq_flds_g2x_fluxes)
!    seq_flds_x2g_fields = trim(seq_flds_x2g_states)//":"//trim(seq_flds_x2g_fluxes)
!    seq_flds_s2x_fields = trim(seq_flds_s2x_states)//":"//trim(seq_flds_s2x_fluxes)
!    seq_flds_x2s_fields = trim(seq_flds_x2s_states)//":"//trim(seq_flds_x2s_fluxes)
!    seq_flds_xao_fields = trim(seq_flds_xao_states)//":"//trim(seq_flds_xao_fluxes)
!    seq_flds_r2x_fields = trim(seq_flds_r2x_fluxes)
!!   seq_flds_r2x_fields = trim(seq_flds_r2x_states)//":"//trim(seq_flds_r2x_fluxes)

    call seq_flds_catFields(seq_flds_dom_fields,seq_flds_dom_coord ,seq_flds_dom_other )
    call seq_flds_catFields(seq_flds_a2x_fields,seq_flds_a2x_states,seq_flds_a2x_fluxes)
    call seq_flds_catFields(seq_flds_x2a_fields,seq_flds_x2a_states,seq_flds_x2a_fluxes)
    call seq_flds_catFields(seq_flds_i2x_fields,seq_flds_i2x_states,seq_flds_i2x_fluxes)
    call seq_flds_catFields(seq_flds_x2i_fields,seq_flds_x2i_states,seq_flds_x2i_fluxes)
    call seq_flds_catFields(seq_flds_l2x_fields,seq_flds_l2x_states,seq_flds_l2x_fluxes)
    call seq_flds_catFields(seq_flds_x2l_fields,seq_flds_x2l_states,seq_flds_x2l_fluxes)
    call seq_flds_catFields(seq_flds_o2x_fields,seq_flds_o2x_states,seq_flds_o2x_fluxes)
    call seq_flds_catFields(seq_flds_x2o_fields,seq_flds_x2o_states,seq_flds_x2o_fluxes)
    call seq_flds_catFields(seq_flds_g2x_fields,seq_flds_g2x_states,seq_flds_g2x_fluxes)
    call seq_flds_catFields(seq_flds_x2g_fields,seq_flds_x2g_states,seq_flds_x2g_fluxes)
    call seq_flds_catFields(seq_flds_s2x_fields,seq_flds_s2x_states,seq_flds_s2x_fluxes)
    call seq_flds_catFields(seq_flds_x2s_fields,seq_flds_x2s_states,seq_flds_x2s_fluxes)
    call seq_flds_catFields(stringtmp          ,seq_flds_xao_albedo,seq_flds_xao_states)
    call seq_flds_catFields(seq_flds_xao_fields,stringtmp          ,seq_flds_xao_fluxes)
    call seq_flds_catFields(seq_flds_r2x_fields,''                 ,seq_flds_r2x_fluxes)
!   call seq_flds_catFields(seq_flds_r2x_fields,seq_flds_r2x_states,seq_flds_r2x_fluxes)
    
end subroutine seq_flds_set

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: seq_flds_catFields
!
! !DESCRIPTION:
!  Returns {\tt nfld} concatentated field lists
!  in the output character string {\tt outfield}.
!
! !REVISION HISTORY:
!  2003-Jan-24  - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine seq_flds_catFields(outfield, str1, str2)

! !USES:

! !INPUT/OUTPUT PARAMETERS:

   character(len=*),intent(inout) :: outfield   ! output field name
   character(len=*),intent(in)    :: str1       ! string1 
   character(len=*),intent(in )   :: str2       ! string2

!EOP

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

  outfield = ' '

  if (len_trim(str1) > 0 .and. len_trim(str2) > 0) then
     if (len_trim(str1) + len_trim(str2) > len(outfield)) then
        call shr_sys_abort('seq_flds_catFields: maximum length of string has been exceeded sum')
     endif
     outfield = trim(str1)//':'//trim(str2)
  else
     if (len_trim(str1) > 0) then
        if (len_trim(str1) > len(outfield)) then
           call shr_sys_abort('seq_flds_catFields: maximum length of string has been exceeded str1')
        endif
        outfield = trim(str1)
     endif
     if (len_trim(str2) > 0) then
        if (len_trim(str2) > len(outfield)) then
           call shr_sys_abort('seq_flds_catFields: maximum length of string has been exceeded str2')
        endif
        outfield = trim(str2)
     endif
  endif

end subroutine seq_flds_catFields

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: seq_flds_getField
!
! !DESCRIPTION:
!  Returns {\tt nfld} element of the colon-delimited string {\tt cstring}
!  in the output character string {\tt outfield}.
!
! !REVISION HISTORY:
!  2003-Jan-24  - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine seq_flds_getField(outfield, nfld, cstring)

! !USES:
   use mct_mod

! !INPUT/OUTPUT PARAMETERS:

   character(len=*),intent(out) :: outfield   ! output field name
   integer         ,intent(in ) :: nfld       ! field number
   character(len=*),intent(in ) :: cstring    ! colon delimited field string

!EOP

  type(mct_list)   :: mctIstr  ! mct list from input cstring
  type(mct_string) :: mctOStr  ! mct string for output outfield

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

  outfield = ' '

  call mct_list_init(mctIstr,cstring)
  call mct_list_get(mctOStr,nfld,mctIstr)
  outfield = mct_string_toChar(mctOStr)
  call mct_list_clean(mctIstr)
  call mct_string_clean(mctOStr)

end subroutine seq_flds_getField

!===============================================================================

subroutine seq_flds_esmf_metadata_get(shortname, longname, stdname, units)

! !USES:

   use shr_string_mod, only : shr_string_lastindex

   implicit none

! !INPUT/OUTPUT PARAMETERS:

  character(len=*), intent(in)  :: shortname 
  character(len=*),optional, intent(out) :: longname
  character(len=*),optional, intent(out) :: stdname 
  character(len=*),optional, intent(out) :: units 

!EOP

  !--- local ---
  integer :: i,n
  integer, parameter :: nmax = 246 
  logical, save :: firstCall = .true.
  character(len=80), dimension(nmax,4), save :: lookup
  character(len=80) :: llongname, lstdname, lunits, lshortname  ! local copies
  character(len=*),parameter :: undef = 'undefined'
  character(len=*),parameter :: unknown = 'unknown'
  logical :: found
  
  !--- define field metadata (name, long_name, standard_name, units) ---
  if (firstCall) then
    lookup(:,:) = trim(undef)
    
    !----------------------------------------------------------------------------

    lookup(1,1) = 'Sa_z'
    lookup(1,2) = 'Height at the lowest model level'
    lookup(1,3) = 'height'
    lookup(1,4) = 'm'

    lookup(2,1) = 'Sa_u'
    lookup(2,2) = 'Zonal wind at the lowest model level'
    lookup(2,3) = 'eastward_wind'
    lookup(2,4) = 'm s-1'

    lookup(3,1) = 'Sa_v'
    lookup(3,2) = 'Meridional wind at the lowest model level'
    lookup(3,3) = 'northward_wind'
    lookup(3,4) = 'm s-1'

    lookup(4,1) = 'Sa_tbot'
    lookup(4,2) = 'Temperature at the lowest model level'
    lookup(4,3) = 'air_temperature'
    lookup(4,4) = 'K'

    lookup(5,1) = 'Sa_ptem'
    lookup(5,2) = 'Potential temperature at the lowest model level'
    lookup(5,3) = 'air_potential_temperature'
    lookup(5,4) = 'K'

    lookup(6,1) = 'Sa_shum'
    lookup(6,2) = 'Specific humidity at the lowest model level'
    lookup(6,3) = 'specific_humidity'
    lookup(6,4) = 'kg kg-1'

    lookup(7,1) = 'Sa_dens'
    lookup(7,2) = 'Air density at the lowest model level'
    lookup(7,3) = 'air_density'
    lookup(7,4) = 'kg m-3'

    lookup(8,1) = 'Sa_pbot'
    lookup(8,2) = 'Pressure at the lowest model level'
    lookup(8,3) = 'air_pressure'
    lookup(8,4) = 'Pa'

    lookup(9,1) = 'Sa_pslv'
    lookup(9,2) = 'Sea level pressure'
    lookup(9,3) = 'air_pressure_at_sea_level'
    lookup(9,4) = 'Pa'

    lookup(10,1) = 'Faxa_lwdn'
    lookup(10,2) = 'Downward longwave heat flux'
    lookup(10,3) = 'downwelling_longwave_flux'
    lookup(10,4) = 'W m-2'

    lookup(11,1) = 'Faxa_rainc'
    lookup(11,2) = 'Convective precipitation rate'
    lookup(11,3) = 'convective_precipitation_flux'
    lookup(11,4) = 'kg m-2 s-1'

    lookup(12,1) = 'Faxa_rainl'
    lookup(12,2) = 'Large-scale (stable) precipitation rate'
    lookup(12,3) = 'large_scale_precipitation_flux'
    lookup(12,4) = 'kg m-2 s-1'

    lookup(13,1) = 'Faxa_snowc'
    lookup(13,2) = 'Convective snow rate (water equivalent)'
    lookup(13,3) = 'convective_snowfall_flux'
    lookup(13,4) = 'kg m-2 s-1'

    lookup(14,1) = 'Faxa_snowl'
    lookup(14,2) = 'Large-scale (stable) snow rate (water equivalent)'
    lookup(14,3) = 'large_scale_snowfall_flux'
    lookup(14,4) = 'kg m-2 s-1'

    lookup(15,1) = 'Faxa_swndr'
    lookup(15,2) = 'Direct near-infrared incident solar radiation'
    lookup(15,3) = 'surface_downward_direct_shortwave_flux_due_to_near_infrared_radiation'
    lookup(15,4) = 'W m-2'

    lookup(16,1) = 'Faxa_swvdr'
    lookup(16,2) = 'Direct visible incident solar radiation'
    lookup(16,3) = 'surface_downward_direct_shortwave_flux_due_to_visible_radiation'
    lookup(16,4) = 'W m-2'

    lookup(17,1) = 'Faxa_swndf'
    lookup(17,2) = 'Diffuse near-infrared incident solar radiation'
    lookup(17,3) = 'surface_downward_diffuse_shortwave_flux_due_to_near_infrared_radiation'
    lookup(17,4) = 'W m-2'

    lookup(18,1) = 'Faxa_swvdf'
    lookup(18,2) = 'Diffuse visible incident solar radiation'
    lookup(18,3) = 'surface_downward_diffuse_shortwave_flux_due_to_visible_radiation'
    lookup(18,4) = 'W m-2'

    lookup(19,1) = 'Faxa_swnet'
    lookup(19,2) = 'Net shortwave radiation'
    lookup(19,3) = 'surface_net_shortwave_flux'
    lookup(19,4) = 'W m-2'

    lookup(20,1) = 'Faxa_bcphidry'
    lookup(20,2) = 'Hydrophylic black carbon dry deposition flux'
    lookup(20,3) = 'dry_deposition_flux_of_hydrophylic_black_carbon'
    lookup(20,4) = 'kg m-2 s-1'

    lookup(21,1) = 'Faxa_bcphodry'
    lookup(21,2) = 'Hydrophobic black carbon dry deposition flux'
    lookup(21,3) = 'dry_deposition_flux_of_hydrophobic_black_carbon'
    lookup(21,4) = 'kg m-2 s-1'

    lookup(22,1) = 'Faxa_bcphiwet'
    lookup(22,2) = 'Hydrophylic black carbon wet deposition flux'
    lookup(22,3) = 'wet_deposition_flux_of_hydrophylic_black_carbon'
    lookup(22,4) = 'kg m-2 s-1'

    lookup(23,1) = 'Faxa_ocphidry'
    lookup(23,2) = 'Hydrophylic organic carbon dry deposition flux'
    lookup(23,3) = 'dry_deposition_flux_of_hydrophylic_organic_carbon'
    lookup(23,4) = 'kg m-2 s-1'

    lookup(24,1) = 'Faxa_ocphodry'
    lookup(24,2) = 'Hydrophobic organic carbon dry deposition flux'
    lookup(24,3) = 'dry_deposition_flux_of_hydrophobic_organic_carbon'
    lookup(24,4) = 'kg m-2 s-1'

    lookup(25,1) = 'Faxa_ocphiwet'
    lookup(25,2) = 'Hydrophylic organic carbon wet deposition flux'
    lookup(25,3) = 'wet_deposition_flux_of_hydrophylic_organic_carbon'
    lookup(25,4) = 'kg m-2 s-1'

    lookup(26,1) = 'Faxa_dstwet1'
    lookup(26,2) = 'Dust wet deposition flux (size 1)'
    lookup(26,3) = 'wet_deposition_flux_of_dust'
    lookup(26,4) = 'kg m-2 s-1'

    lookup(27,1) = 'Faxa_dstwet2'
    lookup(27,2) = 'Dust wet deposition flux (size 2)'
    lookup(27,3) = 'wet_deposition_flux_of_dust'
    lookup(27,4) = 'kg m-2 s-1'

    lookup(28,1) = 'Faxa_dstwet3'
    lookup(28,2) = 'Dust wet deposition flux (size 3)'
    lookup(28,3) = 'wet_deposition_flux_of_dust'
    lookup(28,4) = 'kg m-2 s-1'

    lookup(29,1) = 'Faxa_dstwet4'
    lookup(29,2) = 'Dust wet deposition flux (size 4)'
    lookup(29,3) = 'wet_deposition_flux_of_dust'
    lookup(29,4) = 'kg m-2 s-1'

    lookup(30,1) = 'Faxa_dstdry1'
    lookup(30,2) = 'Dust dry deposition flux (size 1)'
    lookup(30,3) = 'dry_deposition_flux_of_dust'
    lookup(30,4) = 'kg m-2 s-1'

    lookup(31,1) = 'Faxa_dstdry2'
    lookup(31,2) = 'Dust dry deposition flux (size 2)'
    lookup(31,3) = 'dry_deposition_flux_of_dust'
    lookup(31,4) = 'kg m-2 s-1'

    lookup(32,1) = 'Faxa_dstdry3'
    lookup(32,2) = 'Dust dry deposition flux (size 3)'
    lookup(32,3) = 'dry_deposition_flux_of_dust'
    lookup(32,4) = 'kg m-2 s-1'

    lookup(33,1) = 'Faxa_dstdry4'
    lookup(33,2) = 'Dust dry deposition flux (size 4)'
    lookup(33,3) = 'dry_deposition_flux_of_dust'
    lookup(33,4) = 'kg m-2 s-1'

    lookup(34,1) = 'Sx_tref'
    lookup(34,2) = 'Reference temperature at 2 meters'
    lookup(34,3) = 'air_temperature'
    lookup(34,4) = 'K'

    lookup(35,1) = 'Sx_qref'
    lookup(35,2) = 'Reference specific humidity at 2 meters'
    lookup(35,3) = 'specific_humidity'
    lookup(35,4) = 'kg kg-1'

    lookup(36,1) = 'Sx_avsdr'
    lookup(36,2) = 'Direct albedo (visible radiation)'
    lookup(36,3) = 'surface_direct_albedo_due_to_visible_radiation'
    lookup(36,4) = 'unitless'

    lookup(37,1) = 'Sx_anidr'
    lookup(37,2) = 'Direct albedo (near-infrared radiation)'
    lookup(37,3) = 'surface_direct_albedo_due_to_near_infrared_radiation'
    lookup(37,4) = 'unitless'

    lookup(38,1) = 'Sx_avsdf'
    lookup(38,2) = 'Diffuse albedo (visible radiation)'
    lookup(38,3) = 'surface_diffuse_albedo_due_to_visible_radiation'
    lookup(38,4) = 'unitless'

    lookup(39,1) = 'Sx_anidf'
    lookup(39,2) = 'Diffuse albedo (near-infrared radiation)'
    lookup(39,3) = 'surface_diffuse_albedo_due_to_near_infrared_radiation'
    lookup(39,4) = 'unitless'

    lookup(40,1) = 'Sx_t'
    lookup(40,2) = 'Surface temperature'
    lookup(40,3) = 'surface_temperature'
    lookup(40,4) = 'K'

    lookup(41,1) = 'So_t'
    lookup(41,2) = 'Sea surface temperature'
    lookup(41,3) = 'sea_surface_temperature'
    lookup(41,4) = 'K'

    lookup(42,1) = 'Sl_snowh'
    lookup(42,2) = 'Surface snow water equivalent'
    lookup(42,3) = 'surface_snow_water_equivalent'
    lookup(42,4) = 'm'

    lookup(43,1) = 'Sx_lfrac'
    lookup(43,2) = 'Surface land fraction'
    lookup(43,3) = 'land_area_fraction'
    lookup(43,4) = 'unitless'

    lookup(44,1) = 'Sx_ifrac'
    lookup(44,2) = 'Surface ice fraction'
    lookup(44,3) = 'sea_ice_area_fraction'
    lookup(44,4) = 'unitless'

    lookup(45,1) = 'Sx_ofrac'
    lookup(45,2) = 'Surface ocean fraction'
    lookup(45,3) = 'sea_area_fraction'
    lookup(45,4) = 'unitless'

    lookup(46,1) = 'So_ustar'
    lookup(46,2) = 'Surface fraction velocity in ocean'
    lookup(46,3) = 'fraction_velocity'
    lookup(46,4) = 'm s-1'

    lookup(47,1) = 'So_re'
    lookup(47,2) = 'Square of exch. coeff (tracers)'
!   lookup(47,3) = ''
!   lookup(47,4) = ''

    lookup(48,1) = 'So_ssq'
    lookup(48,2) = 'Surface saturation specific humidity in ocean'
    lookup(48,3) = 'specific_humidity_at_saturation'
    lookup(48,4) = 'kg kg-1'

    lookup(49,1) = 'Sl_fv'
    lookup(49,2) = 'Surface fraction velocity in land'
    lookup(49,3) = 'fraction_velocity'
    lookup(49,4) = 'm s-1'

    lookup(50,1) = 'Sl_ram1'
    lookup(50,2) = 'Aerodynamical resistance'
!   lookup(50,3) = ''
!   lookup(50,4) = ''

    lookup(51,1) = 'Faxx_taux'
    lookup(51,2) = 'Zonal surface stress'
    lookup(51,3) = 'surface_downward_eastward_stress'
    lookup(51,4) = 'N m-2'

    lookup(52,1) = 'Faxx_tauy'
    lookup(52,2) = 'Meridional surface stress'
    lookup(52,3) = 'surface_downward_northward_stress'
    lookup(52,4) = 'N m-2'

    lookup(53,1) = 'Faxx_lat'
    lookup(53,2) = 'Surface latent heat flux'
    lookup(53,3) = 'surface_upward_latent_heat_flux'
    lookup(53,4) = 'W m-2'

    lookup(54,1) = 'Faxx_sen'
    lookup(54,2) = 'Surface sensible heat flux'
    lookup(54,3) = 'surface_upward_sensible_heat_flux'
    lookup(54,4) = 'W m-2'

    lookup(55,1) = 'Faxx_lwup'
    lookup(55,2) = 'Surface upward longwave heat flux'
    lookup(55,3) = 'surface_net_upward_longwave_flux'
    lookup(55,4) = 'W m-2'

    lookup(56,1) = 'Faxx_evap'
    lookup(56,2) = 'Evaporation water flux'
    lookup(56,3) = 'water_evaporation_flux'
    lookup(56,4) = 'kg m-2 s-1'

    lookup(57,1) = 'Fall_flxdst1'
    lookup(57,2) = 'Dust flux (particle bin number 1)'
    lookup(57,3) = 'dust_flux'
    lookup(57,4) = 'kg m-2 s-1'

    lookup(58,1) = 'Fall_flxdst2'
    lookup(58,2) = 'Dust flux (particle bin number 2)'
    lookup(58,3) = 'dust_flux'
    lookup(58,4) = 'kg m-2 s-1'

    lookup(59,1) = 'Fall_flxdst3'
    lookup(59,2) = 'Dust flux (particle bin number 3)'
    lookup(59,3) = 'dust_flux'
    lookup(59,4) = 'kg m-2 s-1'

    lookup(60,1) = 'Fall_flxdst4'
    lookup(60,2) = 'Dust flux (particle bin number 4)'
    lookup(60,3) = 'dust_flux'
    lookup(60,4) = 'kg m-2 s-1'

    lookup(61,1) = 'Si_t'
    lookup(61,2) = 'Surface temperature of snow and ice surface'
    lookup(61,3) = 'surface_temperature_due_to_snow_and_ice'
    lookup(61,4) = 'K'

    lookup(62,1) = 'Si_tref'
    lookup(62,2) = 'Reference temperature at 2 meters'
    lookup(62,3) = 'air_temperature'
    lookup(62,4) = 'K'

    lookup(63,1) = 'Si_qref'
    lookup(63,2) = 'Reference specific humidity at 2 meters'
    lookup(63,3) = 'specific_humidity'
    lookup(63,4) = 'kg kg-1'

    lookup(64,1) = 'Si_ifrac'
    lookup(64,2) = 'Fractional ice coverage wrt ocean'
    lookup(64,3) = 'sea_ice_area_fraction'
    lookup(64,4) = 'unitless'

    lookup(65,1) = 'Si_avsdr'
    lookup(65,2) = 'Direct albedo (visible radiation)'
    lookup(65,3) = 'surface_direct_albedo_in_sea_ice_due_to_visible_radiation'
    lookup(65,4) = 'unitless'

    lookup(66,1) = 'Si_anidr'
    lookup(66,2) = 'Direct albedo (near-infrared radiation)'
    lookup(66,3) = 'surface_direct_albedo_in_sea_ice_due_to_near_infrared_radiation'
    lookup(66,4) = 'unitless'

    lookup(67,1) = 'Si_avsdf'
    lookup(67,2) = 'Diffuse albedo (visible radiation)'
    lookup(67,3) = 'surface_diffuse_albedo_in_sea_ice_due_to_visible_radiation'
    lookup(67,4) = 'unitless'

    lookup(68,1) = 'Si_anidf'
    lookup(68,2) = 'Diffuse albedo (near-infrared radiation)'
    lookup(68,3) = 'surface_diffuse_albedo_in_sea_ice_due_to_near_infrared_radiation'
    lookup(68,4) = 'unitless'

    lookup(69,1) = 'Si_u10'
    lookup(69,2) = 'Sea ice 10 meter wind'
    lookup(69,3) = 'sea_ice_10_meter_wind'
    lookup(69,4) = 'm'

    lookup(70,1) = 'Faii_taux'
    lookup(70,2) = 'Zonal surface stress'
    lookup(70,3) = 'surface_downward_eastward_stress'
    lookup(70,4) = 'N m-2'

    lookup(71,1) = 'Faii_tauy'
    lookup(71,2) = 'Meridional surface stress'
    lookup(71,3) = 'surface_downward_northward_stress'
    lookup(71,4) = 'N m-2'

    lookup(72,1) = 'Faii_lat'
    lookup(72,2) = 'Latent heat flux'
    lookup(72,3) = 'surface_upward_latent_heat_flux'
    lookup(72,4) = 'W m-2'

    lookup(73,1) = 'Faii_sen'
    lookup(73,2) = 'Sensible heat flux'
    lookup(73,3) = 'surface_upward_sensible_heat_flux'
    lookup(73,4) = 'W m-2'

    lookup(74,1) = 'Faii_lwup'
    lookup(74,2) = 'Outgoing longwave heat flux'
    lookup(74,3) = 'surface_net_upward_longwave_flux'
    lookup(74,4) = 'W m-2'

    lookup(75,1) = 'Faii_evap'
    lookup(75,2) = 'Evaporative water flux'
    lookup(75,3) = 'water_evaporation_flux'
    lookup(75,4) = 'kg m-2 s-1 '

    lookup(76,1) = 'Faii_swnet'
    lookup(76,2) = 'Net absorbed shortwave radiation in snow/ocean/ice'
    lookup(76,3) = 'surface_net_downward_shortwave_flux'
    lookup(76,4) = 'W m-2'

    lookup(77,1) = 'Fioi_swpen'
    lookup(77,2) = 'Net shortwave radiation penetrating into ice and ocean'
    lookup(77,3) = 'net_downward_shortwave_flux_in_sea_ice_due_to_penetration'
    lookup(77,4) = 'W m-2'

    lookup(78,1) = 'Fioi_melth'
    lookup(78,2) = 'Heat flux from melting ice'
    lookup(78,3) = 'surface_snow_and_ice_melt_heat_flux_in_sea_ice'
    lookup(78,4) = 'W m-2'

    lookup(79,1) = 'Fioi_meltw'
    lookup(79,2) = 'Water flux from melting ice'
    lookup(79,3) = 'water_flux_out_of_sea_ice_due_to_sea_ice_thermodynamics'
    lookup(79,4) = 'kg m-2 s-1'

    lookup(80,1) = 'Fioi_salt'
    lookup(80,2) = 'Salt flux from melting ice'
    lookup(80,3) = 'downward_sea_ice_basal_salt_flux'
    lookup(80,4) = 'kg m-2 s-1'

    lookup(81,1) = 'Fioi_taux'
    lookup(81,2) = 'Zonal surface stress'
    lookup(81,3) = 'surface_downward_eastward_stress'
    lookup(81,4) = 'N m-2'

    lookup(82,1) = 'Fioi_tauy'
    lookup(82,2) = 'Meridional surface stress'
    lookup(82,3) = 'surface_downward_northward_stress'
    lookup(82,4) = 'N m-2'

    lookup(83,1) = 'So_s'
    lookup(83,2) = 'Sea surface salinity'
    lookup(83,3) = 'sea_surface_salinity'
    lookup(83,4) = 'g kg-1'

    lookup(84,1) = 'So_u'
    lookup(84,2) = 'Zonal sea water velocity'
    lookup(84,3) = 'eastward_sea_water_velocity'
    lookup(84,4) = 'm s-1'

    lookup(85,1) = 'So_v'
    lookup(85,2) = 'Meridional sea water velocity'
    lookup(85,3) = 'northward_sea_water_velocity'
    lookup(85,4) = 'm s-1'

    lookup(86,1) = 'So_dhdx'
    lookup(86,2) = 'Zonal sea surface slope'
    lookup(86,3) = 'sea_surface_eastward_slope'
    lookup(86,4) = 'm m-1'

    lookup(87,1) = 'So_dhdy'
    lookup(87,2) = 'Meridional sea surface slope'
    lookup(87,3) = 'sea_surface_northward_slope'
    lookup(87,4) = 'm m-1'

    lookup(88,1) = 'Fioo_q'
    lookup(88,2) = 'Ocean freeze (q>0) or melt (q<0) potential'
    lookup(88,3) = 'surface_snow_and_ice_melt_heat_flux'
    lookup(88,4) = 'W m-2'
    
    lookup(89,1) = 'Faxa_rain'
    lookup(89,2) = 'Precipitation (liquid)'
    lookup(89,3) = 'rainfall_flux'
    lookup(89,4) = 'kg m-2 s-1'

    lookup(90,1) = 'Faxa_snow'
    lookup(90,2) = 'Precipitation (frozen)'
    lookup(90,3) = 'snowfall_flux'
    lookup(90,4) = 'kg m-2 s-1'

    lookup(91,1) = 'Sl_t'
    lookup(91,2) = 'Surface temperature'
    lookup(91,3) = 'surface_temperature'
    lookup(91,4) = 'K'

    lookup(92,1) = 'Sl_tref'
    lookup(92,2) = 'Reference temperature at 2 meters'
    lookup(92,3) = 'air_temperature'
    lookup(92,4) = 'K'

    lookup(93,1) = 'Sl_qref'
    lookup(93,2) = 'Reference specific humidity at 2 meters'
    lookup(93,3) = 'specific_humidity'
    lookup(93,4) = 'kg kg-1'

    lookup(94,1) = 'Sl_avsdr'
    lookup(94,2) = 'Direct albedo (visible radiation)'
    lookup(94,3) = 'surface_direct_albedo_in_land_due_to_visible_radiation'
    lookup(94,4) = 'unitless'

    lookup(95,1) = 'Sl_anidr'
    lookup(95,2) = 'Direct albedo (near-infrared radiation)'
    lookup(95,3) = 'surface_direct_albedo_in_land_due_to_near_infrared_radiation'
    lookup(95,4) = 'unitless'

    lookup(96,1) = 'Sl_avsdf'
    lookup(96,2) = 'Diffuse albedo (visible radiation)'
    lookup(96,3) = 'surface_diffuse_albedo_in_land_due_to_visible_radiation'
    lookup(96,4) = 'unitless'

    lookup(97,1) = 'Sl_anidf'
    lookup(97,2) = 'Diffuse albedo (near-infrared radiation)'
    lookup(97,3) = 'surface_diffuse_albedo_in_land_due_to_near_infrared_radiation'
    lookup(97,4) = 'unitless'

    lookup(98,1) = 'Sl_landfrac'
    lookup(98,2) = 'fractional land'
    lookup(98,3) = 'land_area_fraction'
    lookup(98,4) = 'unitless'

    lookup(99,1) = 'Fall_taux'
    lookup(99,2) = 'Zonal surface stress'
    lookup(99,3) = 'surface_downward_eastward_stress'
    lookup(99,4) = 'N m-2'

    lookup(100,1) = 'Fall_tauy'
    lookup(100,2) = 'Meridional surface stress'
    lookup(100,3) = 'surface_downward_northward_stress'
    lookup(100,4) = 'N m-2'

    lookup(101,1) = 'Fall_lat'
    lookup(101,2) = 'Latent heat flux'
    lookup(101,3) = 'surface_upward_latent_heat_flux'
    lookup(101,4) = 'W m-2'

    lookup(102,1) = 'Fall_sen'
    lookup(102,2) = 'Sensible heat flux'
    lookup(102,3) = 'surface_upward_sensible_heat_flux'
    lookup(102,4) = 'W m-2'

    lookup(103,1) = 'Fall_lwup'
    lookup(103,2) = 'Upward longwave heat flux'
    lookup(103,3) = 'surface_net_upward_longwave_flux'
    lookup(103,4) = 'W m-2'

    lookup(104,1) = 'Fall_evap'
    lookup(104,2) = 'Evaporation water flux'
    lookup(104,3) = 'water_evaporation_flux'
    lookup(104,4) = 'kg m-2 s-1'

    lookup(105,1) = 'Fall_swnet'
    lookup(105,2) = 'Net absorbed shortwave radiation'
    lookup(105,3) = 'surface_net_downward_shortwave_flux' 
    lookup(105,4) = 'W m-2'

    lookup(106,1) = 'Sx_duu10n'
    lookup(106,2) = 'Wind speed squared at 10 meters'
    lookup(106,3) = 'square_of_wind_speed'
    lookup(106,4) = 'm2 s-2'

    lookup(107,1) = 'Foxx_taux'
    lookup(107,2) = 'Zonal surface stress'
    lookup(107,3) = 'surface_downward_eastward_stress'
    lookup(107,4) = 'Pa'

    lookup(108,1) = 'Foxx_tauy'
    lookup(108,2) = 'Meridional surface stress'
    lookup(108,3) = 'surface_downward_northward_stress'
    lookup(108,4) = 'Pa'

    lookup(109,1) = 'Foxx_swnet'
    lookup(109,2) = 'Net shortwave radiation'
    lookup(109,3) = 'surface_net_downward_shortwave_flux'
    lookup(109,4) = 'W m-2'

    lookup(110,1) = 'Foxx_lat'
    lookup(110,2) = 'Downward latent heat flux'
    lookup(110,3) = 'surface_downward_latent_heat_flux'
    lookup(110,4) = 'W m-2'

    lookup(111,1) = 'Foxx_sen'
    lookup(111,2) = 'Downward sensible heat flux'
    lookup(111,3) = 'surface_downward_sensible_heat_flux'
    lookup(111,4) = 'W m-2'

    lookup(112,1) = 'Foxx_lwdn'
    lookup(112,2) = 'Downward longwave heat flux'
    lookup(112,3) = 'downwelling_longwave_flux'
    lookup(112,4) = 'W m-2'

    lookup(113,1) = 'Foxx_lwup'
    lookup(113,2) = 'Upward longwave heat flux'
    lookup(113,3) = 'upwelling_longwave_flux'
    lookup(113,4) = 'W m-2'

    lookup(114,1) = 'Foxx_melth'
    lookup(114,2) = 'Heat flux from melting'
    lookup(114,3) = 'surface_snow_melt_heat_flux'
    lookup(114,4) = 'W m-2'

    lookup(115,1) = 'Foxx_salt'
    lookup(115,2) = 'Salt flux'
    lookup(115,3) = 'virtual_salt_flux_into_sea_water'
    lookup(115,4) = 'kg m-2 s-1'

    lookup(116,1) = 'Foxx_prec'
    lookup(116,2) = 'Water flux (rain+snow)'
    lookup(116,3) = 'precipitation_flux'
    lookup(116,4) = 'kg m-2 s-1'

    lookup(117,1) = 'Foxx_snow'
    lookup(117,2) = 'Water flux due to snow'
    lookup(117,3) = 'surface_snow_melt_flux'
    lookup(117,4) = 'kg m-2 s-1'

    lookup(118,1) = 'Foxx_rain'
    lookup(118,2) = 'Water flux due to rain'
    lookup(118,3) = 'rainfall_flux'
    lookup(118,4) = 'kg m-2 s-1'

    lookup(119,1) = 'Foxx_evap'
    lookup(119,2) = 'Water flux due to evaporation'
    lookup(119,3) = 'water_evaporation_flux'
    lookup(119,4) = 'kg m-2 s-1'

    lookup(120,1) = 'Foxx_meltw'
    lookup(120,2) = 'Water flux due to melting'
    lookup(120,3) = 'surface_snow_melt_flux'
    lookup(120,4) = 'kg m-2 s-1'

    lookup(121,1) = 'Forr_roff'
    lookup(121,2) = 'Water flux due to runoff (liquid)'
    lookup(121,3) = 'water_flux_into_sea_water'
    lookup(121,4) = 'kg m-2 s-1'

    lookup(122,1) = 'Forr_ioff'
    lookup(122,2) = 'Water flux due to runoff (frozen)'
    lookup(122,3) = 'frozen_water_flux_into_sea_water'
    lookup(122,4) = 'kg m-2 s-1'

    lookup(123,1) = 'Foxx_bcphidry'
    lookup(123,2) = 'Hydrophylic black carbon dry deposition flux'
    lookup(123,3) = 'dry_deposition_flux_of_hydrophylic_black_carbon'
    lookup(123,4) = 'kg m-2 s-1'

    lookup(124,1) = 'Foxx_bcphodry'
    lookup(124,2) = 'Hydrophobic black carbon dry deposition flux'
    lookup(124,3) = 'dry_deposition_flux_of_hydrophobic_black_carbon'
    lookup(124,4) = 'kg m-2 s-1'

    lookup(125,1) = 'Foxx_bcphiwet'
    lookup(125,2) = 'Hydrophylic black carbon wet deposition flux'
    lookup(125,3) = 'wet_deposition_flux_of_hydrophylic_black_carbon'
    lookup(125,4) = 'kg m-2 s-1'

    lookup(126,1) = 'Foxx_ocphidry'
    lookup(126,2) = 'Hydrophylic organic carbon dry deposition flux'
    lookup(126,3) = 'dry_deposition_flux_of_hydrophylic_organic_carbon'
    lookup(126,4) = 'kg m-2 s-1'

    lookup(127,1) = 'Foxx_ocphodry'
    lookup(127,2) = 'Hydrophobic organic carbon dry deposition flux'
    lookup(127,3) = 'dry_deposition_flux_of_hydrophobic_organic_carbon'
    lookup(127,4) = 'kg m-2 s-1'

    lookup(128,1) = 'Foxx_ocphiwet'
    lookup(128,2) = 'Hydrophylic organic carbon wet deposition flux'
    lookup(128,3) = 'wet_deposition_flux_of_hydrophylic_organic_carbon'
    lookup(128,4) = 'kg m-2 s-1'

    lookup(129,1) = 'Foxx_dstwet1'
    lookup(129,2) = 'Dust wet deposition flux (size 1)'
    lookup(129,3) = 'wet_deposition_flux_of_dust'
    lookup(129,4) = 'kg m-2 s-1'

    lookup(130,1) = 'Foxx_dstwet2'
    lookup(130,2) = 'Dust wet deposition flux (size 2)'
    lookup(130,3) = 'wet_deposition_flux_of_dust'
    lookup(130,4) = 'kg m-2 s-1'

    lookup(131,1) = 'Foxx_dstwet3'
    lookup(131,2) = 'Dust wet deposition flux (size 3)'
    lookup(131,3) = 'wet_deposition_flux_of_dust'
    lookup(131,4) = 'kg m-2 s-1'

    lookup(132,1) = 'Foxx_dstwet4'
    lookup(132,2) = 'Dust wet deposition flux (size 4)'
    lookup(132,3) = 'wet_deposition_flux_of_dust'
    lookup(132,4) = 'kg m-2 s-1'

    lookup(133,1) = 'Foxx_dstdry1'
    lookup(133,2) = 'Dust dry deposition flux (size 1)'
    lookup(133,3) = 'dry_deposition_flux_of_dust'
    lookup(133,4) = 'kg m-2 s-1'

    lookup(134,1) = 'Foxx_dstdry2'
    lookup(134,2) = 'Dust dry deposition flux (size 2)'
    lookup(134,3) = 'dry_deposition_flux_of_dust'
    lookup(134,4) = 'kg m-2 s-1'

    lookup(135,1) = 'Foxx_dstdry3'
    lookup(135,2) = 'Dust dry deposition flux (size 3)'
    lookup(135,3) = 'dry_deposition_flux_of_dust'
    lookup(135,4) = 'kg m-2 s-1'

    lookup(136,1) = 'Foxx_dstdry4'
    lookup(136,2) = 'Dust dry deposition flux (size 4)'
    lookup(136,3) = 'dry_deposition_flux_of_dust'
    lookup(136,4) = 'kg m-2 s-1'

    lookup(137,1) = 'So_tref'
    lookup(137,2) = 'Reference temperature at 2 meters'
    lookup(137,3) = 'air_temperature'
    lookup(137,4) = 'K'

    lookup(138,1) = 'So_qref'
    lookup(138,2) = 'Reference specific humidity at 2 meters'
    lookup(138,3) = 'specific_humidity'
    lookup(138,4) = 'kg kg-1'

    lookup(139,1) = 'So_avsdr'
    lookup(139,2) = 'Direct albedo (visible radiation)'
    lookup(139,3) = 'surface_direct_albedo_in_sea_water_due_to_visible_radiation'
    lookup(139,4) = 'unitless'

    lookup(140,1) = 'So_anidr'
    lookup(140,2) = 'Direct albedo (near-infrared radiation)'
    lookup(140,3) = 'surface_direct_albedo_in_sea_water_due_to_near_infrared_radiation'
    lookup(140,4) = 'unitless'

    lookup(141,1) = 'So_avsdf'
    lookup(141,2) = 'Diffuse albedo (visible radiation)'
    lookup(141,3) = 'surface_diffuse_albedo_in_sea_water_due_to_visible_radiation'
    lookup(141,4) = 'unitless'

    lookup(142,1) = 'So_anidf'
    lookup(142,2) = 'Diffuse albedo (near-infrared radiation)'
    lookup(142,3) = 'surface_diffuse_albedo_in_sea_water_due_to_near_infrared_radiation'
    lookup(142,4) = 'unitless'

    lookup(143,1) = 'Faox_taux'
    lookup(143,2) = 'Zonal wind stress'
    lookup(143,3) = 'surface_downward_eastward_stress'
    lookup(143,4) = 'N m-2'

    lookup(144,1) = 'Faox_tauy'
    lookup(144,2) = 'Meridional wind stress'
    lookup(144,3) = 'surface_downward_northward_stress'
    lookup(144,4) = 'N m-2'

    lookup(145,1) = 'Faox_lat'
    lookup(145,2) = 'Latent heat flux'
    lookup(145,3) = 'surface_upward_latent_heat_flux'
    lookup(145,4) = 'W m-2'

    lookup(146,1) = 'Faox_sen'
    lookup(146,2) = 'Sensible heat flux'
    lookup(146,3) = 'surface_upward_sensible_heat_flux'
    lookup(146,4) = 'W m-2'

    lookup(147,1) = 'Faox_evap'
    lookup(147,2) = 'Evaporation water flux'
    lookup(147,3) = 'water_evaporation_flux'
    lookup(147,4) = 'kg m-2 s-1'

    lookup(148,1) = 'Faox_lwup'
    lookup(148,2) = 'Upward longwave heat flux'
    lookup(148,3) = 'surface_net_upward_longwave_flux'
    lookup(148,4) = 'W m-2'

    lookup(149,1) = 'Sg_frac01'
    lookup(149,2) = 'Fraction of glacier area (elevation class 1)'
    lookup(149,3) = 'glacier_area_fraction'
    lookup(149,4) = 'unitless'    

    lookup(150,1) = 'Sg_topo01'
    lookup(150,2) = 'Surface height of glacier (elevation class 1)'
    lookup(150,3) = 'height'
    lookup(150,4) = 'm'

    lookup(151,1) = 'Sg_frac02'
    lookup(151,2) = 'Fraction of glacier area (elevation class 2)'
    lookup(151,3) = 'glacier_area_fraction'
    lookup(151,4) = 'unitless'

    lookup(152,1) = 'Sg_topo02'
    lookup(152,2) = 'Surface height of glacier (elevation class 2)'
    lookup(152,3) = 'height'
    lookup(152,4) = 'm'

    lookup(153,1) = 'Sg_frac03'
    lookup(153,2) = 'Fraction of glacier area (elevation class 3)'
    lookup(153,3) = 'glacier_area_fraction'
    lookup(153,4) = 'unitless'

    lookup(154,1) = 'Sg_topo03'
    lookup(154,2) = 'Surface height of glacier (elevation class 3)'
    lookup(154,3) = 'height'
    lookup(154,4) = 'm'

    lookup(155,1) = 'Sg_frac04'
    lookup(155,2) = 'Fraction of glacier area (elevation class 4)'
    lookup(155,3) = 'glacier_area_fraction'
    lookup(155,4) = 'unitless'

    lookup(156,1) = 'Sg_topo04'
    lookup(156,2) = 'Surface height of glacier (elevation class 4)'
    lookup(156,3) = 'height'
    lookup(156,4) = 'm'

    lookup(157,1) = 'Sg_frac05'
    lookup(157,2) = 'Fraction of glacier area (elevation class 5)'
    lookup(157,3) = 'glacier_area_fraction'
    lookup(157,4) = 'unitless'

    lookup(158,1) = 'Sg_topo05'
    lookup(158,2) = 'Surface height of glacier (elevation class 5)'
    lookup(158,3) = 'height'
    lookup(158,4) = 'm'

    lookup(159,1) = 'Sg_frac06'
    lookup(159,2) = 'Fraction of glacier area (elevation class 6)'
    lookup(159,3) = 'glacier_area_fraction'
    lookup(159,4) = 'unitless'

    lookup(160,1) = 'Sg_topo06'
    lookup(160,2) = 'Surface height of glacier (elevation class 6)'
    lookup(160,3) = 'height'
    lookup(160,4) = 'm'

    lookup(161,1) = 'Sg_frac07'
    lookup(161,2) = 'Fraction of glacier area (elevation class 7)'
    lookup(161,3) = 'glacier_area_fraction'
    lookup(161,4) = 'unitless'

    lookup(162,1) = 'Sg_topo07'
    lookup(162,2) = 'Surface height of glacier (elevation class 7)'
    lookup(162,3) = 'height'
    lookup(162,4) = 'm'

    lookup(163,1) = 'Sg_frac08'
    lookup(163,2) = 'Fraction of glacier area (elevation class 8)'
    lookup(163,3) = 'glacier_area_fraction'
    lookup(163,4) = 'unitless'

    lookup(164,1) = 'Sg_topo08'
    lookup(164,2) = 'Surface height of glacier (elevation class 8)'
    lookup(164,3) = 'height'
    lookup(164,4) = 'm'

    lookup(165,1) = 'Sg_frac09'
    lookup(165,2) = 'Fraction of glacier area (elevation class 9)'
    lookup(165,3) = 'glacier_area_fraction'
    lookup(165,4) = 'unitless'

    lookup(166,1) = 'Sg_topo09'
    lookup(166,2) = 'Surface height of glacier (elevation class 9)'
    lookup(166,3) = 'height'
    lookup(166,4) = 'm'

    lookup(167,1) = 'Sg_frac10'
    lookup(167,2) = 'Fraction of glacier area (elevation class 10)'
    lookup(167,3) = 'glacier_area_fraction'
    lookup(167,4) = 'unitless'

    lookup(168,1) = 'Sg_topo10'
    lookup(168,2) = 'Surface height of glacier (elevation class 10)'
    lookup(168,3) = 'height'
    lookup(168,4) = 'm'

    lookup(169,1) = 'Fsgg_rofi01'
    lookup(169,2) = 'Ice runoff flux (elevation class 1)'
    lookup(169,3) = 'ice_runoff_flux_in_glacier'
    lookup(169,4) = 'kg m-2 s-1'

    lookup(170,1) = 'Fsgg_rofl01'
    lookup(170,2) = 'Liquid runoff flux (elevation class 1)'
    lookup(170,3) = 'liquid_runoff_flux_in_glacier'
    lookup(170,4) = 'kg m-2 s-1'

    lookup(171,1) = 'Fsgg_hflx01'
    lookup(171,2) = 'Downward heat flux from glacier interior (elevation class 1)'
    lookup(171,3) = 'downward_heat_flux_in_glacier'
    lookup(171,4) = 'W m-2'    

    lookup(172,1) = 'Fsgg_rofi02'
    lookup(172,2) = 'Ice runoff flux (elevation class 2)'
    lookup(172,3) = 'ice_runoff_flux_in_glacier'
    lookup(172,4) = 'kg m-2 s-1'

    lookup(173,1) = 'Fsgg_rofl02'
    lookup(173,2) = 'Liquid runoff flux (elevation class 2)'
    lookup(173,3) = 'liquid_runoff_flux_in_glacier'
    lookup(173,4) = 'kg m-2 s-1'

    lookup(174,1) = 'Fsgg_hflx02'
    lookup(174,2) = 'Downward heat flux from glacier interior (elevation class 2)'
    lookup(174,3) = 'downward_heat_flux_in_glacier'
    lookup(174,4) = 'W m-2'

    lookup(175,1) = 'Fsgg_rofi03'
    lookup(175,2) = 'Ice runoff flux (elevation class 3)'
    lookup(175,3) = 'ice_runoff_flux_in_glacier'
    lookup(175,4) = 'kg m-2 s-1'

    lookup(176,1) = 'Fsgg_rofl03'
    lookup(176,2) = 'Liquid runoff flux (elevation class 3)'
    lookup(176,3) = 'liquid_runoff_flux_in_glacier'
    lookup(176,4) = 'kg m-2 s-1'

    lookup(177,1) = 'Fsgg_hflx03'
    lookup(177,2) = 'Downward heat flux from glacier interior (elevation class 3)'
    lookup(177,3) = 'downward_heat_flux_in_glacier'
    lookup(177,4) = 'W m-2'

    lookup(178,1) = 'Fsgg_rofi04'
    lookup(178,2) = 'Ice runoff flux (elevation class 4)'
    lookup(178,3) = 'ice_runoff_flux_in_glacier'
    lookup(178,4) = 'kg m-2 s-1'

    lookup(179,1) = 'Fsgg_rofl04'
    lookup(179,2) = 'Liquid runoff flux (elevation class 4)'
    lookup(179,3) = 'liquid_runoff_flux_in_glacier'
    lookup(179,4) = 'kg m-2 s-1'

    lookup(180,1) = 'Fsgg_hflx04'
    lookup(180,2) = 'Downward heat flux from glacier interior (elevation class 4)'
    lookup(180,3) = 'downward_heat_flux_in_glacier'
    lookup(180,4) = 'W m-2'

    lookup(181,1) = 'Fsgg_rofi05'
    lookup(181,2) = 'Ice runoff flux (elevation class 5)'
    lookup(181,3) = 'ice_runoff_flux_in_glacier'
    lookup(181,4) = 'kg m-2 s-1'

    lookup(182,1) = 'Fsgg_rofl05'
    lookup(182,2) = 'Liquid runoff flux (elevation class 5)'
    lookup(182,3) = 'liquid_runoff_flux_in_glacier'
    lookup(182,4) = 'kg m-2 s-1'

    lookup(183,1) = 'Fsgg_hflx05'
    lookup(183,2) = 'Downward heat flux from glacier interior (elevation class 5)'
    lookup(183,3) = 'downward_heat_flux_in_glacier'
    lookup(183,4) = 'W m-2'

    lookup(184,1) = 'Fsgg_rofi06'
    lookup(184,2) = 'Ice runoff flux (elevation class 6)'
    lookup(184,3) = 'ice_runoff_flux_in_glacier'
    lookup(184,4) = 'kg m-2 s-1'

    lookup(185,1) = 'Fsgg_rofl06'
    lookup(185,2) = 'Liquid runoff flux (elevation class 6)'
    lookup(185,3) = 'liquid_runoff_flux_in_glacier'
    lookup(185,4) = 'kg m-2 s-1'

    lookup(186,1) = 'Fsgg_hflx06'
    lookup(186,2) = 'Downward heat flux from glacier interior (elevation class 6)'
    lookup(186,3) = 'downward_heat_flux_in_glacier'
    lookup(186,4) = 'W m-2'

    lookup(187,1) = 'Fsgg_rofi07'
    lookup(187,2) = 'Ice runoff flux (elevation class 7)'
    lookup(187,3) = 'ice_runoff_flux_in_glacier'
    lookup(187,4) = 'kg m-2 s-1'

    lookup(188,1) = 'Fsgg_rofl07'
    lookup(188,2) = 'Liquid runoff flux (elevation class 7)'
    lookup(188,3) = 'liquid_runoff_flux_in_glacier'
    lookup(188,4) = 'kg m-2 s-1'

    lookup(189,1) = 'Fsgg_hflx07'
    lookup(189,2) = 'Downward heat flux from glacier interior (elevation class 7)'
    lookup(189,3) = 'downward_heat_flux_in_glacier'
    lookup(189,4) = 'W m-2'

    lookup(190,1) = 'Fsgg_rofi08'
    lookup(190,2) = 'Ice runoff flux (elevation class 8)'
    lookup(190,3) = 'ice_runoff_flux_in_glacier'
    lookup(190,4) = 'kg m-2 s-1'

    lookup(191,1) = 'Fsgg_rofl08'
    lookup(191,2) = 'Liquid runoff flux (elevation class 8)'
    lookup(191,3) = 'liquid_runoff_flux_in_glacier'
    lookup(191,4) = 'kg m-2 s-1'

    lookup(192,1) = 'Fsgg_hflx08'
    lookup(192,2) = 'Downward heat flux from glacier interior (elevation class 8)'
    lookup(192,3) = 'downward_heat_flux_in_glacier'
    lookup(192,4) = 'W m-2'

    lookup(193,1) = 'Fsgg_rofi09'
    lookup(193,2) = 'Ice runoff flux (elevation class 9)'
    lookup(193,3) = 'ice_runoff_flux_in_glacier'
    lookup(193,4) = 'kg m-2 s-1'

    lookup(194,1) = 'Fsgg_rofl09'
    lookup(194,2) = 'Liquid runoff flux (elevation class 9)'
    lookup(194,3) = 'liquid_runoff_flux_in_glacier'
    lookup(194,4) = 'kg m-2 s-1'

    lookup(195,1) = 'Fsgg_hflx09'
    lookup(195,2) = 'Downward heat flux from glacier interior (elevation class 9)'
    lookup(195,3) = 'downward_heat_flux_in_glacier'
    lookup(195,4) = 'W m-2'

    lookup(196,1) = 'Fsgg_rofi10'
    lookup(196,2) = 'Ice runoff flux (elevation class 10)'
    lookup(196,3) = 'ice_runoff_flux_in_glacier'
    lookup(196,4) = 'kg m-2 s-1'

    lookup(197,1) = 'Fsgg_rofl10'
    lookup(197,2) = 'Liquid runoff flux (elevation class 10)'
    lookup(197,3) = 'liquid_runoff_flux_in_glacier'
    lookup(197,4) = 'kg m-2 s-1'

    lookup(198,1) = 'Fsgg_hflx10'
    lookup(198,2) = 'Downward heat flux from glacier interior (elevation class 10)'
    lookup(198,3) = 'downward_heat_flux_in_glacier'
    lookup(198,4) = 'W m-2'

    lookup(199,1) = 'Ss_tsrf01'
    lookup(199,2) = 'Surface temperature (elevation class 1)'
    lookup(199,3) = 'surface_temperature'
    lookup(199,4) = 'deg C'

    lookup(200,1) = 'Ss_topo01'
    lookup(200,2) = 'Surface height of glacier (elevation class 1)'
    lookup(200,3) = 'height'
    lookup(200,4) = 'm'

    lookup(201,1) = 'Ss_tsrf02'
    lookup(201,2) = 'Surface temperature (elevation class 2)'
    lookup(201,3) = 'surface_temperature'
    lookup(201,4) = 'deg C'

    lookup(202,1) = 'Ss_topo02'
    lookup(202,2) = 'Surface height of glacier (elevation class 2)'
    lookup(202,3) = 'height'
    lookup(202,4) = 'm'

    lookup(203,1) = 'Ss_tsrf03'
    lookup(203,2) = 'Surface temperature (elevation class 3)'
    lookup(203,3) = 'surface_temperature'
    lookup(203,4) = 'deg C'

    lookup(204,1) = 'Ss_topo03'
    lookup(204,2) = 'Surface height of glacier (elevation class 3)'
    lookup(204,3) = 'height'
    lookup(204,4) = 'm'

    lookup(205,1) = 'Ss_tsrf04'
    lookup(205,2) = 'Surface temperature (elevation class 4)'
    lookup(205,3) = 'surface_temperature'
    lookup(205,4) = 'deg C'

    lookup(206,1) = 'Ss_topo04'
    lookup(206,2) = 'Surface height of glacier (elevation class 4)'
    lookup(206,3) = 'height'
    lookup(206,4) = 'm'

    lookup(207,1) = 'Ss_tsrf05'
    lookup(207,2) = 'Surface temperature (elevation class 5)'
    lookup(207,3) = 'surface_temperature'
    lookup(207,4) = 'deg C'

    lookup(208,1) = 'Ss_topo05'
    lookup(208,2) = 'Surface height of glacier (elevation class 5)'
    lookup(208,3) = 'height'
    lookup(208,4) = 'm'

    lookup(209,1) = 'Ss_tsrf06'
    lookup(209,2) = 'Surface temperature (elevation class 6)'
    lookup(209,3) = 'surface_temperature'
    lookup(209,4) = 'deg C'

    lookup(210,1) = 'Ss_topo06'
    lookup(210,2) = 'Surface height of glacier (elevation class 6)'
    lookup(210,3) = 'height'
    lookup(210,4) = 'm'

    lookup(211,1) = 'Ss_tsrf07'
    lookup(211,2) = 'Surface temperature (elevation class 7)'
    lookup(211,3) = 'surface_temperature'
    lookup(211,4) = 'deg C'

    lookup(212,1) = 'Ss_topo07'
    lookup(212,2) = 'Surface height of glacier (elevation class 7)'
    lookup(212,3) = 'height'
    lookup(212,4) = 'm'

    lookup(213,1) = 'Ss_tsrf08'
    lookup(213,2) = 'Surface temperature (elevation class 8)'
    lookup(213,3) = 'surface_temperature'
    lookup(213,4) = 'deg C'

    lookup(214,1) = 'Ss_topo08'
    lookup(214,2) = 'Surface height of glacier (elevation class 8)'
    lookup(214,3) = 'height'
    lookup(214,4) = 'm'

    lookup(215,1) = 'Ss_tsrf09'
    lookup(215,2) = 'Surface temperature (elevation class 9)'
    lookup(215,3) = 'surface_temperature'
    lookup(215,4) = 'deg C'

    lookup(216,1) = 'Ss_topo09'
    lookup(216,2) = 'Surface height of glacier (elevation class 9)'
    lookup(216,3) = 'height'
    lookup(216,4) = 'm'

    lookup(217,1) = 'Ss_tsrf10'
    lookup(217,2) = 'Surface temperature (elevation class 10)'
    lookup(217,3) = 'surface_temperature'
    lookup(217,4) = 'deg C'

    lookup(218,1) = 'Ss_topo10'
    lookup(218,2) = 'Surface height of glacier (elevation class 10)'
    lookup(218,3) = 'height'
    lookup(218,4) = 'm'

    lookup(219,1) = 'Fgss_qice01'
    lookup(219,2) = 'New glacier ice flux (elevation class 1)'
    lookup(219,3) = 'ice_flux_out_of_glacier'
    lookup(219,4) = 'kg m-2 s-1'

    lookup(220,1) = 'Fgss_qice02'
    lookup(220,2) = 'New glacier ice flux (elevation class 2)'
    lookup(220,3) = 'ice_flux_out_of_glacier'
    lookup(220,4) = 'kg m-2 s-1'

    lookup(221,1) = 'Fgss_qice03'
    lookup(221,2) = 'New glacier ice flux (elevation class 3)'
    lookup(221,3) = 'ice_flux_out_of_glacier'
    lookup(221,4) = 'kg m-2 s-1'

    lookup(222,1) = 'Fgss_qice04'
    lookup(222,2) = 'New glacier ice flux (elevation class 4)'
    lookup(222,3) = 'ice_flux_out_of_glacier'
    lookup(222,4) = 'kg m-2 s-1'

    lookup(223,1) = 'Fgss_qice05'
    lookup(223,2) = 'New glacier ice flux (elevation class 5)'
    lookup(223,3) = 'ice_flux_out_of_glacier'
    lookup(223,4) = 'kg m-2 s-1'

    lookup(224,1) = 'Fgss_qice06'
    lookup(224,2) = 'New glacier ice flux (elevation class 6)'
    lookup(224,3) = 'ice_flux_out_of_glacier'
    lookup(224,4) = 'kg m-2 s-1'

    lookup(225,1) = 'Fgss_qice07'
    lookup(225,2) = 'New glacier ice flux (elevation class 7)'
    lookup(225,3) = 'ice_flux_out_of_glacier'
    lookup(225,4) = 'kg m-2 s-1'

    lookup(226,1) = 'Fgss_qice08'
    lookup(226,2) = 'New glacier ice flux (elevation class 8)'
    lookup(226,3) = 'ice_flux_out_of_glacier'
    lookup(226,4) = 'kg m-2 s-1'

    lookup(227,1) = 'Fgss_qice09'
    lookup(227,2) = 'New glacier ice flux (elevation class 9)'
    lookup(227,3) = 'ice_flux_out_of_glacier'
    lookup(227,4) = 'kg m-2 s-1'

    lookup(228,1) = 'Fgss_qice10'
    lookup(228,2) = 'New glacier ice flux (elevation class 10)'
    lookup(228,3) = 'ice_flux_out_of_glacier'
    lookup(228,4) = 'kg m-2 s-1'

    lookup(229,1) = 'Sa_co2prog'
    lookup(229,2) = 'Prognostic CO2 at the lowest model level'
!   lookup(229,3) = ''
    lookup(229,4) = '1e-6 mol/mol'

    lookup(230,1) = 'Sa_co2diag'
    lookup(230,2) = 'Diagnostic CO2 at the lowest model level'
!   lookup(230,3) = ''
    lookup(230,4) = '1e-6 mol/mol'

    lookup(231,1) = 'Faxx_fco2_lnd'
    lookup(231,2) = 'Surface flux of CO2 from land'
    lookup(231,3) = 'surface_upward_flux_of_carbon_dioxide_where_land'
    lookup(231,4) = 'moles m-2 s-1'

    lookup(232,1) = 'Fall_co2_lnd'
    lookup(232,2) = 'Surface flux of CO2 from land'
    lookup(232,3) = 'surface_upward_flux_of_carbon_dioxide_where_land'
    lookup(232,4) = 'moles m-2 s-1'

    lookup(233,1) = 'Faxx_fco2_ocn'
    lookup(233,2) = 'Surface flux of CO2 from ocean'
    lookup(233,3) = 'surface_upward_flux_of_carbon_dioxide_where_open_sea'
    lookup(233,4) = 'moles m-2 s-1'

    lookup(234,1) = 'Faoo_fco2'
    lookup(234,2) = 'Surface flux of CO2 from ocean'
    lookup(234,3) = 'surface_upward_flux_of_carbon_dioxide_where_open_sea'
    lookup(234,4) = 'moles m-2 s-1'

    lookup(235,1) = 'Faxx_fdms'
    lookup(235,2) = 'Surface flux of DMS'
    lookup(235,3) = 'surface_upward_flux_of_dimethyl_sulfide'
    lookup(235,4) = 'moles m-2 s-1'

    lookup(236,1) = 'Faoo_fdms'
    lookup(236,2) = 'Surface flux of DMS from ocean'
    lookup(236,3) = 'surface_upward_flux_of_dimethyl_sulfide_where_open_sea'
    lookup(236,4) = 'moles m-2 s-1'

    lookup(237,1) = 'lat'
!   lookup(237,2) = ''
    lookup(237,3) = 'latitude'
    lookup(237,4) = 'degrees north'

    lookup(238,1) = 'lon'
!   lookup(238,2) = ''
    lookup(238,3) = 'longitude'
    lookup(238,4) = 'degrees east'

    lookup(239,1) = 'area'
!   lookup(239,2) = ''
    lookup(239,3) = 'cell area'
    lookup(239,4) = 'm^2'

    lookup(240,1) = 'aream'
!   lookup(240,2) = ''
    lookup(240,3) = 'cell area from mapping file'
    lookup(240,4) = 'm^2'

    lookup(241,1) = 'mask'
!   lookup(241,2) = ''
    lookup(241,3) = 'mask'
    lookup(241,4) = 'unitless'

    lookup(242,1) = 'frac'
    lookup(242,2) = 'area_fraction'
    lookup(242,3) = 'area fraction'
    lookup(242,4) = 'unitless'

    lookup(243,1) = 'ascale'
!   lookup(243,2) = ''
    lookup(243,3) = 'area scale factor'
    lookup(243,4) = 'unitless'

    lookup(244,1) = 'Si_snowh'
    lookup(244,2) = 'Surface snow depth'
    lookup(244,3) = 'surface_snow_thickness'
    lookup(244,4) = 'm'

    lookup(245,1) = 'Sl_u10'
    lookup(245,2) = 'Land 10m wind'
    lookup(245,3) = 'Land_10m_wind'
    lookup(245,4) = 'm'

    lookup(246,1) = 'So_u10'
    lookup(246,2) = 'Ocean 10m wind'
    lookup(246,3) = 'Ocean_10m_wind'
    lookup(246,4) = 'm'

    lookup(246,1) = 'Sx_u10'
    lookup(246,2) = 'Merged 10m wind'
    lookup(246,3) = 'Merged_10m_wind'
    lookup(246,4) = 'm'
  end if  

  llongname = trim(unknown)
  lstdname  = trim(unknown)
  lunits    = trim(unknown)

  found = .false.

  if (.not.found) then
     i = 1
     do while (i <= nmax .and. .not.found)
        lshortname = trim(shortname)
        if (trim(lshortname) == trim(lookup(i,1))) then
           llongname = trim(lookup(i,2)) 
           lstdname = trim(lookup(i,3))
           lunits = trim(lookup(i,4))    
           found = .true.
        end if
        i = i + 1
     end do
  endif

  if (.not.found) then
     i = 1
     do while (i <= nmax .and. .not.found)
        n = shr_string_lastIndex(shortname,"_")
        lshortname = ""
        if (n < len_trim(shortname)) lshortname = shortname(n+1:len_trim(shortname))
        if (trim(lshortname) == trim(lookup(i,1))) then
           llongname = trim(lookup(i,2)) 
           lstdname = trim(lookup(i,3))
           lunits = trim(lookup(i,4))    
           found = .true.
        end if
        i = i + 1
     end do
  endif

  if (present(longname)) then
     longname = trim(llongname)
  endif
  if (present(stdname))  then
     stdname = trim(lstdname)
  endif
  if (present(units)) then
     units = trim(lunits)
  endif

end subroutine seq_flds_esmf_metadata_get

!===============================================================================

end module seq_flds_mod

