!===============================================================================
!BOP ===========================================================================
!
! !MODULE: seq_flds_mod -- coupler/comp list of exchange fields indices
!
! !DESCRIPTION:
! List of all supported coupler/component list of all possible exchanged 
! The actual experiment will only use those flds in seq_flds_mod
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
! !INTERFACE: ------------------------------------------------------------------

module seq_flds_indices

! !USES:

  use shr_string_mod
  use shr_sys_mod, only: shr_sys_abort
  use seq_flds_mod  
  use seq_drydep_mod, only: drydep_fields_token, lnd_drydep
  use mct_mod         

  implicit none
  save
  public

! !PUBLIC TYPES:

  ! domain

  integer :: nflds_dom

  ! atm -> drv

  integer :: index_a2x_Sa_z            ! bottom atm level height
  integer :: index_a2x_Sa_u            ! bottom atm level zon wind
  integer :: index_a2x_Sa_v            ! bottom atm level mer wind
  integer :: index_a2x_Sa_tbot         ! bottom atm level temp
  integer :: index_a2x_Sa_ptem         ! bottom atm level pot temp
  integer :: index_a2x_Sa_shum         ! bottom atm level spec hum
  integer :: index_a2x_Sa_dens         ! bottom atm level air den
  integer :: index_a2x_Sa_pbot         ! bottom atm level pressure
  integer :: index_a2x_Sa_pslv         ! sea level atm pressure
  integer :: index_a2x_Faxa_lwdn       ! downward lw heat flux
  integer :: index_a2x_Faxa_rainc      ! prec: liquid "convective"
  integer :: index_a2x_Faxa_rainl      ! prec: liquid "large scale"
  integer :: index_a2x_Faxa_snowc      ! prec: frozen "convective"
  integer :: index_a2x_Faxa_snowl      ! prec: frozen "large scale"
  integer :: index_a2x_Faxa_swndr      ! sw: nir direct  downward
  integer :: index_a2x_Faxa_swvdr      ! sw: vis direct  downward
  integer :: index_a2x_Faxa_swndf      ! sw: nir diffuse downward
  integer :: index_a2x_Faxa_swvdf      ! sw: vis diffuse downward
  integer :: index_a2x_Faxa_swnet      ! sw: net
  integer :: index_a2x_Faxa_bcphidry   ! flux: Black Carbon hydrophilic dry deposition
  integer :: index_a2x_Faxa_bcphodry   ! flux: Black Carbon hydrophobic dry deposition
  integer :: index_a2x_Faxa_bcphiwet   ! flux: Black Carbon hydrophilic wet deposition
  integer :: index_a2x_Faxa_ocphidry   ! flux: Organic Carbon hydrophilic dry deposition
  integer :: index_a2x_Faxa_ocphodry   ! flux: Organic Carbon hydrophobic dry deposition
  integer :: index_a2x_Faxa_ocphiwet   ! flux: Organic Carbon hydrophilic dry deposition
  integer :: index_a2x_Faxa_dstwet1    ! flux: Size 1 dust -- wet deposition
  integer :: index_a2x_Faxa_dstwet2    ! flux: Size 2 dust -- wet deposition
  integer :: index_a2x_Faxa_dstwet3    ! flux: Size 3 dust -- wet deposition
  integer :: index_a2x_Faxa_dstwet4    ! flux: Size 4 dust -- wet deposition
  integer :: index_a2x_Faxa_dstdry1    ! flux: Size 1 dust -- dry deposition
  integer :: index_a2x_Faxa_dstdry2    ! flux: Size 2 dust -- dry deposition
  integer :: index_a2x_Faxa_dstdry3    ! flux: Size 3 dust -- dry deposition
  integer :: index_a2x_Faxa_dstdry4    ! flux: Size 4 dust -- dry deposition
  integer :: index_a2x_Sa_co2prog      ! bottom atm level prognostic co2
  integer :: index_a2x_Sa_co2diag      ! bottom atm level diagnostic co2
  integer :: nflds_a2x

  ! drv -> atm

  integer :: index_x2a_Sx_t            ! surface temperature             
  integer :: index_x2a_So_t            ! sea surface temperature         
  integer :: index_x2a_Sx_lfrac        ! surface land fraction           
  integer :: index_x2a_Sx_ifrac        ! surface ice fraction            
  integer :: index_x2a_Sx_ofrac        ! surface ocn fraction            
  integer :: index_x2a_Sx_tref         ! 2m reference temperature        
  integer :: index_x2a_Sx_qref         ! 2m reference specific humidity  
  integer :: index_x2a_Sx_avsdr        ! albedo, visible, direct         
  integer :: index_x2a_Sx_anidr        ! albedo, near-ir, direct         
  integer :: index_x2a_Sx_avsdf        ! albedo, visible, diffuse        
  integer :: index_x2a_Sx_anidf        ! albedo, near-ir, diffuse        
  integer :: index_x2a_Sl_snowh        ! surface snow depth over land
  integer :: index_x2a_Si_snowh        ! surface snow depth over ice
  integer :: index_x2a_Sl_fv           ! friction velocity
  integer :: index_x2a_Sl_ram1         ! aerodynamical resistance
  integer :: index_x2a_Faxx_taux       ! wind stress, zonal              
  integer :: index_x2a_Faxx_tauy       ! wind stress, meridional         
  integer :: index_x2a_Faxx_lat        ! latent          heat flux       
  integer :: index_x2a_Faxx_sen        ! sensible        heat flux       
  integer :: index_x2a_Faxx_lwup       ! upward longwave heat flux       
  integer :: index_x2a_Faxx_evap       ! evaporation    water flux       
  integer :: index_x2a_Fall_flxdst1    ! dust flux size bin 1    
  integer :: index_x2a_Fall_flxdst2    ! dust flux size bin 2    
  integer :: index_x2a_Fall_flxdst3    ! dust flux size bin 3    
  integer :: index_x2a_Fall_flxdst4    ! dust flux size bin 4
  integer :: index_x2a_Faxx_flxvoc1    ! voc flux size bin 1    
  integer :: index_x2a_Faxx_flxvoc2    ! voc flux size bin 2    
  integer :: index_x2a_Faxx_flxvoc3    ! voc flux size bin 3    
  integer :: index_x2a_Faxx_flxvoc4    ! voc flux size bin 4
  integer :: index_x2a_Faxx_flxvoc5    ! voc flux size bin 5
  integer :: index_x2a_So_ustar	
  integer :: index_x2a_So_re
  integer :: index_x2a_So_ssq
  integer :: index_x2a_Sx_ddvel        ! dry deposition velocities
  integer :: index_x2a_Sx_u10          ! 10m wind
  !
  integer :: index_x2a_Si_avsdr        ! albedo, visible, direct         
  integer :: index_x2a_Si_anidr        ! albedo, near-ir, direct         
  integer :: index_x2a_Si_avsdf        ! albedo, visible, diffuse        
  integer :: index_x2a_Si_anidf        ! albedo, near-ir, diffuse        
  integer :: index_x2a_So_avsdr        ! albedo, visible, direct   ! (flux module)
  integer :: index_x2a_So_anidr        ! albedo, near-ir, direct   ! (flux module)
  integer :: index_x2a_So_avsdf        ! albedo, visible, diffuse  ! (flux module)
  integer :: index_x2a_So_anidf        ! albedo, near-ir, diffuse  ! (flux module)
  integer :: index_x2a_Faii_lat        ! latent          heat flux       
  integer :: index_x2a_Faii_sen        ! sensible        heat flux       
  integer :: index_x2a_Faii_lwup       ! upward longwave heat flux       
  integer :: index_x2a_Faox_lat        ! latent          heat flux ! (flux module)
  integer :: index_x2a_Faox_sen        ! sensible        heat flux ! (flux module)      
  integer :: index_x2a_Faox_lwup       ! upward longwave heat flux       
  integer :: index_x2a_Faxx_fco2_lnd   ! co2 flux from land   
  integer :: index_x2a_Faxx_fco2_ocn   ! co2 flux from ocean  
  integer :: index_x2a_Faxx_fdms       ! dms flux
  !
  integer :: nflds_x2a

  ! ice -> drv

  integer :: index_i2x_Si_ifrac        ! fractional ice coverage wrt ocean
  integer :: index_i2x_Si_sicthk       ! sea ice thickness (needed only for cam/som)
  integer :: index_i2x_Si_snowh        ! snow height (m)
  integer :: index_i2x_Si_t            ! temperature                     
  integer :: index_i2x_Si_tref         ! 2m reference temperature        
  integer :: index_i2x_Si_qref         ! 2m reference specific humidity  
  integer :: index_i2x_Si_avsdr        ! albedo: visible, direct         
  integer :: index_i2x_Si_avsdf        ! albedo: near ir, direct         
  integer :: index_i2x_Si_anidr        ! albedo: visible, diffuse        
  integer :: index_i2x_Si_anidf        ! albedo: near ir, diffuse        
  integer :: index_i2x_Si_u10          ! 10m wind
  integer :: index_i2x_Faii_lwup       ! upward longwave heat flux  
  integer :: index_i2x_Faii_lat        ! latent          heat flux  
  integer :: index_i2x_Faii_sen        ! sensible        heat flux      
  integer :: index_i2x_Faii_evap       ! evaporation    water flux      
  integer :: index_i2x_Faii_taux       ! wind stress, zonal            
  integer :: index_i2x_Faii_tauy       ! wind stress, meridional       
  integer :: index_i2x_Faii_swnet      ! sw: net
  integer :: index_i2x_Fioi_swpen      ! sw: net penetrating ice
  integer :: index_i2x_Fioi_melth      ! heat  flux from melting ice (<0)
  integer :: index_i2x_Fioi_meltw      ! water flux from melting ice
  integer :: index_i2x_Fioi_salt       ! salt  flux from meting  ice
  integer :: index_i2x_Fioi_taux       ! ice/ocn stress, zonal
  integer :: index_i2x_Fioi_tauy       ! ice/ocn stress, zonal
  integer :: nflds_i2x

  ! drv -> ice

  integer :: index_x2i_So_t            ! ocn layer temperature
  integer :: index_x2i_So_s            ! ocn salinity
  integer :: index_x2i_So_u            ! ocn u velocity
  integer :: index_x2i_So_v            ! ocn v velocity
  integer :: index_x2i_Sa_z            ! bottom atm level height
  integer :: index_x2i_Sa_u            ! bottom atm level zon wind
  integer :: index_x2i_Sa_v            ! bottom atm level mer wind
  integer :: index_x2i_Sa_tbot         ! bottom atm level temp
  integer :: index_x2i_Sa_pbot         ! bottom atm level pressure
  integer :: index_x2i_Sa_ptem         ! bottom atm level pot temp
  integer :: index_x2i_Sa_shum         ! bottom atm level spec hum
  integer :: index_x2i_Sa_dens         ! bottom atm level air den
  integer :: index_x2i_So_dhdx         ! ocn surface slope, zonal
  integer :: index_x2i_So_dhdy         ! ocn surface slope, meridional
  integer :: index_x2i_Faxa_lwdn       ! downward lw heat flux
  integer :: index_x2i_Faxa_rain       ! prec: liquid 
  integer :: index_x2i_Faxa_snow       ! prec: frozen 
  integer :: index_x2i_Faxa_swndr      ! sw: nir direct  downward
  integer :: index_x2i_Faxa_swvdr      ! sw: vis direct  downward
  integer :: index_x2i_Faxa_swndf      ! sw: nir diffuse downward
  integer :: index_x2i_Faxa_swvdf      ! sw: vis diffuse downward
  integer :: index_x2i_Faxa_swnet      ! sw: net
  integer :: index_x2i_Fioo_q          ! ocn freeze or melt heat  
  integer :: index_x2i_Faxa_bcphidry   ! flux: Black Carbon hydrophilic dry deposition
  integer :: index_x2i_Faxa_bcphodry   ! flux: Black Carbon hydrophobic dry deposition
  integer :: index_x2i_Faxa_bcphiwet   ! flux: Black Carbon hydrophilic wet deposition
  integer :: index_x2i_Faxa_ocphidry   ! flux: Organic Carbon hydrophilic dry deposition
  integer :: index_x2i_Faxa_ocphodry   ! flux: Organic Carbon hydrophobic dry deposition
  integer :: index_x2i_Faxa_ocphiwet   ! flux: Organic Carbon hydrophilic dry deposition
  integer :: index_x2i_Faxa_dstwet1    ! flux: Size 1 dust -- wet deposition
  integer :: index_x2i_Faxa_dstwet2    ! flux: Size 2 dust -- wet deposition
  integer :: index_x2i_Faxa_dstwet3    ! flux: Size 3 dust -- wet deposition
  integer :: index_x2i_Faxa_dstwet4    ! flux: Size 4 dust -- wet deposition
  integer :: index_x2i_Faxa_dstdry1    ! flux: Size 1 dust -- dry deposition
  integer :: index_x2i_Faxa_dstdry2    ! flux: Size 2 dust -- dry deposition
  integer :: index_x2i_Faxa_dstdry3    ! flux: Size 3 dust -- dry deposition
  integer :: index_x2i_Faxa_dstdry4    ! flux: Size 4 dust -- dry deposition
  integer :: index_x2i_Sa_co2prog      ! bottom atm level prognostic co2
  integer :: nflds_x2i

  ! hub atm/ocn flufxes and states (computed by flux module)	

  integer :: index_xao_So_tref    
  integer :: index_xao_So_qref    
  integer :: index_xao_So_avsdr   
  integer :: index_xao_So_avsdf   
  integer :: index_xao_So_anidr   
  integer :: index_xao_So_anidf   
  integer :: index_xao_Sx_duu10n 
  integer :: index_xao_So_u10
  integer :: index_xao_Faox_taux  
  integer :: index_xao_Faox_tauy   
  integer :: index_xao_Faox_lat   
  integer :: index_xao_Faox_sen   
  integer :: index_xao_Faox_evap  
  integer :: index_xao_Faox_lwup
  integer :: index_xao_So_ustar         ! optional
  integer :: index_xao_So_re            ! optional
  integer :: index_xao_So_ssq           ! optional
  integer :: nflds_xao

  ! ocn -> drv

  integer :: index_o2x_So_t      
  integer :: index_o2x_So_u
  integer :: index_o2x_So_v
  integer :: index_o2x_So_s
  integer :: index_o2x_So_dhdx
  integer :: index_o2x_So_dhdy
  integer :: index_o2x_Fioo_q
  integer :: index_o2x_Faoo_fco2
  integer :: index_o2x_Faoo_fdms
  integer :: nflds_o2x	

  ! drv -> ocn

  integer :: index_x2o_Si_ifrac        ! fractional ice wrt ocean
  integer :: index_x2o_Sx_duu10n
  integer :: index_x2o_Sa_pslv
  integer :: index_x2o_Sa_co2prog
  integer :: index_x2o_Sa_co2diag
  integer :: index_x2o_Foxx_taux  
  integer :: index_x2o_Foxx_tauy  
  integer :: index_x2o_Foxx_swnet 
  integer :: index_x2o_Foxx_sen   
  integer :: index_x2o_Foxx_lat   
  integer :: index_x2o_Foxx_lwdn  
  integer :: index_x2o_Foxx_lwup
  integer :: index_x2o_Foxx_melth 
  integer :: index_x2o_Foxx_salt
  integer :: index_x2o_Foxx_prec
  integer :: index_x2o_Foxx_snow
  integer :: index_x2o_Foxx_rain
  integer :: index_x2o_Foxx_evap
  integer :: index_x2o_Foxx_meltw 
  integer :: index_x2o_Forr_roff
  integer :: index_x2o_Forr_ioff
  integer :: index_x2o_Foxx_bcphidry   ! flux: Black Carbon hydrophilic dry deposition
  integer :: index_x2o_Foxx_bcphodry   ! flux: Black Carbon hydrophobic dry deposition
  integer :: index_x2o_Foxx_bcphiwet   ! flux: Black Carbon hydrophilic wet deposition
  integer :: index_x2o_Foxx_ocphidry   ! flux: Organic Carbon hydrophilic dry deposition
  integer :: index_x2o_Foxx_ocphodry   ! flux: Organic Carbon hydrophobic dry deposition
  integer :: index_x2o_Foxx_ocphiwet   ! flux: Organic Carbon hydrophilic dry deposition
  integer :: index_x2o_Foxx_dstwet1    ! flux: Size 1 dust -- wet deposition
  integer :: index_x2o_Foxx_dstwet2    ! flux: Size 2 dust -- wet deposition
  integer :: index_x2o_Foxx_dstwet3    ! flux: Size 3 dust -- wet deposition
  integer :: index_x2o_Foxx_dstwet4    ! flux: Size 4 dust -- wet deposition
  integer :: index_x2o_Foxx_dstdry1    ! flux: Size 1 dust -- dry deposition
  integer :: index_x2o_Foxx_dstdry2    ! flux: Size 2 dust -- dry deposition
  integer :: index_x2o_Foxx_dstdry3    ! flux: Size 3 dust -- dry deposition
  integer :: index_x2o_Foxx_dstdry4    ! flux: Size 4 dust -- dry deposition
  integer :: nflds_x2o

  ! lnd -> drv

  integer :: index_l2x_Sl_landfrac     ! land fraction
  integer :: index_l2x_Sl_t            ! temperature
  integer :: index_l2x_Sl_tref         ! 2m reference temperature
  integer :: index_l2x_Sl_qref         ! 2m reference specific humidity
  integer :: index_l2x_Sl_avsdr        ! albedo: direct , visible
  integer :: index_l2x_Sl_anidr        ! albedo: direct , near-ir
  integer :: index_l2x_Sl_avsdf        ! albedo: diffuse, visible
  integer :: index_l2x_Sl_anidf        ! albedo: diffuse, near-ir
  integer :: index_l2x_Sl_snowh        ! snow height
  integer :: index_l2x_Sl_u10          ! 10m wind
  integer :: index_l2x_Fall_taux       ! wind stress, zonal
  integer :: index_l2x_Fall_tauy       ! wind stress, meridional
  integer :: index_l2x_Fall_lat        ! latent          heat flux
  integer :: index_l2x_Fall_sen        ! sensible        heat flux
  integer :: index_l2x_Fall_lwup       ! upward longwave heat flux
  integer :: index_l2x_Fall_evap       ! evaporation     water flux
  integer :: index_l2x_Fall_swnet      ! heat flux       shortwave net       
  integer :: index_l2x_Fall_fco2_lnd   ! co2 flux **For testing set to 0
  integer :: index_l2x_Sl_fv           ! friction velocity  
  integer :: index_l2x_Sl_ram1         ! aerodynamical resistance
  integer :: index_l2x_Fall_flxdst1    ! dust flux size bin 1    
  integer :: index_l2x_Fall_flxdst2    ! dust flux size bin 2    
  integer :: index_l2x_Fall_flxdst3    ! dust flux size bin 3    
  integer :: index_l2x_Fall_flxdst4    ! dust flux size bin 4
  integer :: index_l2x_Fall_flxvoc1    ! voc flux size bin 1    
  integer :: index_l2x_Fall_flxvoc2    ! voc flux size bin 2    
  integer :: index_l2x_Fall_flxvoc3    ! voc flux size bin 3    
  integer :: index_l2x_Fall_flxvoc4    ! voc flux size bin 4
  integer :: index_l2x_Fall_flxvoc5    ! voc flux size bin 5
  integer :: index_l2x_Sl_ddvel        ! dry deposition velocities
  integer :: nflds_l2x

  ! roff to driver (part of land for now)

  integer :: index_r2x_Forr_roff       ! runoff to ocean
  integer :: index_r2x_Forr_ioff       ! runoff to ocean
  integer :: nflds_r2x

  ! drv -> lnd

  integer :: index_x2l_Sa_z            ! bottom atm level height
  integer :: index_x2l_Sa_u            ! bottom atm level zon wind
  integer :: index_x2l_Sa_v            ! bottom atm level mer wind
  integer :: index_x2l_Sa_ptem         ! bottom atm level pot temp
  integer :: index_x2l_Sa_shum         ! bottom atm level spec hum
  integer :: index_x2l_Sa_pbot         ! bottom atm level pressure
  integer :: index_x2l_Sa_tbot         ! bottom atm level temp
  integer :: index_x2l_Faxa_lwdn       ! downward lw heat flux
  integer :: index_x2l_Faxa_rainc      ! prec: liquid "convective"
  integer :: index_x2l_Faxa_rainl      ! prec: liquid "large scale"
  integer :: index_x2l_Faxa_snowc      ! prec: frozen "convective"
  integer :: index_x2l_Faxa_snowl      ! prec: frozen "large scale"
  integer :: index_x2l_Faxa_swndr      ! sw: nir direct  downward
  integer :: index_x2l_Faxa_swvdr      ! sw: vis direct  downward
  integer :: index_x2l_Faxa_swndf      ! sw: nir diffuse downward
  integer :: index_x2l_Faxa_swvdf      ! sw: vis diffuse downward
  integer :: index_x2l_Sa_co2prog      ! bottom atm level prognostic co2
  integer :: index_x2l_Sa_co2diag      ! bottom atm level diagnostic co2
  integer :: index_x2l_Faxa_bcphidry   ! flux: Black Carbon hydrophilic dry deposition
  integer :: index_x2l_Faxa_bcphodry   ! flux: Black Carbon hydrophobic dry deposition
  integer :: index_x2l_Faxa_bcphiwet   ! flux: Black Carbon hydrophilic wet deposition
  integer :: index_x2l_Faxa_ocphidry   ! flux: Organic Carbon hydrophilic dry deposition
  integer :: index_x2l_Faxa_ocphodry   ! flux: Organic Carbon hydrophobic dry deposition
  integer :: index_x2l_Faxa_ocphiwet   ! flux: Organic Carbon hydrophilic dry deposition
  integer :: index_x2l_Faxa_dstwet1    ! flux: Size 1 dust -- wet deposition
  integer :: index_x2l_Faxa_dstwet2    ! flux: Size 2 dust -- wet deposition
  integer :: index_x2l_Faxa_dstwet3    ! flux: Size 3 dust -- wet deposition
  integer :: index_x2l_Faxa_dstwet4    ! flux: Size 4 dust -- wet deposition
  integer :: index_x2l_Faxa_dstdry1    ! flux: Size 1 dust -- dry deposition
  integer :: index_x2l_Faxa_dstdry2    ! flux: Size 2 dust -- dry deposition
  integer :: index_x2l_Faxa_dstdry3    ! flux: Size 3 dust -- dry deposition
  integer :: index_x2l_Faxa_dstdry4    ! flux: Size 4 dust -- dry deposition
  integer :: nflds_x2l

  ! drv -> glc

  integer :: index_x2g_Ss_tsrf01
  integer :: index_x2g_Ss_topo01
  integer :: index_x2g_Fgss_qice01
  integer :: index_x2g_Ss_tsrf02
  integer :: index_x2g_Ss_topo02
  integer :: index_x2g_Fgss_qice02
  integer :: index_x2g_Ss_tsrf03
  integer :: index_x2g_Ss_topo03
  integer :: index_x2g_Fgss_qice03
  integer :: index_x2g_Ss_tsrf04
  integer :: index_x2g_Ss_topo04
  integer :: index_x2g_Fgss_qice04
  integer :: index_x2g_Ss_tsrf05
  integer :: index_x2g_Ss_topo05
  integer :: index_x2g_Fgss_qice05
  integer :: index_x2g_Ss_tsrf06
  integer :: index_x2g_Ss_topo06
  integer :: index_x2g_Fgss_qice06
  integer :: index_x2g_Ss_tsrf07
  integer :: index_x2g_Ss_topo07
  integer :: index_x2g_Fgss_qice07
  integer :: index_x2g_Ss_tsrf08
  integer :: index_x2g_Ss_topo08
  integer :: index_x2g_Fgss_qice08
  integer :: index_x2g_Ss_tsrf09
  integer :: index_x2g_Ss_topo09
  integer :: index_x2g_Fgss_qice09
  integer :: index_x2g_Ss_tsrf10
  integer :: index_x2g_Ss_topo10
  integer :: index_x2g_Fgss_qice10
  integer :: nflds_x2g

  ! glc -> drv

  integer :: index_g2x_Sg_frac01
  integer :: index_g2x_Sg_topo01
  integer :: index_g2x_Fsgg_rofi01
  integer :: index_g2x_Fsgg_rofl01
  integer :: index_g2x_Fsgg_hflx01
  integer :: index_g2x_Sg_frac02
  integer :: index_g2x_Sg_topo02
  integer :: index_g2x_Fsgg_rofi02
  integer :: index_g2x_Fsgg_rofl02
  integer :: index_g2x_Fsgg_hflx02
  integer :: index_g2x_Sg_frac03
  integer :: index_g2x_Sg_topo03
  integer :: index_g2x_Fsgg_rofi03
  integer :: index_g2x_Fsgg_rofl03
  integer :: index_g2x_Fsgg_hflx03
  integer :: index_g2x_Sg_frac04
  integer :: index_g2x_Sg_topo04
  integer :: index_g2x_Fsgg_rofi04
  integer :: index_g2x_Fsgg_rofl04
  integer :: index_g2x_Fsgg_hflx04
  integer :: index_g2x_Sg_frac05
  integer :: index_g2x_Sg_topo05
  integer :: index_g2x_Fsgg_rofi05
  integer :: index_g2x_Fsgg_rofl05
  integer :: index_g2x_Fsgg_hflx05
  integer :: index_g2x_Sg_frac06
  integer :: index_g2x_Sg_topo06
  integer :: index_g2x_Fsgg_rofi06
  integer :: index_g2x_Fsgg_rofl06
  integer :: index_g2x_Fsgg_hflx06
  integer :: index_g2x_Sg_frac07
  integer :: index_g2x_Sg_topo07
  integer :: index_g2x_Fsgg_rofi07
  integer :: index_g2x_Fsgg_rofl07
  integer :: index_g2x_Fsgg_hflx07
  integer :: index_g2x_Sg_frac08
  integer :: index_g2x_Sg_topo08
  integer :: index_g2x_Fsgg_rofi08
  integer :: index_g2x_Fsgg_rofl08
  integer :: index_g2x_Fsgg_hflx08
  integer :: index_g2x_Sg_frac09
  integer :: index_g2x_Sg_topo09
  integer :: index_g2x_Fsgg_rofi09
  integer :: index_g2x_Fsgg_rofl09
  integer :: index_g2x_Fsgg_hflx09
  integer :: index_g2x_Sg_frac10
  integer :: index_g2x_Sg_topo10
  integer :: index_g2x_Fsgg_rofi10
  integer :: index_g2x_Fsgg_rofl10
  integer :: index_g2x_Fsgg_hflx10
  integer :: nflds_g2x	

  ! sno -> drv

  integer :: index_s2x_Ss_tsrf01
  integer :: index_s2x_Ss_topo01
  integer :: index_s2x_Fgss_qice01
  integer :: index_s2x_Ss_tsrf02
  integer :: index_s2x_Ss_topo02
  integer :: index_s2x_Fgss_qice02
  integer :: index_s2x_Ss_tsrf03
  integer :: index_s2x_Ss_topo03
  integer :: index_s2x_Fgss_qice03
  integer :: index_s2x_Ss_tsrf04
  integer :: index_s2x_Ss_topo04
  integer :: index_s2x_Fgss_qice04
  integer :: index_s2x_Ss_tsrf05
  integer :: index_s2x_Ss_topo05
  integer :: index_s2x_Fgss_qice05
  integer :: index_s2x_Ss_tsrf06
  integer :: index_s2x_Ss_topo06
  integer :: index_s2x_Fgss_qice06
  integer :: index_s2x_Ss_tsrf07
  integer :: index_s2x_Ss_topo07
  integer :: index_s2x_Fgss_qice07
  integer :: index_s2x_Ss_tsrf08
  integer :: index_s2x_Ss_topo08
  integer :: index_s2x_Fgss_qice08
  integer :: index_s2x_Ss_tsrf09
  integer :: index_s2x_Ss_topo09
  integer :: index_s2x_Fgss_qice09
  integer :: index_s2x_Ss_tsrf10
  integer :: index_s2x_Ss_topo10
  integer :: index_s2x_Fgss_qice10
  integer :: nflds_s2x

  ! drv -> sno

  integer :: index_x2s_Sg_frac01
  integer :: index_x2s_Sg_topo01
  integer :: index_x2s_Fsgg_rofi01
  integer :: index_x2s_Fsgg_rofl01
  integer :: index_x2s_Fsgg_hflx01
  integer :: index_x2s_Sg_frac02
  integer :: index_x2s_Sg_topo02
  integer :: index_x2s_Fsgg_rofi02
  integer :: index_x2s_Fsgg_rofl02
  integer :: index_x2s_Fsgg_hflx02
  integer :: index_x2s_Sg_frac03
  integer :: index_x2s_Sg_topo03
  integer :: index_x2s_Fsgg_rofi03
  integer :: index_x2s_Fsgg_rofl03
  integer :: index_x2s_Fsgg_hflx03
  integer :: index_x2s_Sg_frac04
  integer :: index_x2s_Sg_topo04
  integer :: index_x2s_Fsgg_rofi04
  integer :: index_x2s_Fsgg_rofl04
  integer :: index_x2s_Fsgg_hflx04
  integer :: index_x2s_Sg_frac05
  integer :: index_x2s_Sg_topo05
  integer :: index_x2s_Fsgg_rofi05
  integer :: index_x2s_Fsgg_rofl05
  integer :: index_x2s_Fsgg_hflx05
  integer :: index_x2s_Sg_frac06
  integer :: index_x2s_Sg_topo06
  integer :: index_x2s_Fsgg_rofi06
  integer :: index_x2s_Fsgg_rofl06
  integer :: index_x2s_Fsgg_hflx06
  integer :: index_x2s_Sg_frac07
  integer :: index_x2s_Sg_topo07
  integer :: index_x2s_Fsgg_rofi07
  integer :: index_x2s_Fsgg_rofl07
  integer :: index_x2s_Fsgg_hflx07
  integer :: index_x2s_Sg_frac08
  integer :: index_x2s_Sg_topo08
  integer :: index_x2s_Fsgg_rofi08
  integer :: index_x2s_Fsgg_rofl08
  integer :: index_x2s_Fsgg_hflx08
  integer :: index_x2s_Sg_frac09
  integer :: index_x2s_Sg_topo09
  integer :: index_x2s_Fsgg_rofi09
  integer :: index_x2s_Fsgg_rofl09
  integer :: index_x2s_Fsgg_hflx09
  integer :: index_x2s_Sg_frac10
  integer :: index_x2s_Sg_topo10
  integer :: index_x2s_Fsgg_rofi10
  integer :: index_x2s_Fsgg_rofl10
  integer :: index_x2s_Fsgg_hflx10
  integer :: nflds_x2s	

! !PUBLIC MEMBER FUNCTIONS:

   public :: seq_flds_indices_set 

!----------------------------------------------------------------------------
contains
!----------------------------------------------------------------------------

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: seq_flds_indices_set
!
! !DESCRIPTION:
! Sets the values for all component state/flux fields
!
! !REVISION HISTORY:
! M. Vertenstein - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine seq_flds_indices_set

! !USES:

! !INPUT/OUTPUT PARAMETERS:

!EOP
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

    !-------------------------------------------------------------
    ! domain
    !-------------------------------------------------------------

    nflds_dom = shr_string_listGetNum(seq_flds_dom_fields)

    !-------------------------------------------------------------
    ! atm -> drv
    !-------------------------------------------------------------

    index_a2x_Sa_z          = shr_string_listGetIndexF(seq_flds_a2x_fields,'Sa_z')
    if (index_a2x_Sa_z == 0) call shr_sys_abort('index_a2x_Sa_z is zero')

    index_a2x_Sa_u          = shr_string_listGetIndexF(seq_flds_a2x_fields,'Sa_u')
    if (index_a2x_Sa_u == 0) call shr_sys_abort('index_a2x_Sa_u is zero')

    index_a2x_Sa_v          = shr_string_listGetIndexF(seq_flds_a2x_fields,'Sa_v')
    if (index_a2x_Sa_v == 0) call shr_sys_abort('index_a2x_Sa_v is zero')

    index_a2x_Sa_tbot       = shr_string_listGetIndexF(seq_flds_a2x_fields,'Sa_tbot')
    if (index_a2x_Sa_tbot == 0) call shr_sys_abort('index_a2x_Sa_tbot is zero')

    index_a2x_Sa_ptem       = shr_string_listGetIndexF(seq_flds_a2x_fields,'Sa_ptem')
    if (index_a2x_Sa_ptem == 0) call shr_sys_abort('index_a2x_Sa_ptem is zero')

    index_a2x_Sa_pbot       = shr_string_listGetIndexF(seq_flds_a2x_fields,'Sa_pbot')
    if (index_a2x_Sa_pbot == 0) call shr_sys_abort('index_a2x_Sa_pbot is zero')

    index_a2x_Sa_pslv       = shr_string_listGetIndexF(seq_flds_a2x_fields,'Sa_pslv')
    if (index_a2x_Sa_pslv == 0) call shr_sys_abort('index_a2x_Sa_pslv is zero')

    index_a2x_Sa_shum       = shr_string_listGetIndexF(seq_flds_a2x_fields,'Sa_shum')
    if (index_a2x_Sa_shum == 0) call shr_sys_abort('index_a2x_Sa_shum is zero')

    index_a2x_Sa_dens       = shr_string_listGetIndexF(seq_flds_a2x_fields,'Sa_dens')
    if (index_a2x_Sa_dens == 0) call shr_sys_abort('index_a2x_Sa_dens is zero')

    index_a2x_Faxa_swnet    = shr_string_listGetIndexF(seq_flds_a2x_fields,'Faxa_swnet')
    if (index_a2x_Faxa_swnet == 0) call shr_sys_abort('index_a2x_Faxa_swnet is zero')

    index_a2x_Faxa_lwdn     = shr_string_listGetIndexF(seq_flds_a2x_fields,'Faxa_lwdn')
    if (index_a2x_Faxa_lwdn == 0) call shr_sys_abort('index_a2x_Faxa_lwdn is zero')

    index_a2x_Faxa_rainc    = shr_string_listGetIndexF(seq_flds_a2x_fields,'Faxa_rainc')
    if (index_a2x_Faxa_rainc == 0) call shr_sys_abort('index_a2x_Faxa_rainc is zero')

    index_a2x_Faxa_rainl    = shr_string_listGetIndexF(seq_flds_a2x_fields,'Faxa_rainl')
    if (index_a2x_Faxa_rainl == 0) call shr_sys_abort('index_a2x_Faxa_rainl is zero')

    index_a2x_Faxa_snowc    = shr_string_listGetIndexF(seq_flds_a2x_fields,'Faxa_snowc')
    if (index_a2x_Faxa_snowc == 0) call shr_sys_abort('index_a2x_Faxa_snowc is zero')

    index_a2x_Faxa_snowl    = shr_string_listGetIndexF(seq_flds_a2x_fields,'Faxa_snowl')
    if (index_a2x_Faxa_snowl == 0) call shr_sys_abort('index_a2x_Faxa_snowl is zero')

    index_a2x_Faxa_swndr    = shr_string_listGetIndexF(seq_flds_a2x_fields,'Faxa_swndr')
    if (index_a2x_Faxa_swndr == 0) call shr_sys_abort('index_a2x_Faxa_swndr is zero')

    index_a2x_Faxa_swvdr    = shr_string_listGetIndexF(seq_flds_a2x_fields,'Faxa_swvdr')
    if (index_a2x_Faxa_swvdr == 0) call shr_sys_abort('index_a2x_Faxa_swvdr is zero')

    index_a2x_Faxa_swndf    = shr_string_listGetIndexF(seq_flds_a2x_fields,'Faxa_swndf')
    if (index_a2x_Faxa_swndf == 0) call shr_sys_abort('index_a2x_Faxa_swndf is zero')

    index_a2x_Faxa_swvdf    = shr_string_listGetIndexF(seq_flds_a2x_fields,'Faxa_swvdf')
    if (index_a2x_Faxa_swvdf == 0) call shr_sys_abort('index_a2x_Faxa_swvdf is zero')

    index_a2x_Faxa_bcphidry = shr_string_listGetIndexF(seq_flds_a2x_fields,'Faxa_bcphidry')
    if (index_a2x_Faxa_bcphidry == 0) call shr_sys_abort('index_a2x_Faxa_bcphidry is zero')

    index_a2x_Faxa_bcphodry = shr_string_listGetIndexF(seq_flds_a2x_fields,'Faxa_bcphodry')
    if (index_a2x_Faxa_bcphodry == 0) call shr_sys_abort('index_a2x_Faxa_bcphodry is zero')

    index_a2x_Faxa_bcphiwet = shr_string_listGetIndexF(seq_flds_a2x_fields,'Faxa_bcphiwet')
    if (index_a2x_Faxa_bcphiwet== 0) call shr_sys_abort('index_a2x_Faxa_bcphiwet is zero')

    index_a2x_Faxa_ocphidry = shr_string_listGetIndexF(seq_flds_a2x_fields,'Faxa_ocphidry')
    if (index_a2x_Faxa_ocphidry == 0) call shr_sys_abort('index_a2x_Faxa_ocphidry is zero')

    index_a2x_Faxa_ocphodry = shr_string_listGetIndexF(seq_flds_a2x_fields,'Faxa_ocphodry')
    if (index_a2x_Faxa_ocphodry == 0) call shr_sys_abort('index_a2x_Faxa_ocphodry is zero')

    index_a2x_Faxa_ocphiwet = shr_string_listGetIndexF(seq_flds_a2x_fields,'Faxa_ocphiwet')
    if (index_a2x_Faxa_ocphiwet == 0) call shr_sys_abort('index_a2x_Faxa_ocphiwet is zero')

    index_a2x_Faxa_dstdry1  = shr_string_listGetIndexF(seq_flds_a2x_fields,'Faxa_dstdry1')
    if (index_a2x_Faxa_dstdry1 == 0) call shr_sys_abort('index_a2x_Faxa_dstdry1 is zero')

    index_a2x_Faxa_dstdry2  = shr_string_listGetIndexF(seq_flds_a2x_fields,'Faxa_dstdry2')
    if (index_a2x_Faxa_dstdry2 == 0) call shr_sys_abort('index_a2x_Faxa_dstdry2 is zero')

    index_a2x_Faxa_dstdry3  = shr_string_listGetIndexF(seq_flds_a2x_fields,'Faxa_dstdry3')
    if (index_a2x_Faxa_dstdry3 == 0) call shr_sys_abort('index_a2x_Faxa_dstdry3 is zero')

    index_a2x_Faxa_dstdry4  = shr_string_listGetIndexF(seq_flds_a2x_fields,'Faxa_dstdry4')
    if (index_a2x_Faxa_dstdry4 == 0) call shr_sys_abort('index_a2x_Faxa_dstdry4 is zero')

    index_a2x_Faxa_dstwet1  = shr_string_listGetIndexF(seq_flds_a2x_fields,'Faxa_dstwet1')
    if (index_a2x_Faxa_dstwet1 == 0) call shr_sys_abort('index_a2x_Faxa_dstwet1 is zero')

    index_a2x_Faxa_dstwet2  = shr_string_listGetIndexF(seq_flds_a2x_fields,'Faxa_dstwet2')
    if (index_a2x_Faxa_dstwet2 == 0) call shr_sys_abort('index_a2x_Faxa_dstwet2 is zero')

    index_a2x_Faxa_dstwet3  = shr_string_listGetIndexF(seq_flds_a2x_fields,'Faxa_dstwet3')
    if (index_a2x_Faxa_dstwet3 == 0) call shr_sys_abort('index_a2x_Faxa_dstwet3 is zero')

    index_a2x_Faxa_dstwet4  = shr_string_listGetIndexF(seq_flds_a2x_fields,'Faxa_dstwet4')
    if (index_a2x_Faxa_dstwet4 == 0) call shr_sys_abort('index_a2x_Faxa_dstwet4 is zero')

    ! Optional fields

    index_a2x_Sa_co2prog    = shr_string_listGetIndexF(seq_flds_a2x_fields,'Sa_co2prog')
    index_a2x_Sa_co2diag    = shr_string_listGetIndexF(seq_flds_a2x_fields,'Sa_co2diag')

    nflds_a2x = shr_string_listGetNum(seq_flds_a2x_fields)

    !-------------------------------------------------------------
    ! drv -> atm
    !-------------------------------------------------------------

    index_x2a_Sx_avsdr      = shr_string_listGetIndexF(seq_flds_x2a_fields,'Sx_avsdr')
    if (index_x2a_Sx_avsdr == 0) call shr_sys_abort('index_x2a_Sx_avsdr is zero')

    index_x2a_Sx_anidr      = shr_string_listGetIndexF(seq_flds_x2a_fields,'Sx_anidr')
    if (index_x2a_Sx_anidr == 0) call shr_sys_abort('index_x2a_Sx_anidr is zero')

    index_x2a_Sx_avsdf      = shr_string_listGetIndexF(seq_flds_x2a_fields,'Sx_avsdf')
    if (index_x2a_Sx_avsdf == 0) call shr_sys_abort('index_x2a_Sx_avsdf is zero')

    index_x2a_Sx_anidf      = shr_string_listGetIndexF(seq_flds_x2a_fields,'Sx_anidf')
    if (index_x2a_Sx_anidf == 0) call shr_sys_abort('index_x2a_Sx_anidf is zero')

    index_x2a_Sx_t          = shr_string_listGetIndexF(seq_flds_x2a_fields,'Sx_t')
    if (index_x2a_Sx_t == 0) call shr_sys_abort('index_x2a_Sx_t is zero')

    index_x2a_So_t          = shr_string_listGetIndexF(seq_flds_x2a_fields,'So_t')
    if (index_x2a_So_t == 0) call shr_sys_abort('index_x2a_So_t is zero')

    index_x2a_Sl_snowh      = shr_string_listGetIndexF(seq_flds_x2a_fields,'Sl_snowh')
    if (index_x2a_Sl_snowh == 0) call shr_sys_abort('index_x2a_Sl_snowh is zero')

    index_x2a_Si_snowh      = shr_string_listGetIndexF(seq_flds_x2a_fields,'Si_snowh')
    if (index_x2a_Si_snowh == 0) call shr_sys_abort('index_x2a_Si_snowh is zero')

    index_x2a_Sx_tref       = shr_string_listGetIndexF(seq_flds_x2a_fields,'Sx_tref')
    if (index_x2a_Sx_tref == 0) call shr_sys_abort('index_x2a_Sx_tref is zero')

    index_x2a_Sx_qref       = shr_string_listGetIndexF(seq_flds_x2a_fields,'Sx_qref')
    if (index_x2a_Sx_qref == 0) call shr_sys_abort('index_x2a_Sx_qref is zero')

    index_x2a_Sx_ifrac      = shr_string_listGetIndexF(seq_flds_x2a_fields,'Sx_ifrac')
    if (index_x2a_Sx_ifrac == 0) call shr_sys_abort('index_x2a_Sx_ifrac is zero')

    index_x2a_Sx_ofrac      = shr_string_listGetIndexF(seq_flds_x2a_fields,'Sx_ofrac')
    if (index_x2a_Sx_ofrac == 0) call shr_sys_abort('index_x2a_Sx_ofrac is zero')

    index_x2a_Sx_lfrac      = shr_string_listGetIndexF(seq_flds_x2a_fields,'Sx_lfrac')
    if (index_x2a_Sx_lfrac == 0) call shr_sys_abort('index_x2a_Sx_lfrac is zero')

    index_x2a_Sx_u10       = shr_string_listGetIndexF(seq_flds_x2a_fields,'Sx_u10')
    if (index_x2a_Sx_u10 == 0) call shr_sys_abort('index_x2a_Sx_u10 is zero')

    index_x2a_Faxx_taux     = shr_string_listGetIndexF(seq_flds_x2a_fields,'Faxx_taux')
    if (index_x2a_Faxx_taux == 0) call shr_sys_abort('index_x2a_Faxx_taux is zero')

    index_x2a_Faxx_tauy     = shr_string_listGetIndexF(seq_flds_x2a_fields,'Faxx_tauy')
    if (index_x2a_Faxx_tauy == 0) call shr_sys_abort('index_x2a_Faxx_tauy is zero')

    index_x2a_Faxx_lat      = shr_string_listGetIndexF(seq_flds_x2a_fields,'Faxx_lat')
    if (index_x2a_Faxx_lat == 0) call shr_sys_abort('index_x2a_Faxx_lat is zero')

    index_x2a_Faxx_sen      = shr_string_listGetIndexF(seq_flds_x2a_fields,'Faxx_sen')
    if (index_x2a_Faxx_sen == 0) call shr_sys_abort('index_x2a_Faxx_sen is zero')

    index_x2a_Faxx_lwup     = shr_string_listGetIndexF(seq_flds_x2a_fields,'Faxx_lwup')
    if (index_x2a_Faxx_lwup == 0) call shr_sys_abort('index_x2a_Faxx_lwup is zero')

    index_x2a_Faxx_evap     = shr_string_listGetIndexF(seq_flds_x2a_fields,'Faxx_evap')
    if (index_x2a_Faxx_evap == 0) call shr_sys_abort('index_x2a_Faxx_evap is zero')

    ! fields needed to calculate water isotopes to ocean evaporation processes

    index_x2a_So_ustar     = shr_string_listGetIndexF(seq_flds_x2a_fields,'So_ustar')
    if (index_x2a_So_ustar == 0) call shr_sys_abort('index_x2a_So_ustar is zero')
    
    index_x2a_So_re     = shr_string_listGetIndexF(seq_flds_x2a_fields,'So_re')
    if (index_x2a_So_re == 0) call shr_sys_abort('index_x2a_So_re is zero')

    index_x2a_So_ssq     = shr_string_listGetIndexF(seq_flds_x2a_fields,'So_ssq')
    if (index_x2a_So_ssq == 0) call shr_sys_abort('index_x2a_So_ssq is zero')

    ! Optional fields

    index_x2a_Si_avsdr      = shr_string_listGetIndexF(seq_flds_x2a_fields,'Si_avsdr')
    index_x2a_si_anidr      = shr_string_listGetIndexF(seq_flds_x2a_fields,'Si_anidr')
    index_x2a_si_avsdf      = shr_string_listGetIndexF(seq_flds_x2a_fields,'Si_avsdf')
    index_x2a_si_anidf      = shr_string_listGetIndexF(seq_flds_x2a_fields,'Si_anidf')
    index_x2a_So_avsdr      = shr_string_listGetIndexF(seq_flds_x2a_fields,'So_avsdr')
    index_x2a_So_anidr      = shr_string_listGetIndexF(seq_flds_x2a_fields,'So_anidr')
    index_x2a_So_avsdf      = shr_string_listGetIndexF(seq_flds_x2a_fields,'So_avsdf')
    index_x2a_So_anidf      = shr_string_listGetIndexF(seq_flds_x2a_fields,'So_anidf')
    index_x2a_Faox_lat      = shr_string_listGetIndexF(seq_flds_x2a_fields,'Faox_lat')
    index_x2a_Faox_sen      = shr_string_listGetIndexF(seq_flds_x2a_fields,'Faox_sen')
    index_x2a_Faox_lwup     = shr_string_listGetIndexF(seq_flds_x2a_fields,'Faox_lwup')
    index_x2a_Faii_lat      = shr_string_listGetIndexF(seq_flds_x2a_fields,'Faii_lat')
    index_x2a_Faii_sen      = shr_string_listGetIndexF(seq_flds_x2a_fields,'Faii_sen')
    index_x2a_Faii_lwup     = shr_string_listGetIndexF(seq_flds_x2a_fields,'Faii_lwup')
    index_x2a_Sl_fv         = shr_string_listGetIndexF(seq_flds_x2a_fields,'Sl_fv')
    index_x2a_Sl_ram1       = shr_string_listGetIndexF(seq_flds_x2a_fields,'Sl_ram1')
    index_x2a_Fall_flxdst1  = shr_string_listGetIndexF(seq_flds_x2a_fields,'Fall_flxdst1' )
    index_x2a_Fall_flxdst2  = shr_string_listGetIndexF(seq_flds_x2a_fields,'Fall_flxdst2' )
    index_x2a_Fall_flxdst3  = shr_string_listGetIndexF(seq_flds_x2a_fields,'Fall_flxdst3' )
    index_x2a_Fall_flxdst4  = shr_string_listGetIndexF(seq_flds_x2a_fields,'Fall_flxdst4' )
    index_x2a_Faxx_fco2_lnd = shr_string_listGetIndexF(seq_flds_x2a_fields,'Faxx_fco2_lnd')
    index_x2a_Faxx_fco2_ocn = shr_string_listGetIndexF(seq_flds_x2a_fields,'Faxx_fco2_ocn')
    index_x2a_Faxx_fdms     = shr_string_listGetIndexF(seq_flds_x2a_fields,'Faxx_fdms'    )

   ! dry deposition velocities
    if ( lnd_drydep )then
       index_x2a_Sx_ddvel   = shr_string_listGetIndexF(seq_flds_x2a_fields, trim(drydep_fields_token))
    else
       index_x2a_Sx_ddvel   = 0
    end if

    index_x2a_Faxx_flxvoc1  = shr_string_listGetIndexF(seq_flds_x2a_fields,'Fall_flxvoc1')
    index_x2a_Faxx_flxvoc2  = shr_string_listGetIndexF(seq_flds_x2a_fields,'Fall_flxvoc2')
    index_x2a_Faxx_flxvoc3  = shr_string_listGetIndexF(seq_flds_x2a_fields,'Fall_flxvoc3')
    index_x2a_Faxx_flxvoc4  = shr_string_listGetIndexF(seq_flds_x2a_fields,'Fall_flxvoc4')
    index_x2a_Faxx_flxvoc5  = shr_string_listGetIndexF(seq_flds_x2a_fields,'Fall_flxvoc5')

    nflds_x2a = shr_string_listGetNum(seq_flds_x2a_fields)

    !-------------------------------------------------------------
    ! ocn -> drv
    !-------------------------------------------------------------

    index_o2x_So_t          = shr_string_listGetIndexF(seq_flds_o2x_fields,'So_t')
    if (index_o2x_So_t == 0) call shr_sys_abort('index_o2x_So_t is zero')

    index_o2x_So_u          = shr_string_listGetIndexF(seq_flds_o2x_fields,'So_u')
    if (index_o2x_So_u == 0) call shr_sys_abort('index_o2x_So_u is zero')

    index_o2x_So_v          = shr_string_listGetIndexF(seq_flds_o2x_fields,'So_v')
    if (index_o2x_So_v == 0) call shr_sys_abort('index_o2x_So_v is zero')

    index_o2x_So_s          = shr_string_listGetIndexF(seq_flds_o2x_fields,'So_s')
    if (index_o2x_So_s == 0) call shr_sys_abort('index_o2x_So_s is zero')

    index_o2x_So_dhdx          = shr_string_listGetIndexF(seq_flds_o2x_fields,'So_dhdx')
    if (index_o2x_So_dhdx == 0) call shr_sys_abort('index_o2x_So_dhdx is zero')

    index_o2x_So_dhdy          = shr_string_listGetIndexF(seq_flds_o2x_fields,'So_dhdy')
    if (index_o2x_So_dhdy == 0) call shr_sys_abort('index_o2x_So_dhdy is zero')

    index_o2x_Fioo_q        = shr_string_listGetIndexF(seq_flds_o2x_fields,'Fioo_q')
    if (index_o2x_Fioo_q == 0) call shr_sys_abort('index_o2x_Fioo_q is zero')

    index_o2x_Faoo_fco2     = shr_string_listGetIndexF(seq_flds_o2x_fields,'Faoo_fco2')
    index_o2x_Faoo_fdms     = shr_string_listGetIndexF(seq_flds_o2x_fields,'Faoo_fdms')

    nflds_o2x = shr_string_listGetNum(seq_flds_o2x_fields)

    !-------------------------------------------------------------
    ! hub atm/ocn fluxes/states
    !-------------------------------------------------------------

    index_xao_So_tref       = shr_string_listGetIndexF(seq_flds_xao_fields,'So_tref')
    if (index_xao_So_tref == 0) call shr_sys_abort('index_xao_So_tref is zero')

    index_xao_So_qref       = shr_string_listGetIndexF(seq_flds_xao_fields,'So_qref')
    if (index_xao_So_qref == 0) call shr_sys_abort('index_xao_So_qref is zero')

    index_xao_So_avsdr      = shr_string_listGetIndexF(seq_flds_xao_fields,'So_avsdr')
    if (index_xao_So_avsdr == 0) call shr_sys_abort('index_xao_So_avsdr is zero')

    index_xao_So_anidr      = shr_string_listGetIndexF(seq_flds_xao_fields,'So_anidr')
    if (index_xao_So_anidr == 0) call shr_sys_abort('index_xao_So_anidr is zero')

    index_xao_So_avsdf      = shr_string_listGetIndexF(seq_flds_xao_fields,'So_avsdf')
    if (index_xao_So_avsdf == 0) call shr_sys_abort('index_xao_So_avsdf is zero')

    index_xao_So_anidf      = shr_string_listGetIndexF(seq_flds_xao_fields,'So_anidf')
    if (index_xao_So_anidf == 0) call shr_sys_abort('index_xao_So_anidf is zero')

    index_xao_Sx_duu10n      = shr_string_listGetIndexF(seq_flds_xao_fields,'Sx_duu10n')
    if (index_xao_Sx_duu10n == 0) call shr_sys_abort('index_xao_Sx_duu10n is zero')

    index_xao_So_u10         = shr_string_listGetIndexF(seq_flds_xao_fields,'So_u10')
    if (index_xao_So_u10 == 0) call shr_sys_abort('index_xao_So_u10 is zero')

    index_xao_So_ustar      = shr_string_listGetIndexF(seq_flds_xao_fields,'So_ustar')
    if (index_xao_So_ustar == 0) call shr_sys_abort('index_xao_So_ustar is zero')

    index_xao_So_re      = shr_string_listGetIndexF(seq_flds_xao_fields,'So_re')
    if (index_xao_So_re == 0) call shr_sys_abort('index_xao_So_re is zero')

    index_xao_So_ssq      = shr_string_listGetIndexF(seq_flds_xao_fields,'So_ssq')
    if (index_xao_So_ssq == 0) call shr_sys_abort('index_xao_So_ssq is zero')

    index_xao_Faox_taux     = shr_string_listGetIndexF(seq_flds_xao_fields,'Faox_taux')
    if (index_xao_Faox_taux == 0) call shr_sys_abort('index_xao_Faox_taux is zero')

    index_xao_Faox_tauy     = shr_string_listGetIndexF(seq_flds_xao_fields,'Faox_tauy')
    if (index_xao_Faox_tauy == 0) call shr_sys_abort('index_xao_Faox_tauy is zero')

    index_xao_Faox_lat      = shr_string_listGetIndexF(seq_flds_xao_fields,'Faox_lat')
    if (index_xao_Faox_lat == 0) call shr_sys_abort('index_xao_Faox_lat is zero')

    index_xao_Faox_sen      = shr_string_listGetIndexF(seq_flds_xao_fields,'Faox_sen')
    if (index_xao_Faox_sen == 0) call shr_sys_abort('index_xao_Faox_sen is zero')

    index_xao_Faox_evap     = shr_string_listGetIndexF(seq_flds_xao_fields,'Faox_evap')
    if (index_xao_Faox_evap == 0) call shr_sys_abort('index_xao_Faox_evap is zero')

    index_xao_Faox_lwup     = shr_string_listGetIndexF(seq_flds_xao_fields,'Faox_lwup')
    if (index_xao_Faox_lwup == 0) call shr_sys_abort('index_xao_Faox_lwup is zero')

    nflds_xao = shr_string_listGetNum(seq_flds_xao_fields)

    !-------------------------------------------------------------
    ! drv -> ocn
    !-------------------------------------------------------------

    index_x2o_Si_ifrac      = shr_string_listGetIndexF(seq_flds_x2o_fields,'Si_ifrac')
    if (index_x2o_Si_ifrac == 0) call shr_sys_abort('index_x2o_Si_ifrac is zero')

    index_x2o_Sa_pslv      = shr_string_listGetIndexF(seq_flds_x2o_fields,'Sa_pslv')
    if (index_x2o_Sa_pslv == 0) call shr_sys_abort('index_x2o_Sa_pslv is zero')

    index_x2o_Sx_duu10n      = shr_string_listGetIndexF(seq_flds_x2o_fields,'Sx_duu10n')
    if (index_x2o_Sx_duu10n == 0) call shr_sys_abort('index_x2o_Sx_duu10n is zero')

    index_x2o_Foxx_tauy     = shr_string_listGetIndexF(seq_flds_x2o_fields,'Foxx_tauy')
    if (index_x2o_Foxx_tauy == 0) call shr_sys_abort('index_x2o_Foxx_tauy is zero')

    index_x2o_Foxx_taux     = shr_string_listGetIndexF(seq_flds_x2o_fields,'Foxx_taux')
    if (index_x2o_Foxx_taux == 0) call shr_sys_abort('index_x2o_Foxx_taux is zero')

    index_x2o_Foxx_swnet    = shr_string_listGetIndexF(seq_flds_x2o_fields,'Foxx_swnet')
    if (index_x2o_Foxx_swnet == 0) call shr_sys_abort('index_x2o_Foxx_swnet is zero')

    index_x2o_Foxx_lat      = shr_string_listGetIndexF(seq_flds_x2o_fields,'Foxx_lat')
    if (index_x2o_Foxx_lat == 0) call shr_sys_abort('index_x2o_Foxx_lat is zero')

    index_x2o_Foxx_sen      = shr_string_listGetIndexF(seq_flds_x2o_fields,'Foxx_sen')
    if (index_x2o_Foxx_sen == 0) call shr_sys_abort('index_x2o_Foxx_sen is zero')

    index_x2o_Foxx_lwdn     = shr_string_listGetIndexF(seq_flds_x2o_fields,'Foxx_lwdn')
    if (index_x2o_Foxx_lwdn == 0) call shr_sys_abort('index_x2o_Foxx_lwdn is zero')

    index_x2o_Foxx_lwup     = shr_string_listGetIndexF(seq_flds_x2o_fields,'Foxx_lwup')
    if (index_x2o_Foxx_lwup == 0) call shr_sys_abort('index_x2o_Foxx_lwup is zero')

    index_x2o_Foxx_melth    = shr_string_listGetIndexF(seq_flds_x2o_fields,'Foxx_melth')   
    if (index_x2o_Foxx_melth == 0) call shr_sys_abort('index_x2o_Foxx_melth is zero')

    index_x2o_Foxx_salt    = shr_string_listGetIndexF(seq_flds_x2o_fields,'Foxx_salt')   
    if (index_x2o_Foxx_salt == 0) call shr_sys_abort('index_x2o_Foxx_salt is zero')

    index_x2o_Foxx_prec    = shr_string_listGetIndexF(seq_flds_x2o_fields,'Foxx_prec')   
    if (index_x2o_Foxx_prec == 0) call shr_sys_abort('index_x2o_Foxx_prec is zero')

    index_x2o_Foxx_snow    = shr_string_listGetIndexF(seq_flds_x2o_fields,'Foxx_snow')   
    if (index_x2o_Foxx_snow == 0) call shr_sys_abort('index_x2o_Foxx_snow is zero')

    index_x2o_Foxx_rain    = shr_string_listGetIndexF(seq_flds_x2o_fields,'Foxx_rain')   
    if (index_x2o_Foxx_rain == 0) call shr_sys_abort('index_x2o_Foxx_rain is zero')

    index_x2o_Foxx_evap    = shr_string_listGetIndexF(seq_flds_x2o_fields,'Foxx_evap')
    if (index_x2o_Foxx_evap == 0) call shr_sys_abort('index_x2o_Foxx_evap is zero')

    index_x2o_Foxx_meltw    = shr_string_listGetIndexF(seq_flds_x2o_fields,'Foxx_meltw')
    if (index_x2o_Foxx_meltw == 0) call shr_sys_abort('index_x2o_Foxx_meltw is zero')

    index_x2o_Forr_roff    = shr_string_listGetIndexF(seq_flds_x2o_fields,'Forr_roff')
    if (index_x2o_Forr_roff == 0) call shr_sys_abort('index_x2o_Forr_roff is zero')

    index_x2o_Forr_ioff    = shr_string_listGetIndexF(seq_flds_x2o_fields,'Forr_ioff')
    if (index_x2o_Forr_ioff == 0) call shr_sys_abort('index_x2o_Forr_ioff is zero')

    index_x2o_Foxx_bcphidry = shr_string_listGetIndexF(seq_flds_x2o_fields,'Foxx_bcphidry')
    if (index_x2o_Foxx_bcphidry == 0) call shr_sys_abort('index_x2o_Foxx_bcphidry is zero')

    index_x2o_Foxx_bcphodry = shr_string_listGetIndexF(seq_flds_x2o_fields,'Foxx_bcphodry')
    if (index_x2o_Foxx_bcphodry == 0) call shr_sys_abort('index_x2o_Foxx_bcphodry is zero')

    index_x2o_Foxx_bcphiwet = shr_string_listGetIndexF(seq_flds_x2o_fields,'Foxx_bcphiwet')
    if (index_x2o_Foxx_bcphiwet== 0) call shr_sys_abort('index_x2o_Foxx_bcphiwet is zero')

    index_x2o_Foxx_ocphidry = shr_string_listGetIndexF(seq_flds_x2o_fields,'Foxx_ocphidry')
    if (index_x2o_Foxx_ocphidry == 0) call shr_sys_abort('index_x2o_Foxx_ocphidry is zero')

    index_x2o_Foxx_ocphodry = shr_string_listGetIndexF(seq_flds_x2o_fields,'Foxx_ocphodry')
    if (index_x2o_Foxx_ocphodry == 0) call shr_sys_abort('index_x2o_Foxx_ocphodry is zero')

    index_x2o_Foxx_ocphiwet = shr_string_listGetIndexF(seq_flds_x2o_fields,'Foxx_ocphiwet')
    if (index_x2o_Foxx_ocphiwet == 0) call shr_sys_abort('index_x2o_Foxx_ocphiwet is zero')

    index_x2o_Foxx_dstdry1  = shr_string_listGetIndexF(seq_flds_x2o_fields,'Foxx_dstdry1')
    if (index_x2o_Foxx_dstdry1 == 0) call shr_sys_abort('index_x2o_Foxx_dstdry1 is zero')

    index_x2o_Foxx_dstdry2  = shr_string_listGetIndexF(seq_flds_x2o_fields,'Foxx_dstdry2')
    if (index_x2o_Foxx_dstdry2 == 0) call shr_sys_abort('index_x2o_Foxx_dstdry2 is zero')

    index_x2o_Foxx_dstdry3  = shr_string_listGetIndexF(seq_flds_x2o_fields,'Foxx_dstdry3')
    if (index_x2o_Foxx_dstdry3 == 0) call shr_sys_abort('index_x2o_Foxx_dstdry3 is zero')

    index_x2o_Foxx_dstdry4  = shr_string_listGetIndexF(seq_flds_x2o_fields,'Foxx_dstdry4')
    if (index_x2o_Foxx_dstdry4 == 0) call shr_sys_abort('index_x2o_Foxx_dstdry4 is zero')

    index_x2o_Foxx_dstwet1  = shr_string_listGetIndexF(seq_flds_x2o_fields,'Foxx_dstwet1')
    if (index_x2o_Foxx_dstwet1 == 0) call shr_sys_abort('index_x2o_Foxx_dstwet1 is zero')

    index_x2o_Foxx_dstwet2  = shr_string_listGetIndexF(seq_flds_x2o_fields,'Foxx_dstwet2')
    if (index_x2o_Foxx_dstwet2 == 0) call shr_sys_abort('index_x2o_Foxx_dstwet2 is zero')

    index_x2o_Foxx_dstwet3  = shr_string_listGetIndexF(seq_flds_x2o_fields,'Foxx_dstwet3')
    if (index_x2o_Foxx_dstwet3 == 0) call shr_sys_abort('index_x2o_Foxx_dstwet3 is zero')

    index_x2o_Foxx_dstwet4  = shr_string_listGetIndexF(seq_flds_x2o_fields,'Foxx_dstwet4')
    if (index_x2o_Foxx_dstwet4 == 0) call shr_sys_abort('index_x2o_Foxx_dstwet4 is zero')

    ! Optional fields

    index_x2o_Sa_co2prog    = shr_string_listGetIndexF(seq_flds_x2o_fields,'Sa_co2prog')
    index_x2o_Sa_co2diag    = shr_string_listGetIndexF(seq_flds_x2o_fields,'Sa_co2diag')

    nflds_x2o = shr_string_listGetNum(seq_flds_x2o_fields)
    
    !-------------------------------------------------------------
    ! ice -> drv
    !-------------------------------------------------------------

    index_i2x_Si_t          = shr_string_listGetIndexF(seq_flds_i2x_fields,'Si_t')
    if (index_i2x_Si_t == 0) call shr_sys_abort('index_i2x_Si_t is zero')

    index_i2x_Si_tref       = shr_string_listGetIndexF(seq_flds_i2x_fields,'Si_tref')
    if (index_i2x_Si_tref == 0) call shr_sys_abort('index_i2x_Si_tref is zero')

    index_i2x_Si_qref       = shr_string_listGetIndexF(seq_flds_i2x_fields,'Si_qref')
    if (index_i2x_Si_qref == 0) call shr_sys_abort('index_i2x_Si_qref is zero')

    index_i2x_Si_ifrac      = shr_string_listGetIndexF(seq_flds_i2x_fields,'Si_ifrac')
    if (index_i2x_Si_ifrac == 0) call shr_sys_abort('index_i2x_Si_ifrac is zero')

    index_i2x_Si_avsdr      = shr_string_listGetIndexF(seq_flds_i2x_fields,'Si_avsdr')
    if (index_i2x_Si_avsdr == 0) call shr_sys_abort('index_i2x_Si_avsdr is zero')

    index_i2x_Si_anidr      = shr_string_listGetIndexF(seq_flds_i2x_fields,'Si_anidr')
    if (index_i2x_Si_anidr == 0) call shr_sys_abort('index_i2x_Si_anidr is zero')

    index_i2x_Si_avsdf      = shr_string_listGetIndexF(seq_flds_i2x_fields,'Si_avsdf')
    if (index_i2x_Si_avsdf == 0) call shr_sys_abort('index_i2x_Si_avsdf is zero')

    index_i2x_Si_anidf      = shr_string_listGetIndexF(seq_flds_i2x_fields,'Si_anidf')
    if (index_i2x_Si_anidf == 0) call shr_sys_abort('index_i2x_Si_anidf is zero')

    index_i2x_Si_snowh      = shr_string_listGetIndexF(seq_flds_i2x_fields,'Si_snowh')
    if (index_i2x_Si_snowh == 0) call shr_sys_abort('index_i2x_Si_snowh is zero')

    index_i2x_Si_u10       = shr_string_listGetIndexF(seq_flds_i2x_fields,'Si_u10')
    if (index_i2x_Si_u10 == 0) call shr_sys_abort('index_i2x_Si_u10 is zero')

    index_i2x_Faii_taux     = shr_string_listGetIndexF(seq_flds_i2x_fields,'Faii_taux')
    if (index_i2x_Faii_taux == 0) call shr_sys_abort('index_i2x_Faii_taux is zero')

    index_i2x_Faii_tauy     = shr_string_listGetIndexF(seq_flds_i2x_fields,'Faii_tauy')
    if (index_i2x_Faii_tauy == 0) call shr_sys_abort('index_i2x_Faii_tauy is zero')

    index_i2x_Faii_lat      = shr_string_listGetIndexF(seq_flds_i2x_fields,'Faii_lat')
    if (index_i2x_Faii_lat == 0) call shr_sys_abort('index_i2x_Faii_lat is zero')

    index_i2x_Faii_sen      = shr_string_listGetIndexF(seq_flds_i2x_fields,'Faii_sen')
    if (index_i2x_Faii_sen == 0) call shr_sys_abort('index_i2x_Faii_sen is zero')

    index_i2x_Faii_lwup     = shr_string_listGetIndexF(seq_flds_i2x_fields,'Faii_lwup')
    if (index_i2x_Faii_lwup == 0) call shr_sys_abort('index_i2x_Faii_lwup is zero')

    index_i2x_Faii_evap     = shr_string_listGetIndexF(seq_flds_i2x_fields,'Faii_evap')
    if (index_i2x_Faii_evap == 0) call shr_sys_abort('index_i2x_Faii_evap is zero')

    index_i2x_Faii_swnet    = shr_string_listGetIndexF(seq_flds_i2x_fields,'Faii_swnet')
    if (index_i2x_Faii_swnet == 0) call shr_sys_abort('index_i2x_Faii_swnet is zero')

    index_i2x_Fioi_swpen    = shr_string_listGetIndexF(seq_flds_i2x_fields,'Fioi_swpen')
    if (index_i2x_Fioi_swpen == 0) call shr_sys_abort('index_i2x_Fioi_swpen is zero')

    index_i2x_Fioi_melth    = shr_string_listGetIndexF(seq_flds_i2x_fields,'Fioi_melth')
    if (index_i2x_Fioi_melth == 0) call shr_sys_abort('index_i2x_Fioi_melth is zero')

    index_i2x_Fioi_meltw    = shr_string_listGetIndexF(seq_flds_i2x_fields,'Fioi_meltw')
    if (index_i2x_Fioi_meltw == 0) call shr_sys_abort('index_i2x_Fioi_meltw is zero')

    index_i2x_Fioi_salt     = shr_string_listGetIndexF(seq_flds_i2x_fields,'Fioi_salt')
    if (index_i2x_Fioi_salt == 0) call shr_sys_abort('index_i2x_Fioi_salt is zero')

    index_i2x_Fioi_taux     = shr_string_listGetIndexF(seq_flds_i2x_fields,'Fioi_taux')
    if (index_i2x_Fioi_taux == 0) call shr_sys_abort('index_i2x_Fioi_taux is zero')

    index_i2x_Fioi_tauy     = shr_string_listGetIndexF(seq_flds_i2x_fields,'Fioi_tauy')
    if (index_i2x_Fioi_tauy == 0) call shr_sys_abort('index_i2x_Fioi_tauy is zero')

    nflds_i2x               = shr_string_listGetNum(seq_flds_i2x_fields)

    !-------------------------------------------------------------
    ! drv -> ice
    !-------------------------------------------------------------

    index_x2i_So_t = shr_string_listGetIndexF(seq_flds_x2i_fields,'So_t')
    if (index_x2i_So_t == 0) call shr_sys_abort('index_x2i_So_t is zero')

    index_x2i_So_s = shr_string_listGetIndexF(seq_flds_x2i_fields,'So_s')
    if (index_x2i_So_s == 0) call shr_sys_abort('index_x2i_So_s is zero')

    index_x2i_So_u = shr_string_listGetIndexF(seq_flds_x2i_fields,'So_u')
    if (index_x2i_So_u == 0) call shr_sys_abort('index_x2i_So_u is zero')

    index_x2i_So_v = shr_string_listGetIndexF(seq_flds_x2i_fields,'So_v')
    if (index_x2i_So_v == 0) call shr_sys_abort('index_x2i_So_v is zero')

    index_x2i_Sa_z = shr_string_listGetIndexF(seq_flds_x2i_fields,'Sa_z')
    if (index_x2i_Sa_z == 0) call shr_sys_abort('index_x2i_Sa_z is zero')

    index_x2i_Sa_u = shr_string_listGetIndexF(seq_flds_x2i_fields,'Sa_u')
    if (index_x2i_Sa_u == 0) call shr_sys_abort('index_x2i_Sa_u is zero')

    index_x2i_Sa_v = shr_string_listGetIndexF(seq_flds_x2i_fields,'Sa_v')
    if (index_x2i_Sa_v == 0) call shr_sys_abort('index_x2i_Sa_v is zero')

    index_x2i_Sa_tbot = shr_string_listGetIndexF(seq_flds_x2i_fields,'Sa_tbot')
    if (index_x2i_Sa_tbot == 0) call shr_sys_abort('index_x2i_Sa_tbot is zero')

    index_x2i_Sa_ptem = shr_string_listGetIndexF(seq_flds_x2i_fields,'Sa_ptem')
    if (index_x2i_Sa_ptem == 0) call shr_sys_abort('index_x2i_Sa_ptem is zero')

    index_x2i_Sa_pbot = shr_string_listGetIndexF(seq_flds_x2i_fields,'Sa_pbot')
    if (index_x2i_Sa_pbot == 0) call shr_sys_abort('index_x2i_Sa_pbot is zero')

    index_x2i_Sa_shum = shr_string_listGetIndexF(seq_flds_x2i_fields,'Sa_shum')
    if (index_x2i_Sa_shum == 0) call shr_sys_abort('index_x2i_Sa_shum is zero')

    index_x2i_Sa_dens = shr_string_listGetIndexF(seq_flds_x2i_fields,'Sa_dens')
    if (index_x2i_Sa_dens == 0) call shr_sys_abort('index_x2i_Sa_dens is zero')

    index_x2i_So_dhdx = shr_string_listGetIndexF(seq_flds_x2i_fields,'So_dhdx')
    if (index_x2i_So_dhdx == 0) call shr_sys_abort('index_x2i_So_dhdx is zero')

    index_x2i_So_dhdy = shr_string_listGetIndexF(seq_flds_x2i_fields,'So_dhdy')
    if (index_x2i_So_dhdy == 0) call shr_sys_abort('index_x2i_So_dhdy is zero')

    index_x2i_Faxa_lwdn = shr_string_listGetIndexF(seq_flds_x2i_fields,'Faxa_lwdn')
    if (index_x2i_Faxa_lwdn == 0) call shr_sys_abort('index_x2i_Faxa_lwdn is zero')

    index_x2i_Faxa_rain = shr_string_listGetIndexF(seq_flds_x2i_fields,'Faxa_rain')
    if (index_x2i_Faxa_rain == 0) call shr_sys_abort('index_x2i_Faxa_rain is zero')

    index_x2i_Faxa_snow = shr_string_listGetIndexF(seq_flds_x2i_fields,'Faxa_snow')
    if (index_x2i_Faxa_snow == 0) call shr_sys_abort('index_x2i_Faxa_snow is zero')

    index_x2i_Faxa_swndr = shr_string_listGetIndexF(seq_flds_x2i_fields,'Faxa_swndr')
    if (index_x2i_Faxa_swndr == 0) call shr_sys_abort('index_x2i_Faxa_swndr is zero')

    index_x2i_Faxa_swvdr = shr_string_listGetIndexF(seq_flds_x2i_fields,'Faxa_swvdr')
    if (index_x2i_Faxa_swvdr == 0) call shr_sys_abort('index_x2i_Faxa_swvdr is zero')

    index_x2i_Faxa_swndf = shr_string_listGetIndexF(seq_flds_x2i_fields,'Faxa_swndf')
    if (index_x2i_Faxa_swndf == 0) call shr_sys_abort('index_x2i_Faxa_swndf is zero')

    index_x2i_Faxa_swvdf    = shr_string_listGetIndexF(seq_flds_x2i_fields,'Faxa_swvdf')
    if (index_x2i_Faxa_swvdf == 0) call shr_sys_abort('index_x2i_Faxa_swvdf is zero')

    index_x2i_Fioo_q = shr_string_listGetIndexF(seq_flds_x2i_fields,'Fioo_q')
    if (index_x2i_Fioo_q == 0) call shr_sys_abort('index_x2i_Fioo_q is zero')

    index_x2i_Faxa_bcphidry = shr_string_listGetIndexF(seq_flds_x2i_fields,'Faxa_bcphidry')
    if (index_x2i_Faxa_bcphidry == 0) call shr_sys_abort('index_x2i_Faxa_bcphidry is zero')

    index_x2i_Faxa_bcphodry = shr_string_listGetIndexF(seq_flds_x2i_fields,'Faxa_bcphodry')
    if (index_x2i_Faxa_bcphodry == 0) call shr_sys_abort('index_x2i_Faxa_bcphodry is zero')

    index_x2i_Faxa_bcphiwet = shr_string_listGetIndexF(seq_flds_x2i_fields,'Faxa_bcphiwet')
    if (index_x2i_Faxa_bcphiwet== 0) call shr_sys_abort('index_x2i_Faxa_bcphiwet is zero')

    index_x2i_Faxa_ocphidry = shr_string_listGetIndexF(seq_flds_x2i_fields,'Faxa_ocphidry')
    if (index_x2i_Faxa_ocphidry == 0) call shr_sys_abort('index_x2i_Faxa_ocphidry is zero')

    index_x2i_Faxa_ocphodry = shr_string_listGetIndexF(seq_flds_x2i_fields,'Faxa_ocphodry')
    if (index_x2i_Faxa_ocphodry == 0) call shr_sys_abort('index_x2i_Faxa_ocphodry is zero')

    index_x2i_Faxa_ocphiwet = shr_string_listGetIndexF(seq_flds_x2i_fields,'Faxa_ocphiwet')
    if (index_x2i_Faxa_ocphiwet == 0) call shr_sys_abort('index_x2i_Faxa_ocphiwet is zero')

    index_x2i_Faxa_dstdry1  = shr_string_listGetIndexF(seq_flds_x2i_fields,'Faxa_dstdry1')
    if (index_x2i_Faxa_dstdry1 == 0) call shr_sys_abort('index_x2i_Faxa_dstdry1 is zero')

    index_x2i_Faxa_dstdry2  = shr_string_listGetIndexF(seq_flds_x2i_fields,'Faxa_dstdry2')
    if (index_x2i_Faxa_dstdry2 == 0) call shr_sys_abort('index_x2i_Faxa_dstdry2 is zero')

    index_x2i_Faxa_dstdry3  = shr_string_listGetIndexF(seq_flds_x2i_fields,'Faxa_dstdry3')
    if (index_x2i_Faxa_dstdry3 == 0) call shr_sys_abort('index_x2i_Faxa_dstdry3 is zero')

    index_x2i_Faxa_dstdry4  = shr_string_listGetIndexF(seq_flds_x2i_fields,'Faxa_dstdry4')
    if (index_x2i_Faxa_dstdry4 == 0) call shr_sys_abort('index_x2i_Faxa_dstdry4 is zero')

    index_x2i_Faxa_dstwet1  = shr_string_listGetIndexF(seq_flds_x2i_fields,'Faxa_dstwet1')
    if (index_x2i_Faxa_dstwet1 == 0) call shr_sys_abort('index_x2i_Faxa_dstwet1 is zero')

    index_x2i_Faxa_dstwet2  = shr_string_listGetIndexF(seq_flds_x2i_fields,'Faxa_dstwet2')
    if (index_x2i_Faxa_dstwet2 == 0) call shr_sys_abort('index_x2i_Faxa_dstwet2 is zero')

    index_x2i_Faxa_dstwet3  = shr_string_listGetIndexF(seq_flds_x2i_fields,'Faxa_dstwet3')
    if (index_x2i_Faxa_dstwet3 == 0) call shr_sys_abort('index_x2i_Faxa_dstwet3 is zero')

    index_x2i_Faxa_dstwet4  = shr_string_listGetIndexF(seq_flds_x2i_fields,'Faxa_dstwet4')
    if (index_x2i_Faxa_dstwet4 == 0) call shr_sys_abort('index_x2i_Faxa_dstwet4 is zero')

    nflds_x2i = shr_string_listGetNum(seq_flds_x2i_fields)

    !-------------------------------------------------------------
    ! lnd -> drv 
    !-------------------------------------------------------------

    index_l2x_Sl_landfrac   = shr_string_listGetIndexF(seq_flds_l2x_fields,'Sl_landfrac')
    if (index_l2x_Sl_landfrac == 0) call shr_sys_abort('index_l2x_Sl_landfrac is zero')

    index_l2x_Sl_t          = shr_string_listGetIndexF(seq_flds_l2x_fields,'Sl_t')
    if (index_l2x_Sl_t == 0) call shr_sys_abort('index_l2x_Sl_t is zero')

    index_l2x_Sl_snowh      = shr_string_listGetIndexF(seq_flds_l2x_fields,'Sl_snowh')
    if (index_l2x_Sl_snowh == 0) call shr_sys_abort('index_l2x_Sl_snowh is zero')

    index_l2x_Sl_avsdr      = shr_string_listGetIndexF(seq_flds_l2x_fields,'Sl_avsdr')
    if (index_l2x_Sl_avsdr == 0) call shr_sys_abort('index_l2x_Sl_avsdr is zero')

    index_l2x_Sl_anidr      = shr_string_listGetIndexF(seq_flds_l2x_fields,'Sl_anidr')
    if (index_l2x_Sl_anidr == 0) call shr_sys_abort('index_l2x_Sl_anidr is zero')

    index_l2x_Sl_avsdf      = shr_string_listGetIndexF(seq_flds_l2x_fields,'Sl_avsdf')
    if (index_l2x_Sl_avsdf == 0) call shr_sys_abort('index_l2x_Sl_avsdf is zero')

    index_l2x_Sl_anidf      = shr_string_listGetIndexF(seq_flds_l2x_fields,'Sl_anidf')
    if (index_l2x_Sl_anidf == 0) call shr_sys_abort('index_l2x_Sl_anidf is zero')

    index_l2x_Sl_tref       = shr_string_listGetIndexF(seq_flds_l2x_fields,'Sl_tref')
    if (index_l2x_Sl_tref == 0) call shr_sys_abort('index_l2x_Sl_tref is zero')

    index_l2x_Sl_qref       = shr_string_listGetIndexF(seq_flds_l2x_fields,'Sl_qref')
    if (index_l2x_Sl_qref == 0) call shr_sys_abort('index_l2x_Sl_qref is zero')

    index_l2x_Sl_u10        = shr_string_listGetIndexF(seq_flds_l2x_fields,'Sl_u10')
    if (index_l2x_Sl_u10 == 0) call shr_sys_abort('index_l2x_Sl_u10 is zero')

    index_l2x_Fall_taux     = shr_string_listGetIndexF(seq_flds_l2x_fields,'Fall_taux')
    if (index_l2x_Fall_taux == 0) call shr_sys_abort('index_l2x_Fall_taux is zero')

    index_l2x_Fall_tauy     = shr_string_listGetIndexF(seq_flds_l2x_fields,'Fall_tauy')
    if (index_l2x_Fall_tauy == 0) call shr_sys_abort('index_l2x_Fall_tauy is zero')

    index_l2x_Fall_lat      = shr_string_listGetIndexF(seq_flds_l2x_fields,'Fall_lat')
    if (index_l2x_Fall_lat == 0) call shr_sys_abort('index_l2x_Fall_lat is zero')

    index_l2x_Fall_sen      = shr_string_listGetIndexF(seq_flds_l2x_fields,'Fall_sen')
    if (index_l2x_Fall_sen == 0) call shr_sys_abort('index_l2x_Fall_sen is zero')

    index_l2x_Fall_lwup     = shr_string_listGetIndexF(seq_flds_l2x_fields,'Fall_lwup')
    if (index_l2x_Fall_lwup == 0) call shr_sys_abort('index_l2x_Fall_lwup is zero')

    index_l2x_Fall_evap     = shr_string_listGetIndexF(seq_flds_l2x_fields,'Fall_evap')
    if (index_l2x_Fall_evap == 0) call shr_sys_abort('index_l2x_Fall_evap is zero')

    index_l2x_Fall_swnet    = shr_string_listGetIndexF(seq_flds_l2x_fields,'Fall_swnet')
    if (index_l2x_Fall_swnet == 0) call shr_sys_abort('index_l2x_Fall_swnet is zero')

    index_l2x_Sl_ram1      = shr_string_listGetIndexF(seq_flds_l2x_fields,'Sl_ram1')
    index_l2x_Sl_fv        = shr_string_listGetIndexF(seq_flds_l2x_fields,'Sl_fv')
    index_l2x_Fall_fco2_lnd= shr_string_listGetIndexF(seq_flds_l2x_fields,'Fall_fco2_lnd')
    index_l2x_Fall_flxdst1 = shr_string_listGetIndexF(seq_flds_l2x_fields,'Fall_flxdst1')
    index_l2x_Fall_flxdst2 = shr_string_listGetIndexF(seq_flds_l2x_fields,'Fall_flxdst2')
    index_l2x_Fall_flxdst3 = shr_string_listGetIndexF(seq_flds_l2x_fields,'Fall_flxdst3')
    index_l2x_Fall_flxdst4 = shr_string_listGetIndexF(seq_flds_l2x_fields,'Fall_flxdst4')

    ! dry deposition velocities
    if ( lnd_drydep )then
       index_l2x_Sl_ddvel  = shr_string_listGetIndexF(seq_flds_l2x_fields, trim(drydep_fields_token))
    else
       index_l2x_Sl_ddvel  = 0
    end if

    index_l2x_Fall_flxvoc1 = shr_string_listGetIndexF(seq_flds_l2x_fields,'Fall_flxvoc1')
    index_l2x_Fall_flxvoc2 = shr_string_listGetIndexF(seq_flds_l2x_fields,'Fall_flxvoc2')
    index_l2x_Fall_flxvoc3 = shr_string_listGetIndexF(seq_flds_l2x_fields,'Fall_flxvoc3')
    index_l2x_Fall_flxvoc4 = shr_string_listGetIndexF(seq_flds_l2x_fields,'Fall_flxvoc4')
    index_l2x_Fall_flxvoc5 = shr_string_listGetIndexF(seq_flds_l2x_fields,'Fall_flxvoc5')

    nflds_l2x              = shr_string_listGetNum(seq_flds_l2x_fields)

    !-------------------------------------------------------------
    ! runoff
    !-------------------------------------------------------------

    index_r2x_Forr_roff  = shr_string_listGetIndexF(seq_flds_r2x_fields,'Forr_roff')
    if (index_r2x_Forr_roff == 0) call shr_sys_abort('index_r2x_Forr_roff is zero')

    index_r2x_Forr_ioff  = shr_string_listGetIndexF(seq_flds_r2x_fields,'Forr_ioff')
    if (index_r2x_Forr_ioff == 0) call shr_sys_abort('index_r2x_Forr_ioff is zero')

    nflds_r2x           = shr_string_listGetNum(seq_flds_r2x_fields)

    !-------------------------------------------------------------
    ! drv -> lnd
    !-------------------------------------------------------------

    index_x2l_Sa_z = shr_string_listGetIndexF(seq_flds_x2l_fields,'Sa_z')
    if (index_x2l_Sa_z == 0) call shr_sys_abort('index_x2l_Sa_z is zero')

    index_x2l_Sa_u = shr_string_listGetIndexF(seq_flds_x2l_fields,'Sa_u')
    if (index_x2l_Sa_u == 0) call shr_sys_abort('index_x2l_Sa_u is zero')

    index_x2l_Sa_v = shr_string_listGetIndexF(seq_flds_x2l_fields,'Sa_v')
    if (index_x2l_Sa_v == 0) call shr_sys_abort('index_x2l_Sa_v is zero')

    index_x2l_Sa_ptem = shr_string_listGetIndexF(seq_flds_x2l_fields,'Sa_ptem')
    if (index_x2l_Sa_ptem == 0) call shr_sys_abort('index_x2l_Sa_ptem is zero')

    index_x2l_Sa_pbot = shr_string_listGetIndexF(seq_flds_x2l_fields,'Sa_pbot')
    if (index_x2l_Sa_pbot == 0) call shr_sys_abort('index_x2l_Sa_pbot is zero')

    index_x2l_Sa_tbot = shr_string_listGetIndexF(seq_flds_x2l_fields,'Sa_tbot')
    if (index_x2l_Sa_tbot == 0) call shr_sys_abort('index_x2l_Sa_tbot is zero')

    index_x2l_Sa_shum = shr_string_listGetIndexF(seq_flds_x2l_fields,'Sa_shum')
    if (index_x2l_Sa_shum == 0) call shr_sys_abort('index_x2l_Sa_shum is zero')

    index_x2l_Faxa_lwdn = shr_string_listGetIndexF(seq_flds_x2l_fields,'Faxa_lwdn')
    if (index_x2l_Faxa_lwdn == 0) call shr_sys_abort('index_x2l_Faxa_lwdn is zero')

    index_x2l_Faxa_rainc = shr_string_listGetIndexF(seq_flds_x2l_fields,'Faxa_rainc')
    if (index_x2l_Faxa_rainc == 0) call shr_sys_abort('index_x2l_Faxa_rainc is zero')

    index_x2l_Faxa_rainl = shr_string_listGetIndexF(seq_flds_x2l_fields,'Faxa_rainl')
    if (index_x2l_Faxa_rainl == 0) call shr_sys_abort('index_x2l_Faxa_rainl is zero')

    index_x2l_Faxa_snowc = shr_string_listGetIndexF(seq_flds_x2l_fields,'Faxa_snowc')
    if (index_x2l_Faxa_snowc == 0) call shr_sys_abort('index_x2l_Faxa_snowc is zero')

    index_x2l_Faxa_snowl = shr_string_listGetIndexF(seq_flds_x2l_fields,'Faxa_snowl')
    if (index_x2l_Faxa_snowl == 0) call shr_sys_abort('index_x2l_Faxa_snowl is zero')

    index_x2l_Faxa_swndr = shr_string_listGetIndexF(seq_flds_x2l_fields,'Faxa_swndr')
    if (index_x2l_Faxa_swndr == 0) call shr_sys_abort('index_x2l_Faxa_swndr is zero')

    index_x2l_Faxa_swvdr = shr_string_listGetIndexF(seq_flds_x2l_fields,'Faxa_swvdr')
    if (index_x2l_Faxa_swvdr == 0) call shr_sys_abort('index_x2l_Faxa_swvdr is zero')

    index_x2l_Faxa_swndf = shr_string_listGetIndexF(seq_flds_x2l_fields,'Faxa_swndf')
    if (index_x2l_Faxa_swndf == 0) call shr_sys_abort('index_x2l_Faxa_swndf is zero')

    index_x2l_Faxa_swvdf = shr_string_listGetIndexF(seq_flds_x2l_fields,'Faxa_swvdf')
    if (index_x2l_Faxa_swvdf == 0) call shr_sys_abort('index_x2l_Faxa_swvdf is zero')

    index_x2l_Faxa_bcphidry = shr_string_listGetIndexF(seq_flds_x2l_fields,'Faxa_bcphidry')
    if (index_x2l_Faxa_bcphidry == 0) call shr_sys_abort('index_x2l_Faxa_bcphidry is zero')

    index_x2l_Faxa_bcphodry = shr_string_listGetIndexF(seq_flds_x2l_fields,'Faxa_bcphodry')
    if (index_x2l_Faxa_bcphodry == 0) call shr_sys_abort('index_x2l_Faxa_bcphodry is zero')

    index_x2l_Faxa_bcphiwet = shr_string_listGetIndexF(seq_flds_x2l_fields,'Faxa_bcphiwet')
    if (index_x2l_Faxa_bcphiwet== 0) call shr_sys_abort('index_x2l_Faxa_bcphiwet is zero')

    index_x2l_Faxa_ocphidry = shr_string_listGetIndexF(seq_flds_x2l_fields,'Faxa_ocphidry')
    if (index_x2l_Faxa_ocphidry == 0) call shr_sys_abort('index_x2l_Faxa_ocphidry is zero')

    index_x2l_Faxa_ocphodry = shr_string_listGetIndexF(seq_flds_x2l_fields,'Faxa_ocphodry')
    if (index_x2l_Faxa_ocphodry == 0) call shr_sys_abort('index_x2l_Faxa_ocphodry is zero')

    index_x2l_Faxa_ocphiwet = shr_string_listGetIndexF(seq_flds_x2l_fields,'Faxa_ocphiwet')
    if (index_x2l_Faxa_ocphiwet == 0) call shr_sys_abort('index_x2l_Faxa_ocphiwet is zero')

    index_x2l_Faxa_dstdry1  = shr_string_listGetIndexF(seq_flds_x2l_fields,'Faxa_dstdry1')
    if (index_x2l_Faxa_dstdry1 == 0) call shr_sys_abort('index_x2l_Faxa_dstdry1 is zero')

    index_x2l_Faxa_dstdry2  = shr_string_listGetIndexF(seq_flds_x2l_fields,'Faxa_dstdry2')
    if (index_x2l_Faxa_dstdry2 == 0) call shr_sys_abort('index_x2l_Faxa_dstdry2 is zero')

    index_x2l_Faxa_dstdry3  = shr_string_listGetIndexF(seq_flds_x2l_fields,'Faxa_dstdry3')
    if (index_x2l_Faxa_dstdry3 == 0) call shr_sys_abort('index_x2l_Faxa_dstdry3 is zero')

    index_x2l_Faxa_dstdry4  = shr_string_listGetIndexF(seq_flds_x2l_fields,'Faxa_dstdry4')
    if (index_x2l_Faxa_dstdry4 == 0) call shr_sys_abort('index_x2l_Faxa_dstdry4 is zero')

    index_x2l_Faxa_dstwet1  = shr_string_listGetIndexF(seq_flds_x2l_fields,'Faxa_dstwet1')
    if (index_x2l_Faxa_dstwet1 == 0) call shr_sys_abort('index_x2l_Faxa_dstwet1 is zero')

    index_x2l_Faxa_dstwet2  = shr_string_listGetIndexF(seq_flds_x2l_fields,'Faxa_dstwet2')
    if (index_x2l_Faxa_dstwet2 == 0) call shr_sys_abort('index_x2l_Faxa_dstwet2 is zero')

    index_x2l_Faxa_dstwet3  = shr_string_listGetIndexF(seq_flds_x2l_fields,'Faxa_dstwet3')
    if (index_x2l_Faxa_dstwet3 == 0) call shr_sys_abort('index_x2l_Faxa_dstwet3 is zero')

    index_x2l_Faxa_dstwet4  = shr_string_listGetIndexF(seq_flds_x2l_fields,'Faxa_dstwet4')
    if (index_x2l_Faxa_dstwet4 == 0) call shr_sys_abort('index_x2l_Faxa_dstwet4 is zero')

    ! Optional fields

    index_x2l_Sa_co2prog = shr_string_listGetIndexF(seq_flds_x2l_fields,'Sa_co2prog')
    index_x2l_Sa_co2diag = shr_string_listGetIndexF(seq_flds_x2l_fields,'Sa_co2diag')

    nflds_x2l = shr_string_listGetNum(seq_flds_x2l_fields)

    !-------------------------------------------------------------
    ! glc -> drv
    !-------------------------------------------------------------

    index_g2x_Sg_frac01 = shr_string_listGetIndexF(seq_flds_g2x_fields,'Sg_frac01')
    index_g2x_Sg_topo01 = shr_string_listGetIndexF(seq_flds_g2x_fields,'Sg_topo01')
    index_g2x_Fsgg_rofi01 = shr_string_listGetIndexF(seq_flds_g2x_fields,'Fsgg_rofi01')
    index_g2x_Fsgg_rofl01 = shr_string_listGetIndexF(seq_flds_g2x_fields,'Fsgg_rofl01')
    index_g2x_Fsgg_hflx01 = shr_string_listGetIndexF(seq_flds_g2x_fields,'Fsgg_hflx01')
    index_g2x_Sg_frac02 = shr_string_listGetIndexF(seq_flds_g2x_fields,'Sg_frac02')
    index_g2x_Sg_topo02 = shr_string_listGetIndexF(seq_flds_g2x_fields,'Sg_topo02')
    index_g2x_Fsgg_rofi02 = shr_string_listGetIndexF(seq_flds_g2x_fields,'Fsgg_rofi02')
    index_g2x_Fsgg_rofl02 = shr_string_listGetIndexF(seq_flds_g2x_fields,'Fsgg_rofl02')
    index_g2x_Fsgg_hflx02 = shr_string_listGetIndexF(seq_flds_g2x_fields,'Fsgg_hflx02')
    index_g2x_Sg_frac03 = shr_string_listGetIndexF(seq_flds_g2x_fields,'Sg_frac03')
    index_g2x_Sg_topo03 = shr_string_listGetIndexF(seq_flds_g2x_fields,'Sg_topo03')
    index_g2x_Fsgg_rofi03 = shr_string_listGetIndexF(seq_flds_g2x_fields,'Fsgg_rofi03')
    index_g2x_Fsgg_rofl03 = shr_string_listGetIndexF(seq_flds_g2x_fields,'Fsgg_rofl03')
    index_g2x_Fsgg_hflx03 = shr_string_listGetIndexF(seq_flds_g2x_fields,'Fsgg_hflx03')
    index_g2x_Sg_frac04 = shr_string_listGetIndexF(seq_flds_g2x_fields,'Sg_frac04')
    index_g2x_Sg_topo04 = shr_string_listGetIndexF(seq_flds_g2x_fields,'Sg_topo04')
    index_g2x_Fsgg_rofi04 = shr_string_listGetIndexF(seq_flds_g2x_fields,'Fsgg_rofi04')
    index_g2x_Fsgg_rofl04 = shr_string_listGetIndexF(seq_flds_g2x_fields,'Fsgg_rofl04')
    index_g2x_Fsgg_hflx04 = shr_string_listGetIndexF(seq_flds_g2x_fields,'Fsgg_hflx04')
    index_g2x_Sg_frac05 = shr_string_listGetIndexF(seq_flds_g2x_fields,'Sg_frac05')
    index_g2x_Sg_topo05 = shr_string_listGetIndexF(seq_flds_g2x_fields,'Sg_topo05')
    index_g2x_Fsgg_rofi05 = shr_string_listGetIndexF(seq_flds_g2x_fields,'Fsgg_rofi05')
    index_g2x_Fsgg_rofl05 = shr_string_listGetIndexF(seq_flds_g2x_fields,'Fsgg_rofl05')
    index_g2x_Fsgg_hflx05 = shr_string_listGetIndexF(seq_flds_g2x_fields,'Fsgg_hflx05')
    index_g2x_Sg_frac06 = shr_string_listGetIndexF(seq_flds_g2x_fields,'Sg_frac06')
    index_g2x_Sg_topo06 = shr_string_listGetIndexF(seq_flds_g2x_fields,'Sg_topo06')
    index_g2x_Fsgg_rofi06 = shr_string_listGetIndexF(seq_flds_g2x_fields,'Fsgg_rofi06')
    index_g2x_Fsgg_rofl06 = shr_string_listGetIndexF(seq_flds_g2x_fields,'Fsgg_rofl06')
    index_g2x_Fsgg_hflx06 = shr_string_listGetIndexF(seq_flds_g2x_fields,'Fsgg_hflx06')
    index_g2x_Sg_frac07 = shr_string_listGetIndexF(seq_flds_g2x_fields,'Sg_frac07')
    index_g2x_Sg_topo07 = shr_string_listGetIndexF(seq_flds_g2x_fields,'Sg_topo07')
    index_g2x_Fsgg_rofi07 = shr_string_listGetIndexF(seq_flds_g2x_fields,'Fsgg_rofi07')
    index_g2x_Fsgg_rofl07 = shr_string_listGetIndexF(seq_flds_g2x_fields,'Fsgg_rofl07')
    index_g2x_Fsgg_hflx07 = shr_string_listGetIndexF(seq_flds_g2x_fields,'Fsgg_hflx07')
    index_g2x_Sg_frac08 = shr_string_listGetIndexF(seq_flds_g2x_fields,'Sg_frac08')
    index_g2x_Sg_topo08 = shr_string_listGetIndexF(seq_flds_g2x_fields,'Sg_topo08')
    index_g2x_Fsgg_rofi08 = shr_string_listGetIndexF(seq_flds_g2x_fields,'Fsgg_rofi08')
    index_g2x_Fsgg_rofl08 = shr_string_listGetIndexF(seq_flds_g2x_fields,'Fsgg_rofl08')
    index_g2x_Fsgg_hflx08 = shr_string_listGetIndexF(seq_flds_g2x_fields,'Fsgg_hflx08')
    index_g2x_Sg_frac09 = shr_string_listGetIndexF(seq_flds_g2x_fields,'Sg_frac09')
    index_g2x_Sg_topo09 = shr_string_listGetIndexF(seq_flds_g2x_fields,'Sg_topo09')
    index_g2x_Fsgg_rofi09 = shr_string_listGetIndexF(seq_flds_g2x_fields,'Fsgg_rofi09')
    index_g2x_Fsgg_rofl09 = shr_string_listGetIndexF(seq_flds_g2x_fields,'Fsgg_rofl09')
    index_g2x_Fsgg_hflx09 = shr_string_listGetIndexF(seq_flds_g2x_fields,'Fsgg_hflx09')
    index_g2x_Sg_frac10 = shr_string_listGetIndexF(seq_flds_g2x_fields,'Sg_frac10')
    index_g2x_Sg_topo10 = shr_string_listGetIndexF(seq_flds_g2x_fields,'Sg_topo10')
    index_g2x_Fsgg_rofi10 = shr_string_listGetIndexF(seq_flds_g2x_fields,'Fsgg_rofi10')
    index_g2x_Fsgg_rofl10 = shr_string_listGetIndexF(seq_flds_g2x_fields,'Fsgg_rofl10')
    index_g2x_Fsgg_hflx10 = shr_string_listGetIndexF(seq_flds_g2x_fields,'Fsgg_hflx10')

    nflds_g2x = shr_string_listGetNum(seq_flds_g2x_fields)

    !-------------------------------------------------------------
    ! drv -> glc
    !-------------------------------------------------------------

    index_x2g_Ss_tsrf01 = shr_string_listGetIndexF(seq_flds_x2g_fields,'Ss_tsrf01')
    index_x2g_Ss_topo01 = shr_string_listGetIndexF(seq_flds_x2g_fields,'Ss_topo01')
    index_x2g_Fgss_qice01 = shr_string_listGetIndexF(seq_flds_x2g_fields,'Fgss_qice01')
    index_x2g_Ss_tsrf02 = shr_string_listGetIndexF(seq_flds_x2g_fields,'Ss_tsrf02')
    index_x2g_Ss_topo02 = shr_string_listGetIndexF(seq_flds_x2g_fields,'Ss_topo02')
    index_x2g_Fgss_qice02 = shr_string_listGetIndexF(seq_flds_x2g_fields,'Fgss_qice02')
    index_x2g_Ss_tsrf03 = shr_string_listGetIndexF(seq_flds_x2g_fields,'Ss_tsrf03')
    index_x2g_Ss_topo03 = shr_string_listGetIndexF(seq_flds_x2g_fields,'Ss_topo03')
    index_x2g_Fgss_qice03 = shr_string_listGetIndexF(seq_flds_x2g_fields,'Fgss_qice03')
    index_x2g_Ss_tsrf04 = shr_string_listGetIndexF(seq_flds_x2g_fields,'Ss_tsrf04')
    index_x2g_Ss_topo04 = shr_string_listGetIndexF(seq_flds_x2g_fields,'Ss_topo04')
    index_x2g_Fgss_qice04 = shr_string_listGetIndexF(seq_flds_x2g_fields,'Fgss_qice04')
    index_x2g_Ss_tsrf05 = shr_string_listGetIndexF(seq_flds_x2g_fields,'Ss_tsrf05')
    index_x2g_Ss_topo05 = shr_string_listGetIndexF(seq_flds_x2g_fields,'Ss_topo05')
    index_x2g_Fgss_qice05 = shr_string_listGetIndexF(seq_flds_x2g_fields,'Fgss_qice05')
    index_x2g_Ss_tsrf06 = shr_string_listGetIndexF(seq_flds_x2g_fields,'Ss_tsrf06')
    index_x2g_Ss_topo06 = shr_string_listGetIndexF(seq_flds_x2g_fields,'Ss_topo06')
    index_x2g_Fgss_qice06 = shr_string_listGetIndexF(seq_flds_x2g_fields,'Fgss_qice06')
    index_x2g_Ss_tsrf07 = shr_string_listGetIndexF(seq_flds_x2g_fields,'Ss_tsrf07')
    index_x2g_Ss_topo07 = shr_string_listGetIndexF(seq_flds_x2g_fields,'Ss_topo07')
    index_x2g_Fgss_qice07 = shr_string_listGetIndexF(seq_flds_x2g_fields,'Fgss_qice07')
    index_x2g_Ss_tsrf08 = shr_string_listGetIndexF(seq_flds_x2g_fields,'Ss_tsrf08')
    index_x2g_Ss_topo08 = shr_string_listGetIndexF(seq_flds_x2g_fields,'Ss_topo08')
    index_x2g_Fgss_qice08 = shr_string_listGetIndexF(seq_flds_x2g_fields,'Fgss_qice08')
    index_x2g_Ss_tsrf09 = shr_string_listGetIndexF(seq_flds_x2g_fields,'Ss_tsrf09')
    index_x2g_Ss_topo09 = shr_string_listGetIndexF(seq_flds_x2g_fields,'Ss_topo09')
    index_x2g_Fgss_qice09 = shr_string_listGetIndexF(seq_flds_x2g_fields,'Fgss_qice09')
    index_x2g_Ss_tsrf10 = shr_string_listGetIndexF(seq_flds_x2g_fields,'Ss_tsrf10')
    index_x2g_Ss_topo10 = shr_string_listGetIndexF(seq_flds_x2g_fields,'Ss_topo10')
    index_x2g_Fgss_qice10 = shr_string_listGetIndexF(seq_flds_x2g_fields,'Fgss_qice10')

    nflds_x2g = shr_string_listGetNum(seq_flds_x2g_fields)

    !-------------------------------------------------------------
    ! drv -> sno
    !-------------------------------------------------------------

    index_x2s_Sg_frac01 = shr_string_listGetIndexF(seq_flds_x2s_fields,'Sg_frac01')
    index_x2s_Sg_topo01 = shr_string_listGetIndexF(seq_flds_x2s_fields,'Sg_topo01')
    index_x2s_Fsgg_rofi01 = shr_string_listGetIndexF(seq_flds_x2s_fields,'Fsgg_rofi01')
    index_x2s_Fsgg_rofl01 = shr_string_listGetIndexF(seq_flds_x2s_fields,'Fsgg_rofl01')
    index_x2s_Fsgg_hflx01 = shr_string_listGetIndexF(seq_flds_x2s_fields,'Fsgg_hflx01')
    index_x2s_Sg_frac02 = shr_string_listGetIndexF(seq_flds_x2s_fields,'Sg_frac02')
    index_x2s_Sg_topo02 = shr_string_listGetIndexF(seq_flds_x2s_fields,'Sg_topo02')
    index_x2s_Fsgg_rofi02 = shr_string_listGetIndexF(seq_flds_x2s_fields,'Fsgg_rofi02')
    index_x2s_Fsgg_rofl02 = shr_string_listGetIndexF(seq_flds_x2s_fields,'Fsgg_rofl02')
    index_x2s_Fsgg_hflx02 = shr_string_listGetIndexF(seq_flds_x2s_fields,'Fsgg_hflx02')
    index_x2s_Sg_frac03 = shr_string_listGetIndexF(seq_flds_x2s_fields,'Sg_frac03')
    index_x2s_Sg_topo03 = shr_string_listGetIndexF(seq_flds_x2s_fields,'Sg_topo03')
    index_x2s_Fsgg_rofi03 = shr_string_listGetIndexF(seq_flds_x2s_fields,'Fsgg_rofi03')
    index_x2s_Fsgg_rofl03 = shr_string_listGetIndexF(seq_flds_x2s_fields,'Fsgg_rofl03')
    index_x2s_Fsgg_hflx03 = shr_string_listGetIndexF(seq_flds_x2s_fields,'Fsgg_hflx03')
    index_x2s_Sg_frac04 = shr_string_listGetIndexF(seq_flds_x2s_fields,'Sg_frac04')
    index_x2s_Sg_topo04 = shr_string_listGetIndexF(seq_flds_x2s_fields,'Sg_topo04')
    index_x2s_Fsgg_rofi04 = shr_string_listGetIndexF(seq_flds_x2s_fields,'Fsgg_rofi04')
    index_x2s_Fsgg_rofl04 = shr_string_listGetIndexF(seq_flds_x2s_fields,'Fsgg_rofl04')
    index_x2s_Fsgg_hflx04 = shr_string_listGetIndexF(seq_flds_x2s_fields,'Fsgg_hflx04')
    index_x2s_Sg_frac05 = shr_string_listGetIndexF(seq_flds_x2s_fields,'Sg_frac05')
    index_x2s_Sg_topo05 = shr_string_listGetIndexF(seq_flds_x2s_fields,'Sg_topo05')
    index_x2s_Fsgg_rofi05 = shr_string_listGetIndexF(seq_flds_x2s_fields,'Fsgg_rofi05')
    index_x2s_Fsgg_rofl05 = shr_string_listGetIndexF(seq_flds_x2s_fields,'Fsgg_rofl05')
    index_x2s_Fsgg_hflx05 = shr_string_listGetIndexF(seq_flds_x2s_fields,'Fsgg_hflx05')
    index_x2s_Sg_frac06 = shr_string_listGetIndexF(seq_flds_x2s_fields,'Sg_frac06')
    index_x2s_Sg_topo06 = shr_string_listGetIndexF(seq_flds_x2s_fields,'Sg_topo06')
    index_x2s_Fsgg_rofi06 = shr_string_listGetIndexF(seq_flds_x2s_fields,'Fsgg_rofi06')
    index_x2s_Fsgg_rofl06 = shr_string_listGetIndexF(seq_flds_x2s_fields,'Fsgg_rofl06')
    index_x2s_Fsgg_hflx06 = shr_string_listGetIndexF(seq_flds_x2s_fields,'Fsgg_hflx06')
    index_x2s_Sg_frac07 = shr_string_listGetIndexF(seq_flds_x2s_fields,'Sg_frac07')
    index_x2s_Sg_topo07 = shr_string_listGetIndexF(seq_flds_x2s_fields,'Sg_topo07')
    index_x2s_Fsgg_rofi07 = shr_string_listGetIndexF(seq_flds_x2s_fields,'Fsgg_rofi07')
    index_x2s_Fsgg_rofl07 = shr_string_listGetIndexF(seq_flds_x2s_fields,'Fsgg_rofl07')
    index_x2s_Fsgg_hflx07 = shr_string_listGetIndexF(seq_flds_x2s_fields,'Fsgg_hflx07')
    index_x2s_Sg_frac08 = shr_string_listGetIndexF(seq_flds_x2s_fields,'Sg_frac08')
    index_x2s_Sg_topo08 = shr_string_listGetIndexF(seq_flds_x2s_fields,'Sg_topo08')
    index_x2s_Fsgg_rofi08 = shr_string_listGetIndexF(seq_flds_x2s_fields,'Fsgg_rofi08')
    index_x2s_Fsgg_rofl08 = shr_string_listGetIndexF(seq_flds_x2s_fields,'Fsgg_rofl08')
    index_x2s_Fsgg_hflx08 = shr_string_listGetIndexF(seq_flds_x2s_fields,'Fsgg_hflx08')
    index_x2s_Sg_frac09 = shr_string_listGetIndexF(seq_flds_x2s_fields,'Sg_frac09')
    index_x2s_Sg_topo09 = shr_string_listGetIndexF(seq_flds_x2s_fields,'Sg_topo09')
    index_x2s_Fsgg_rofi09 = shr_string_listGetIndexF(seq_flds_x2s_fields,'Fsgg_rofi09')
    index_x2s_Fsgg_rofl09 = shr_string_listGetIndexF(seq_flds_x2s_fields,'Fsgg_rofl09')
    index_x2s_Fsgg_hflx09 = shr_string_listGetIndexF(seq_flds_x2s_fields,'Fsgg_hflx09')
    index_x2s_Sg_frac10 = shr_string_listGetIndexF(seq_flds_x2s_fields,'Sg_frac10')
    index_x2s_Sg_topo10 = shr_string_listGetIndexF(seq_flds_x2s_fields,'Sg_topo10')
    index_x2s_Fsgg_rofi10 = shr_string_listGetIndexF(seq_flds_x2s_fields,'Fsgg_rofi10')
    index_x2s_Fsgg_rofl10 = shr_string_listGetIndexF(seq_flds_x2s_fields,'Fsgg_rofl10')
    index_x2s_Fsgg_hflx10 = shr_string_listGetIndexF(seq_flds_x2s_fields,'Fsgg_hflx10')

    nflds_x2s = shr_string_listGetNum(seq_flds_x2s_fields)

    !-------------------------------------------------------------
    ! sno -> drv
    !-------------------------------------------------------------

    index_s2x_Ss_tsrf01 = shr_string_listGetIndexF(seq_flds_s2x_fields,'Ss_tsrf01')
    index_s2x_Ss_topo01 = shr_string_listGetIndexF(seq_flds_s2x_fields,'Ss_topo01')
    index_s2x_Fgss_qice01 = shr_string_listGetIndexF(seq_flds_s2x_fields,'Fgss_qice01')
    index_s2x_Ss_tsrf02 = shr_string_listGetIndexF(seq_flds_s2x_fields,'Ss_tsrf02')
    index_s2x_Ss_topo02 = shr_string_listGetIndexF(seq_flds_s2x_fields,'Ss_topo02')
    index_s2x_Fgss_qice02 = shr_string_listGetIndexF(seq_flds_s2x_fields,'Fgss_qice02')
    index_s2x_Ss_tsrf03 = shr_string_listGetIndexF(seq_flds_s2x_fields,'Ss_tsrf03')
    index_s2x_Ss_topo03 = shr_string_listGetIndexF(seq_flds_s2x_fields,'Ss_topo03')
    index_s2x_Fgss_qice03 = shr_string_listGetIndexF(seq_flds_s2x_fields,'Fgss_qice03')
    index_s2x_Ss_tsrf04 = shr_string_listGetIndexF(seq_flds_s2x_fields,'Ss_tsrf04')
    index_s2x_Ss_topo04 = shr_string_listGetIndexF(seq_flds_s2x_fields,'Ss_topo04')
    index_s2x_Fgss_qice04 = shr_string_listGetIndexF(seq_flds_s2x_fields,'Fgss_qice04')
    index_s2x_Ss_tsrf05 = shr_string_listGetIndexF(seq_flds_s2x_fields,'Ss_tsrf05')
    index_s2x_Ss_topo05 = shr_string_listGetIndexF(seq_flds_s2x_fields,'Ss_topo05')
    index_s2x_Fgss_qice05 = shr_string_listGetIndexF(seq_flds_s2x_fields,'Fgss_qice05')
    index_s2x_Ss_tsrf06 = shr_string_listGetIndexF(seq_flds_s2x_fields,'Ss_tsrf06')
    index_s2x_Ss_topo06 = shr_string_listGetIndexF(seq_flds_s2x_fields,'Ss_topo06')
    index_s2x_Fgss_qice06 = shr_string_listGetIndexF(seq_flds_s2x_fields,'Fgss_qice06')
    index_s2x_Ss_tsrf07 = shr_string_listGetIndexF(seq_flds_s2x_fields,'Ss_tsrf07')
    index_s2x_Ss_topo07 = shr_string_listGetIndexF(seq_flds_s2x_fields,'Ss_topo07')
    index_s2x_Fgss_qice07 = shr_string_listGetIndexF(seq_flds_s2x_fields,'Fgss_qice07')
    index_s2x_Ss_tsrf08 = shr_string_listGetIndexF(seq_flds_s2x_fields,'Ss_tsrf08')
    index_s2x_Ss_topo08 = shr_string_listGetIndexF(seq_flds_s2x_fields,'Ss_topo08')
    index_s2x_Fgss_qice08 = shr_string_listGetIndexF(seq_flds_s2x_fields,'Fgss_qice08')
    index_s2x_Ss_tsrf09 = shr_string_listGetIndexF(seq_flds_s2x_fields,'Ss_tsrf09')
    index_s2x_Ss_topo09 = shr_string_listGetIndexF(seq_flds_s2x_fields,'Ss_topo09')
    index_s2x_Fgss_qice09 = shr_string_listGetIndexF(seq_flds_s2x_fields,'Fgss_qice09')
    index_s2x_Ss_tsrf10 = shr_string_listGetIndexF(seq_flds_s2x_fields,'Ss_tsrf10')
    index_s2x_Ss_topo10 = shr_string_listGetIndexF(seq_flds_s2x_fields,'Ss_topo10')
    index_s2x_Fgss_qice10 = shr_string_listGetIndexF(seq_flds_s2x_fields,'Fgss_qice10')

    nflds_s2x = shr_string_listGetNum(seq_flds_s2x_fields)


end subroutine seq_flds_indices_set

!===============================================================================

end module seq_flds_indices
