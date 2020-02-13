module microbeMod
#ifdef MICROBE
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: microbeMod
!
! !DESCRIPTION:
! Module holding routines to calculate methane fluxes
! The driver averages up to gridcell, weighting by finundated, and checks for balance errors.
! Sources, sinks, "competition" for CH4 & O2, & transport are resolved in ch4_tran.
!
! !USES:
	use shr_kind_mod, only : r8 => shr_kind_r8
	use clm_varpar, only : nlevsoi, ngases, nlevsno
	use clm_varcon, only : denh2o, denice, tfrz, grav, spval, rgas
	use clm_varcon, only : catomw, s_con, d_con_w, d_con_g, c_h_inv, kh_theta, kh_tbase, Fick_D_w, Henry_C_w, Henry_kHpc_w
	use clm_atmlnd, only : clm_a2l, clm_l2a
	use clm_time_manager, only : get_step_size, get_nstep
	use clm_varctl, only : iulog
	use abortutils, only : endrun
	USE microbevarcon
	use clm_varcon  , only : secspday

implicit none
	save
	private
	real(r8) :: f_sat = 0.99_r8 ! volumetric soil water defining top of water table or where production is allowed

! Non-tunable constants
	real(r8) :: rgasm  ! J/mol.K; rgas / 1000; will be set below
	real(r8), parameter :: rgasLatm = 0.0821_r8 ! L.atm/mol.K
! !PUBLIC MEMBER FUNCTIONS:
	public :: microbech4
	public :: microben2o
	public :: microbeCN
	public :: update_finundated
!	public :: seasonality
	private :: get_waterhead
#if (defined HUM_HOL)
	private :: lateral_bgc
#endif
!
! !REVISION HISTORY:
! 2013 Created by Xiaofeng Xu

!EOP
!-----------------------------------------------------------------------
contains
subroutine microbeCN (lbg, ubg, lbl, ubl, lbc, ubc, lbp, ubp, num_soilc, filter_soilc, &
                num_soilp, filter_soilp)
	!
! !USES:
	use clmtype
	use subgridAveMod , only : p2c, c2g
	use clm_varpar, only : nlevgrnd, nlevdecomp
	use pftvarcon,  only : noveg
	use clm_varcon, only : secspday
	use microbevarcon 
!	use nanmod   ,  only : nan, bigint
	use clmtype
!
! !ARGUMENTS:
implicit none
	integer, intent(in) :: lbg, ubg           			! grid-index bounds
	integer, intent(in) :: lbl, ubl           			! landunit-level index bounds
	integer, intent(in) :: lbc, ubc           			! column-index bounds
	integer, intent(in) :: lbp, ubp           			! pft-level index bounds
	integer, intent(in) :: num_soilc          			! number of column soil points in column filter
	integer, intent(in) :: filter_soilc(ubc-lbc+1)    	! column filter for soil points
	integer, intent(in) :: num_soilp          			! number of soil points in pft filter
	integer, intent(in) :: filter_soilp(ubp-lbp+1)      	! pft filter for soil points

end subroutine microbeCN

subroutine microben2o (lbg, ubg, lbl, ubl, lbc, ubc, lbp, ubp, num_soilc, filter_soilc, &
                num_soilp, filter_soilp)
	!
! !USES:
	use clmtype
	use subgridAveMod , only : p2c, c2g
	use clm_varpar, only : nlevgrnd, nlevdecomp
	use pftvarcon,  only : noveg
	use clm_varcon, only : secspday
	use microbevarcon 
!	use nanmod   ,  only : nan, bigint
	use clmtype
!
! !ARGUMENTS:
implicit none
	integer, intent(in) :: lbg, ubg           			! grid-index bounds
	integer, intent(in) :: lbl, ubl           			! landunit-level index bounds
	integer, intent(in) :: lbc, ubc           			! column-index bounds
	integer, intent(in) :: lbp, ubp           			! pft-level index bounds
	integer, intent(in) :: num_soilc          			! number of column soil points in column filter
	integer, intent(in) :: filter_soilc(ubc-lbc+1)    	! column filter for soil points
	integer, intent(in) :: num_soilp          			! number of soil points in pft filter
	integer, intent(in) :: filter_soilp(ubp-lbp+1)      	! pft filter for soil points
end subroutine microben2o
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ch4
!
! !INTERFACE:
subroutine microbech4 (lbg, ubg, lbl, ubl, lbc, ubc, lbp, ubp, num_soilc, filter_soilc, &
                num_soilp, filter_soilp)
! !DESCRIPTION:
! Driver for the methane emissions model
!
! !USES:
	use shr_kind_mod , only : r8 => shr_kind_r8
	use clmtype
	use subgridAveMod , only : p2c, c2g
	use clm_varpar, only : nlevgrnd, nlevdecomp, i_bacteria, i_fungi, i_dom, cn_dom
	use pftvarcon,  only : noveg
	use clm_varcon, only : secspday, istwet, istsoil, istdlak, spval, istcrop
	use microbevarcon 
	use clm_time_manager, only : get_step_size, get_nstep
!	use nanmod,  only : nan, bigint
	use clmtype
	use filterMod,  only : filter
	use shr_const_mod, only: SHR_CONST_TKFRZ, SHR_CONST_RGAS
!
! !ARGUMENTS:
implicit none
	integer, intent(in) :: lbg, ubg           			! grid-index bounds
	integer, intent(in) :: lbl, ubl           			! landunit-level index bounds
	integer, intent(in) :: lbc, ubc           			! column-index bounds
	integer, intent(in) :: lbp, ubp           			! pft-level index bounds
	integer, intent(in) :: num_soilc          			! number of column soil points in column filter
	integer, intent(in) :: filter_soilc(ubc-lbc+1)    	! column filter for soil points
	integer, intent(in) :: num_soilp          			! number of soil points in pft filter
	integer, intent(in) :: filter_soilp(ubp-lbp+1) 		! pft filter for soil points
! !CALLED FROM:
! driver.F90
!
! !REVISION HISTORY:
! !LOCAL VARIABLES:
! local pointers to implicit in variables
	real(r8), pointer :: gmicbios(:)				! grid-level biomass of microbes molC/m2
	real(r8), pointer :: gdocs(:)					! grid-level doc concentration molC/m2
	real(r8), pointer :: gaces(:)					! grid-level acetate concentration molC/m2
	real(r8), pointer :: gacebios(:)				! grid-level biomass of methanogen based on acetate molC/m2
	real(r8), pointer :: gco2bios(:)				! grid-level biomass of methanogen based on CO2 and H2 molC/m2
	real(r8), pointer :: gaerch4bios(:)				! grid-level biomass of aerobix methanotrophy molC/m2
	real(r8), pointer :: ganaerch4bios(:)			! grid-level biomass of anaerobic methanotrophy molC/m2
   
	real(r8), pointer :: cmicbiocs(:,:)				! column-level biomass of all microbes molC/m3
	real(r8), pointer :: cdocs_pre(:,:)				! column-level concentration of DOC molC/m3 
	real(r8), pointer :: cdocs(:,:)					! column-level concentration of DOC molC/m3 
	real(r8), pointer :: cdocs_unsat(:,:)			! column-level concentration of DOC molC/m3 in unsaturated fraction
	real(r8), pointer :: cdocs_sat(:,:)				! column-level concentration of DOC molC/m3 in saturated fraction
	real(r8), pointer :: cmicbions(:,:)				! column-level biomass of all microbes molN/m3
	real(r8), pointer :: cdons(:,:)					! column-level concentration of DON molN/m3 
	real(r8), pointer :: cdons_unsat(:,:)			! column-level concentration of DON molN/m3 in unsaturated fraction
	real(r8), pointer :: cdons_sat(:,:)				! column-level concentration of DON molN/m3 in saturated fraction
	real(r8), pointer :: cdons_min(:,:)				! minteralization of DON
	real(r8), pointer :: caces(:,:)					! column-level concentration of acetate molC/m3
	real(r8), pointer :: cacebios(:,:)				! column-level biomass of methanogen based on acetate molC/m3
	real(r8), pointer :: cco2bios(:,:)				! column-level biomass of methanogen based on CO2/H2 molC/m3
	real(r8), pointer :: caerch4bios(:,:)			! column-level biomass of aerobix methanotrophy molC/m3
	real(r8), pointer :: canaerch4bios(:,:)			! column-level biomass of anaerobic methanotrophy molC/m3
    
	real(r8), pointer :: cacebios_unsat(:,:)			! column-level biomass of methanogen from acetate in unsaturated fraction molC/m3
	real(r8), pointer :: cacebios_sat(:,:)			! column-level biomass of methanogen from acetate in saturated fraction molC/m3
	real(r8), pointer :: cco2bios_unsat(:,:)			! column-level biomass of methanogen from CO2 and H2 in unsaturated fraction molC/m3
	real(r8), pointer :: cco2bios_sat(:,:)			! column-level biomass of methanogen from CO2 and H2 in saturated fraction molC/m3
	real(r8), pointer :: caerch4bios_unsat(:,:)		! column-level biomass of aerobic methanotrophy in unsaturated fraction molC/m3
	real(r8), pointer :: caerch4bios_sat(:,:)			! column-level biomass of aerobic methanotrophy in saturated fraction molC/m3
	real(r8), pointer :: canaerch4bios_unsat(:,:)		! column-level biomass of anaerobic methanotrophy in unsaturated fraction molC/m3
	real(r8), pointer :: canaerch4bios_sat(:,:)		! column-level biomass of anaerobic methanotrophy in saturated fraction molC/m3
	real(r8), pointer :: caces_unsat(:,:)			! column-level acetate in unsaturated fraction       molC/m3
	real(r8), pointer :: caces_sat(:,:)				! column-level acetate in saturated fraction         molC/m3
	
	real(r8), pointer :: caces_prod(:,:)				! acetate acid production
	real(r8), pointer :: caces_unsat_prod(:,:)		! acetate acid produciton in unsaturation fraction of soil column
	real(r8), pointer :: caces_sat_prod(:,:)			! acetate acid produciton in saturation fraction of soil column

	real(r8), pointer :: caces_prod_h2(:,:)				! acetogenesis
	real(r8), pointer :: caces_unsat_prod_h2(:,:)		! acetogenesis in unsaturation fraction of soil column
	real(r8), pointer :: caces_sat_prod_h2(:,:)			! acetogenesis in saturation fraction of soil column
  
	real(r8), pointer :: ccon_ch4s(:,:)				! column-level concentration of CH4 mol C/m3
	real(r8), pointer :: ccon_o2s(:,:)				! column-level concentration of O2 mol O2/m3
	real(r8), pointer :: ccon_co2s(:,:)				! column-level concentration of CO2 mol C/m3
	real(r8), pointer :: ccon_h2s(:,:)				! column-level concentrtation of H2 mol H/m3
	real(r8), pointer :: ccon_ch4s_unsat(:,:)		! column-level concentration of CH4 in unsaturated fraction mol C/m3
	real(r8), pointer :: ccon_ch4s_sat(:,:)			! column-level concentration of CH4 in saturated fraction mol C/m3
	real(r8), pointer :: ccon_o2s_unsat(:,:)			! column-level concentration of O2 in unsaturated fraction mol O2/m3
	real(r8), pointer :: ccon_o2s_sat(:,:)			! column-level concentration of O2 in saturated fraction mol O2/m3
	real(r8), pointer :: ccon_co2s_unsat(:,:)			! column-level concentration of CO2 in unsaturated fraction mol C/m3
	real(r8), pointer :: ccon_co2s_sat(:,:)			! column-level concentration of CO2 in saturated fraction mol C/m3
	real(r8), pointer :: ccon_h2s_unsat(:,:)			! column-level concentration of H2 in unsaturated fraction mol H2/m3
	real(r8), pointer :: ccon_h2s_sat(:,:)			! column-level concentration of H2 in saturated fraction  mol H2/m3
	
	integer, pointer :: ltype(:)					! type of land unit
	integer, pointer :: ptype(:)					! type of plant functional type
	integer, pointer :: clandunit(:)				! index of land nnit in column
	real(r8), pointer :: soilpH_unsat(:,:)			! soil pH for unsaturated fraction soil column
	real(r8), pointer :: soilpH_sat(:,:)				! soil pH for saturated fraction soil column
	real(r8), pointer :: soiltemp(:,:)				! soil temperature
	real(r8), pointer :: sand(:)					! sand
	real(r8), pointer :: silt(:)					! silt
!	real(r8), pointer :: psi(:)						! water potential
!	real(r8), pointer :: psisat(:)						! water potential at saturated
	real(r8), pointer :: vwc(:,:)					! volumetic water content
	real(r8), pointer :: vwcsat(:,:)				! volumetic water content at saturation

	real(r8), pointer :: fsat_pre(:)    				! finundated from previous timestep
	real(r8), pointer :: c_atm(:,:)     			! CH4, O2, CO2 atmospheric conc  (mol/m3)
	real(r8), pointer :: flux_ch4(:)    			! gridcell CH4 flux to atm. (kg C/m**2/s)

	real(r8), pointer :: nee_ch4(:)				! gridcell average net methane correction to CO2 flux (g C/m^2/s)
	real(r8), pointer :: ch4_dfsat_flux(:)			! CH4 flux to atm due to decreasing finundated (kg C/m^2/s) [+]
	real(r8), pointer :: ch4prodg(:)				! gridcell average CH4 production (g C/m^2/s)
	real(r8), pointer :: zwt_ch4_unsat(:) 			! depth of water table for unsaturated fraction (m)
	real(r8), pointer :: rootfr_vr(:,:)				! fraction of roots in each soil layer  (nlevgrnd)
	real(r8), pointer :: finundated(:)				! fractional inundated area in soil column (excluding dedicated wetland columns)
!	real(r8), pointer :: fsat(:)					! fractional inundated area in soil column (excluding dedicated wetland columns)
   
!clm-microbe
	real(r8), pointer :: origionalsoilph(:)			! original soil pH value
	integer , pointer :: cgridcell(:)   				! gridcell of corresponding column
	real(r8), pointer :: forc_t(:)      				! atmospheric temperature (Kelvin)
	real(r8), pointer :: forc_pbot(:)   			! atmospheric pressure (Pa)
	real(r8), pointer :: forc_pco2(:)   			! CO2 partial pressure (Pa)
	real(r8), pointer :: forc_pch4(:)   			! CH4 partial pressure (Pa)
	real(r8), pointer :: forc_po2(:)    			! O2 partial pressure (Pa)
	real(r8), pointer :: forc_ph2(:)    			! H2 partial pressure (Pa)
	real(r8), pointer :: dz(:,:)        				! layer thickness (m)  (-nlevsno+1:nlevsoi)
	real(r8), pointer :: zi(:,:)        				! interface level below a "z" level (m)
	real(r8), pointer :: z(:,:)         				! layer depth (m) (-nlevsno+1:nlevsoi)
	integer , pointer :: pcolumn(:)     			! index into column level quantities
	real(r8), pointer :: wtcol(:)     				! weight of each plant funcitonal type in one column
!   real(r8), pointer :: zwt0(:)        				! decay factor for finundated (m)
!   real(r8), pointer :: f0(:)          				! maximum gridcell fractional inundated area
!   real(r8), pointer :: p3(:)          				! coefficient for qflx_surf_lag for finunated (s/mm)
!   real(r8), pointer :: zwt(:)         				! water table depth (m) (from SoilHydrology)
   
!	real(r8), pointer :: tempavg_bgnpp(:)   				! temporary average below-ground NPP (gC/m2/s)
!	real(r8), pointer :: annavg_bgnpp(:)    				! annual average below-ground NPP (gC/m2/s)
	integer :: jwaterhead_unsat(lbc:ubc)				! layer of the water head in unsaturated fraction
	real(r8), pointer :: waterhead_unsat(:)        		! water table depth (m) (from SoilHydrology)
   
!	real(r8), pointer :: bgnpp_timestep(:) 
!	real(r8), pointer :: bgnpp_avg(:) 
   
   	real(r8), pointer :: ch4_prod_ace_depth_unsat(:,:)
	real(r8), pointer :: ch4_prod_co2_depth_unsat(:,:)
	real(r8), pointer :: ch4_oxid_o2_depth_unsat(:,:)
	real(r8), pointer :: ch4_oxid_aom_depth_unsat(:,:)
	real(r8), pointer :: ch4_aere_depth_unsat(:,:)
	real(r8), pointer :: ch4_dif_depth_unsat(:,:)
	real(r8), pointer :: ch4_ebul_depth_unsat(:,:)
	real(r8), pointer :: co2_prod_ace_depth_unsat(:,:)
	real(r8), pointer :: co2_decomp_depth_unsat(:,:)
	real(r8), pointer :: co2_cons_depth_unsat(:,:)
	real(r8), pointer :: co2_ebul_depth_unsat(:,:)
	real(r8), pointer :: co2_aere_depth_unsat(:,:)
	real(r8), pointer :: co2_dif_depth_unsat(:,:)
	real(r8), pointer :: o2_cons_depth_unsat(:,:)
	real(r8), pointer :: o2_aere_depth_unsat(:,:)
	real(r8), pointer :: o2_aere_oxid_depth_unsat(:,:)
	real(r8), pointer :: o2_decomp_depth_unsat(:,:)
	real(r8), pointer :: o2_dif_depth_unsat(:,:)
	real(r8), pointer :: h2_prod_depth_unsat(:,:)
	real(r8), pointer :: h2_cons_depth_unsat(:,:)
	real(r8), pointer :: h2_aere_depth_unsat(:,:)
	real(r8), pointer :: h2_diff_depth_unsat(:,:)
	real(r8), pointer :: h2_ebul_depth_unsat(:,:)
	real(r8), pointer :: ch4_surf_aere_unsat(:)
	real(r8), pointer :: ch4_surf_ebul_unsat(:)
	real(r8), pointer :: ch4_surf_dif_unsat(:)
	real(r8), pointer :: ch4_surf_netflux_unsat(:)
	real(r8), pointer :: co2_surf_aere_unsat(:)
	real(r8), pointer :: co2_surf_ebul_unsat(:)
	real(r8), pointer :: co2_surf_dif_unsat(:)
	real(r8), pointer :: co2_surf_netflux_unsat(:)
	real(r8), pointer :: o2_surf_aere_unsat(:)
	real(r8), pointer :: o2_surf_dif_unsat(:)
	real(r8), pointer :: o2_surf_netflux_unsat(:)
	real(r8), pointer :: h2_surf_aere_unsat(:)
	real(r8), pointer :: h2_surf_ebul_unsat(:)
	real(r8), pointer :: h2_surf_dif_unsat(:)
	real(r8), pointer :: h2_surf_netflux_unsat(:)

	real(r8), pointer :: ch4_prod_ace_depth_sat(:,:)
	real(r8), pointer :: ch4_prod_co2_depth_sat(:,:)
	real(r8), pointer :: ch4_oxid_o2_depth_sat(:,:)
	real(r8), pointer :: ch4_oxid_aom_depth_sat(:,:)
	real(r8), pointer :: ch4_aere_depth_sat(:,:)
	real(r8), pointer :: ch4_dif_depth_sat(:,:)
	real(r8), pointer :: ch4_ebul_depth_sat(:,:)
	real(r8), pointer :: co2_prod_ace_depth_sat(:,:)
	real(r8), pointer :: co2_decomp_depth_sat(:,:)
	real(r8), pointer :: co2_cons_depth_sat(:,:)
	real(r8), pointer :: co2_ebul_depth_sat(:,:)
	real(r8), pointer :: co2_aere_depth_sat(:,:)
	real(r8), pointer :: co2_dif_depth_sat(:,:)
	real(r8), pointer :: o2_cons_depth_sat(:,:)
	real(r8), pointer :: o2_aere_depth_sat(:,:)
	real(r8), pointer :: o2_aere_oxid_depth_sat(:,:)
	real(r8), pointer :: o2_decomp_depth_sat(:,:)
	real(r8), pointer :: o2_dif_depth_sat(:,:)
	real(r8), pointer :: h2_prod_depth_sat(:,:)
	real(r8), pointer :: h2_cons_depth_sat(:,:)
	real(r8), pointer :: h2_aere_depth_sat(:,:)
	real(r8), pointer :: h2_diff_depth_sat(:,:)
	real(r8), pointer :: h2_ebul_depth_sat(:,:)
	real(r8), pointer :: ch4_surf_aere_sat(:)
	real(r8), pointer :: ch4_surf_ebul_sat(:)
	real(r8), pointer :: ch4_surf_dif_sat(:)
	real(r8), pointer :: ch4_surf_netflux_sat(:)
	real(r8), pointer :: co2_surf_aere_sat(:)
	real(r8), pointer :: co2_surf_ebul_sat(:)
	real(r8), pointer :: co2_surf_dif_sat(:)
	real(r8), pointer :: co2_surf_netflux_sat(:)
	real(r8), pointer :: o2_surf_aere_sat(:)
	real(r8), pointer :: o2_surf_dif_sat(:)
	real(r8), pointer :: o2_surf_netflux_sat(:)
	real(r8), pointer :: h2_surf_aere_sat(:)
	real(r8), pointer :: h2_surf_ebul_sat(:)
	real(r8), pointer :: h2_surf_dif_sat(:)
	real(r8), pointer :: h2_surf_netflux_sat(:)

	real(r8), pointer :: ch4_prod_ace_depth(:,:)
	real(r8), pointer :: ch4_prod_co2_depth(:,:)
	real(r8), pointer :: ch4_oxid_o2_depth(:,:)
	real(r8), pointer :: ch4_oxid_aom_depth(:,:)
	real(r8), pointer :: ch4_aere_depth(:,:)
	real(r8), pointer :: ch4_dif_depth(:,:)
	real(r8), pointer :: ch4_ebul_depth(:,:)
	real(r8), pointer :: co2_prod_ace_depth(:,:)
	real(r8), pointer :: co2_decomp_depth(:,:)
	real(r8), pointer :: co2_cons_depth(:,:)
	real(r8), pointer :: co2_ebul_depth(:,:)
	real(r8), pointer :: co2_aere_depth(:,:)
	real(r8), pointer :: co2_dif_depth(:,:)
	real(r8), pointer :: o2_cons_depth(:,:)
	real(r8), pointer :: o2_aere_depth(:,:)
	real(r8), pointer :: o2_aere_oxid_depth(:,:)
	real(r8), pointer :: o2_decomp_depth(:,:)
	real(r8), pointer :: o2_dif_depth(:,:)
	real(r8), pointer :: h2_prod_depth(:,:)
	real(r8), pointer :: h2_cons_depth(:,:)
	real(r8), pointer :: h2_aere_depth(:,:)
	real(r8), pointer :: h2_diff_depth(:,:)
	real(r8), pointer :: h2_ebul_depth(:,:)
	real(r8), pointer :: ch4_surf_aere(:)
	real(r8), pointer :: ch4_surf_ebul(:)
	real(r8), pointer :: ch4_surf_dif(:)
	real(r8), pointer :: ch4_surf_netflux(:)
	real(r8), pointer :: co2_surf_aere(:)
	real(r8), pointer :: co2_surf_ebul(:)
	real(r8), pointer :: co2_surf_dif(:)
	real(r8), pointer :: co2_surf_netflux(:)
	real(r8), pointer :: o2_surf_aere(:)
	real(r8), pointer :: o2_surf_dif(:)
	real(r8), pointer :: o2_surf_netflux(:)
	real(r8), pointer :: h2_surf_aere(:)
	real(r8), pointer :: h2_surf_ebul(:)
	real(r8), pointer :: h2_surf_dif(:)
	real(r8), pointer :: h2_surf_netflux(:)
	
	real(r8), pointer :: hr_vr(:,:)
	real(r8), pointer :: roothr(:)
	real(r8), pointer :: froot_mr(:)
	real(r8), pointer :: froot_r(:,:)
	
	real(r8), pointer :: sminn_vr(:,:)

	real(r8), pointer :: sucsat(:,:)        		! minimum soil suction (mm)
	real(r8), pointer :: soilpsi(:,:)        		! soil water potential in each soil layer (MPa)
   
	real(r8) :: caces_unsat_temp(lbc:ubc,1:nlevgrnd)	! temporary array
	real(r8) :: caces_sat_temp(lbc:ubc,1:nlevgrnd)	! temporary array
	real(r8) :: cdocs_unsat_temp(lbc:ubc,1:nlevgrnd)	! temporary array
	real(r8) :: cdocs_sat_temp(lbc:ubc,1:nlevgrnd)	! temporary array
	
	real(r8), pointer :: annsum_npp(:)   		! annual sum NPP (gC/m2/yr)
	real(r8), pointer :: annavg_agnpp(:) 		! (gC/m2/s) annual average aboveground NPP
	real(r8), pointer :: annavg_bgnpp(:) 		! (gC/m2/s) annual average belowground NPP
	real(r8), pointer :: col_npp(:) 			! (gC/m2/s) colunm net primary production
	real(r8), pointer :: col_rr(:) 				! (gC/m2/s) colunm root respiration
	real(r8), pointer :: cannsum_npp(:) 			!  annual nitrogen depsotion rate
	
	real(r8), pointer :: decomp_cpools_vr(:,:,:)    ! (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) c pools
	real(r8), pointer :: decomp_npools_vr(:,:,:)    ! (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools

#if (defined HUM_HOL)
!	real(r8), pointer :: h2osoi_vol(:,:)   ! volumetric soil water(0<=h2osoi_vol<=watsat) [m3/m3]
	real(r8), pointer :: qflx_lat_aqu(:)   		! lateral flow in aquifer (mm/s)
	real(r8), pointer :: wa(:)   			! lateral flow in aquifer (mm/s)
	real(r8), pointer :: zwt(:)            		! water table depth (m)
	real(r8) :: jwt(1:2)
#endif

	real(r8) :: nppratio(lbc:ubc)
	
	real(r8) :: roothr_vr(lbc:ubc,1:nlevdecomp)
	real(r8) :: rootfraction(lbc:ubc,1:nlevdecomp)
	
	real(r8) :: AceProd = 0.   				! acetate production
	real(r8) :: ACH2Prod = 0.				! h2 production from dissolved organic carbon mineralization
	real(r8) :: ACCO2Prod = 0.				! co2 production from dissolved organic carbon mineralization
	real(r8) :: H2CH4Prod = 0.				! ch4 production from h2
	real(r8) :: H2AceProd = 0.				! ace produciton from h2
	real(r8) :: H2Cons = 0.				! h2 consumption
	real(r8) :: CH4PlantFlux = 0.			! plant transport of ch4
	real(r8) :: CH4Ebull = 0.				! ebbulltion of ch4
	real(r8) :: PlantO2Cons = 0.				! o2 consumption when transport from atmosphere to the root
	real(r8) :: AerO2Cons = 0.				! Aerobic consumption of O2
	real(r8) :: CH4O2Cons = 0.				! CH4 oxidation with O2
	real(r8) :: dtempratio = 0.				! temporary variable
	real(r8) :: dtemph2 = 0.				! temporary variable for h2
	real(r8) :: H2PlantFlux = 0.				! h2 transport through plant
	real(r8) :: CH4Prod = 0.				! CH4 proudction
	real(r8) :: CH4Oxid = 0.				! CH4 oxidation
	real(r8) :: O2PlantFlux = 0.				! O2 trasnport through plant
	real(r8) :: CO2Prod = 0.				! CO2 proudction in this module (BGC transition in the CN-Microbe excluded)
	real(r8) :: CO2PlantFlux = 0.			! CO2 transport through plant
	real(r8) :: AOMCH4Oxid = 0.				! Anaerobic oxidation of methane 
	real(r8) :: AceCons = 0.				! Acetate acid consumption
	real(r8) :: H2CO2Cons = 0.				! H2 production CO2
	real(r8) :: AceMethanogenGrowth = 0.		! growing rate of aceclastic methanogen
	real(r8) :: AceMethanogenDying = 0.		! death rate of aceclastic methanogen
	real(r8) :: H2MethanogenGrowth = 0.		! growing rate of hydrogen-based methanogen
	real(r8) :: H2MethanogenDying = 0.		! death rate of hydrogen-based methanogen
	real(r8) :: MethanotrophGrowth = 0.		! growing rate of methanotroph
	real(r8) :: MethanotrophDying = 0.		! death rate of methanotroph
	real(r8) :: AOMMethanotrophGrowth = 0.	! growing rate of the anerobic oxidation of methane
	real(r8) :: AOMMethanotrophDying = 0.		! death rate of anaerobic oxidation of methane 
	
	real(r8) :: tem1, tem2, tem3, tem4		! four temporary variables
	real(r8) :: anpp!, nppratio				! above-ground npp, ratio of daily npp to annual total npp
        real(r8) :: OxidAce2CO2                           ! oxidaiton of acetate to co2 (nromally sufficient oxygen), distinguish from aceteclastic methanogenesis

! !OTHER LOCAL VARIABLES:
	integer :: fc,c,g,j,fp,p,l               ! indices
!	real(r8) :: rootfraction(lbp:ubp, 1:nlevgrnd) 
	real(r8) :: pHeffect = 1.
	real(r8) :: ACConcentration = 0.
	real(r8) :: HCH4Prod = 0.
	real(r8):: sumdoc = 0.
	real(r8):: sumace = 0.
	integer :: IsH2Production = 1
	real(r8):: dt        ! time step (seconds)
!EOP

	real(r8):: minpsi, maxpsi                		! limits for soil water scalar for decomp
	real(r8):: psi                           		! temporary soilpsi for water scalar
	real(r8):: micfinundated = 0.99                ! temporary soilpsi for water scalar
   
#if (defined HUM_HOL)
!	real(r8) :: h2osoi_vol
	real(r8) :: hum_frac
	real(r8) :: hol_frac
	real(r8) :: lat_flux_factor
	real(r8) :: lat_flux_factor1
	real(r8) :: lat_flux_factor2
		
	real(r8) :: lxdomunsat
	real(r8) :: lxdomsat
	real(r8) :: lxaceunsat
	real(r8) :: lxacesat
#endif
	real(r8) :: som_diffus = 1e-4_r8 / (secspday * 365._r8)  ! [m^2/sec] = 1 cm^2 / yr
!	real(r8) :: dom_diffus = 1e8_r8							! times to som diffus 


!-----------------------------------------------------------------------
   ! Gridcell level pointers
	forc_t    					=> clm_a2l%forc_t
	forc_pbot 					=> clm_a2l%forc_pbot
	c_atm     					=> gmic%c_atm
	forc_po2  					=> clm_a2l%forc_po2
	forc_pco2 					=> clm_a2l%forc_pco2
	forc_pch4 					=> clm_a2l%forc_pch4
	forc_ph2 					=> clm_a2l%forc_ph2
   
   !flux_ch4  => clm_l2a%flux_ch4
   !ch4co2f   => gch4%ch4co2f
   !ch4prodg  => gch4%ch4prodg
   !atdeg    => latdeg

	gmicbios 					=> gmic%gmicbios
	gdocs 					=> gmic%gdocs
	gaces 					=> gmic%gaces
	gacebios 					=> gmic%gacebios
	gco2bios 					=> gmic%gco2bios
	gaerch4bios 				=> gmic%gaerch4bios
	ganaerch4bios 				=> gmic%ganaerch4bios   

   ! Column level pointers
	cmicbiocs 					=> cmic%cmicbiocs
	cdocs_pre 					=> cmic%cdocs_pre
	cdocs 					=> cmic%cdocs
	cdocs_unsat 				=> cmic%cdocs_unsat
	cdocs_sat					=> cmic%cdocs_sat
	cmicbions 					=> cmic%cmicbions
	cdons 					=> cmic%cdons
	cdons_unsat 				=> cmic%cdons_unsat
	cdons_sat					=> cmic%cdons_sat
	cdons_min					=> cmic%cdons_min
	
	caces_prod 				=> cmic%caces_prod
	caces_unsat_prod 			=> cmic%caces_unsat_prod
	caces_sat_prod				=> cmic%caces_sat_prod

	caces_prod_h2 				=> cmic%caces_prod_h2
	caces_unsat_prod_h2 		=> cmic%caces_unsat_prod_h2
	caces_sat_prod_h2			=> cmic%caces_sat_prod_h2
		
	caces 					=> cmic%caces
	caces_unsat 				=> cmic%caces_unsat
	caces_sat 					=> cmic%caces_sat
	cacebios 					=> cmic%cacebios
	cacebios_unsat 				=> cmic%cacebios_unsat
	cacebios_sat 				=> cmic%cacebios_sat
	cco2bios 					=> cmic%cco2bios
	cco2bios_unsat 				=> cmic%cco2bios_unsat
	cco2bios_sat 				=> cmic%cco2bios_sat
	caerch4bios 				=> cmic%caerch4bios
	caerch4bios_unsat 			=> cmic%caerch4bios_unsat
	caerch4bios_sat 				=> cmic%caerch4bios_sat
	canaerch4bios 				=> cmic%canaerch4bios
	canaerch4bios_unsat 			=> cmic%canaerch4bios_unsat
	canaerch4bios_sat 			=> cmic%canaerch4bios_sat
     
	ccon_ch4s        				=> cmic%ccon_ch4s
	ccon_co2s         				=> cmic%ccon_co2s
	ccon_o2s          				=> cmic%ccon_o2s
	ccon_h2s          				=> cmic%ccon_h2s
	ccon_ch4s_unsat    			=> cmic%ccon_ch4s_unsat
	ccon_ch4s_sat        			=> cmic%ccon_ch4s_sat
	ccon_co2s_unsat    			=> cmic%ccon_co2s_unsat
	ccon_co2s_sat        			=> cmic%ccon_co2s_sat
	ccon_o2s_unsat      			=> cmic%ccon_o2s_unsat
	ccon_o2s_sat          			=> cmic%ccon_o2s_sat
	ccon_h2s_unsat      			=> cmic%ccon_h2s_unsat
	ccon_h2s_sat          			=> cmic%ccon_h2s_sat
 
! start of the transport code
	ch4_prod_ace_depth_unsat		=> cmic%ch4_prod_ace_depth_unsat
	ch4_prod_co2_depth_unsat 		=> cmic%ch4_prod_co2_depth_unsat
	ch4_oxid_o2_depth_unsat 		=> cmic%ch4_oxid_o2_depth_unsat
	ch4_oxid_aom_depth_unsat		=> cmic%ch4_oxid_aom_depth_unsat
	ch4_aere_depth_unsat 			=> cmic%ch4_aere_depth_unsat
	ch4_dif_depth_unsat 			=> cmic%ch4_dif_depth_unsat
	ch4_ebul_depth_unsat 			=> cmic%ch4_ebul_depth_unsat
	co2_prod_ace_depth_unsat		=> cmic%co2_prod_ace_depth_unsat
	co2_decomp_depth_unsat 		=> cmic%co2_decomp_depth_unsat
	co2_cons_depth_unsat 		=> cmic%co2_cons_depth_unsat
	co2_ebul_depth_unsat 			=> cmic%co2_ebul_depth_unsat
	co2_aere_depth_unsat 			=> cmic%co2_aere_depth_unsat
	co2_dif_depth_unsat 			=> cmic%co2_dif_depth_unsat
	o2_cons_depth_unsat 			=> cmic%o2_cons_depth_unsat
	o2_aere_depth_unsat 			=> cmic%o2_aere_depth_unsat
	o2_aere_oxid_depth_unsat 		=> cmic%o2_aere_oxid_depth_unsat
	o2_decomp_depth_unsat 		=> cmic%o2_decomp_depth_unsat
	o2_cons_depth_unsat 			=> cmic%o2_cons_depth_unsat
	o2_dif_depth_unsat 			=> cmic%o2_dif_depth_unsat
	h2_prod_depth_unsat 			=> cmic%h2_prod_depth_unsat
	h2_cons_depth_unsat 			=> cmic%h2_cons_depth_unsat
	h2_aere_depth_unsat 			=> cmic%h2_aere_depth_unsat
	h2_diff_depth_unsat 			=> cmic%h2_diff_depth_unsat
	h2_ebul_depth_unsat 			=> cmic%h2_ebul_depth_unsat
	ch4_surf_aere_unsat			=> cmic%ch4_surf_aere_unsat
	ch4_surf_ebul_unsat 			=> cmic%ch4_surf_ebul_unsat
	ch4_surf_dif_unsat 			=> cmic%ch4_surf_dif_unsat
	ch4_surf_netflux_unsat 		=> cmic%ch4_surf_netflux_unsat
	co2_surf_aere_unsat			=> cmic%co2_surf_aere_unsat
	co2_surf_ebul_unsat 			=> cmic%co2_surf_ebul_unsat
	co2_surf_dif_unsat 			=> cmic%co2_surf_dif_unsat
	co2_surf_netflux_unsat 		=> cmic%co2_surf_netflux_unsat
	o2_surf_aere_unsat			=> cmic%o2_surf_aere_unsat
	o2_surf_dif_unsat 			=> cmic%o2_surf_dif_unsat
	o2_surf_netflux_unsat 			=> cmic%o2_surf_netflux_unsat
	h2_surf_aere_unsat			=> cmic%h2_surf_aere_unsat
	h2_surf_ebul_unsat 			=> cmic%h2_surf_ebul_unsat
	h2_surf_dif_unsat 			=> cmic%h2_surf_dif_unsat
	h2_surf_netflux_unsat 			=> cmic%h2_surf_netflux_unsat

	ch4_prod_ace_depth_sat 		=> cmic%ch4_prod_ace_depth_sat
	ch4_prod_co2_depth_sat 		=> cmic%ch4_prod_co2_depth_sat
	ch4_oxid_o2_depth_sat 		=> cmic%ch4_oxid_o2_depth_sat
	ch4_oxid_aom_depth_sat 		=> cmic%ch4_oxid_aom_depth_sat
	ch4_aere_depth_sat 			=> cmic%ch4_aere_depth_sat
	ch4_dif_depth_sat 			=> cmic%ch4_dif_depth_sat
	ch4_ebul_depth_sat 			=> cmic%ch4_ebul_depth_sat
	co2_prod_ace_depth_sat 		=> cmic%co2_prod_ace_depth_sat
	co2_decomp_depth_sat 		=> cmic%co2_decomp_depth_sat
	co2_cons_depth_sat 			=> cmic%co2_cons_depth_sat
	co2_ebul_depth_sat 			=> cmic%co2_ebul_depth_sat
	co2_aere_depth_sat 			=> cmic%co2_aere_depth_sat
	co2_dif_depth_sat 			=> cmic%co2_dif_depth_sat
	o2_cons_depth_sat 			=> cmic%o2_cons_depth_sat
	o2_aere_depth_sat 			=> cmic%o2_aere_depth_sat
	o2_aere_oxid_depth_sat 		=> cmic%o2_aere_oxid_depth_sat
	o2_decomp_depth_sat 			=> cmic%o2_decomp_depth_sat
	o2_cons_depth_sat 			=> cmic%o2_cons_depth_sat
	o2_dif_depth_sat 			=> cmic%o2_dif_depth_sat
	h2_prod_depth_sat 			=> cmic%h2_prod_depth_sat
	h2_cons_depth_sat 			=> cmic%h2_cons_depth_sat
	h2_aere_depth_sat 			=> cmic%h2_aere_depth_sat
	h2_diff_depth_sat 			=> cmic%h2_diff_depth_sat
	h2_ebul_depth_sat 			=> cmic%h2_ebul_depth_sat
	ch4_surf_aere_sat			=> cmic%ch4_surf_aere_sat
	ch4_surf_ebul_sat 			=> cmic%ch4_surf_ebul_sat
	ch4_surf_dif_sat 				=> cmic%ch4_surf_dif_sat
	ch4_surf_netflux_sat 			=> cmic%ch4_surf_netflux_sat
	co2_surf_aere_sat			=> cmic%co2_surf_aere_sat
	co2_surf_ebul_sat 			=> cmic%co2_surf_ebul_sat
	co2_surf_dif_sat 				=> cmic%co2_surf_dif_sat
	co2_surf_netflux_sat 			=> cmic%co2_surf_netflux_sat
	o2_surf_aere_sat				=> cmic%o2_surf_aere_sat
	o2_surf_dif_sat 				=> cmic%o2_surf_dif_sat
	o2_surf_netflux_sat 			=> cmic%o2_surf_netflux_sat
	h2_surf_aere_sat				=> cmic%h2_surf_aere_sat
	h2_surf_ebul_sat 			=> cmic%h2_surf_ebul_sat
	h2_surf_dif_sat 				=> cmic%h2_surf_dif_sat
	h2_surf_netflux_sat 			=> cmic%h2_surf_netflux_sat

	ch4_prod_ace_depth 			=> cmic%ch4_prod_ace_depth
	ch4_prod_co2_depth 			=> cmic%ch4_prod_co2_depth
	ch4_oxid_o2_depth 			=> cmic%ch4_oxid_o2_depth
	ch4_oxid_aom_depth 			=> cmic%ch4_oxid_aom_depth
	ch4_aere_depth 				=> cmic%ch4_aere_depth
	ch4_dif_depth 				=> cmic%ch4_dif_depth
	ch4_ebul_depth 				=> cmic%ch4_ebul_depth
	co2_prod_ace_depth 			=> cmic%co2_prod_ace_depth
	co2_decomp_depth 			=> cmic%co2_decomp_depth
	co2_cons_depth 				=> cmic%co2_cons_depth
	co2_ebul_depth 				=> cmic%co2_ebul_depth
	co2_aere_depth 				=> cmic%co2_aere_depth
	co2_dif_depth 				=> cmic%co2_dif_depth
	o2_cons_depth 				=> cmic%o2_cons_depth
	o2_aere_depth 				=> cmic%o2_aere_depth
	o2_aere_oxid_depth 			=> cmic%o2_aere_oxid_depth
	o2_decomp_depth 			=> cmic%o2_decomp_depth
	o2_cons_depth 				=> cmic%o2_cons_depth
	o2_dif_depth 				=> cmic%o2_dif_depth
	h2_prod_depth 				=> cmic%h2_prod_depth
	h2_cons_depth 				=> cmic%h2_cons_depth
	h2_aere_depth 				=> cmic%h2_aere_depth
	h2_diff_depth 				=> cmic%h2_diff_depth
	h2_ebul_depth 				=> cmic%h2_ebul_depth
	ch4_surf_aere				=> cmic%ch4_surf_aere
	ch4_surf_ebul 				=> cmic%ch4_surf_ebul
	ch4_surf_dif 				=> cmic%ch4_surf_dif
	ch4_surf_netflux 				=> cmic%ch4_surf_netflux
	co2_surf_aere				=> cmic%co2_surf_aere
	co2_surf_ebul 				=> cmic%co2_surf_ebul
	co2_surf_dif 				=> cmic%co2_surf_dif
	co2_surf_netflux 				=> cmic%co2_surf_netflux
	o2_surf_aere				=> cmic%o2_surf_aere
	o2_surf_dif 				=> cmic%o2_surf_dif
	o2_surf_netflux 				=> cmic%o2_surf_netflux
	h2_surf_aere				=> cmic%h2_surf_aere
	h2_surf_ebul 				=> cmic%h2_surf_ebul
	h2_surf_dif 				=> cmic%h2_surf_dif
	h2_surf_netflux 				=> cmic%h2_surf_netflux
! end of code for gas transport added by xiaofeng 

	origionalsoilph   				=> cps%pH   
	soilpH_unsat               		=> cps%soilpH_unsat     
	soilpH_sat               			=> cps%soilpH_sat
	dz                    				=> cps%dz 
	zi                     				=> cps%zi
	z                      				=> cps%z
	ltype               				=> lun%itype
	ptype               				=> pft%itype
	clandunit       				=> col%landunit
	pcolumn       				=> pft%column
	cgridcell					=> col%gridcell
	wtcol       					=> pft%wtcol
	finundated					=> cws%finundated
!	fsat						=> cws%fsat  
	rootfr_vr          				=> pps%rootfr
	soiltemp 					=> ces%t_soisno
	fsat_pre					=> cmic%fsat_pre
!	tempavg_bgnpp     			=> p%pcf%tempavg_bgnpp
!	annavg_bgnpp      			=> p%pcf%annavg_bgnpp
	waterhead_unsat   			=> cmic%waterhead_unsat
   	
	sminn_vr                       		=> cns%sminn_vr
	hr_vr						=> ccf%hr_vr
	roothr					=> pcf%rr
	froot_mr					=> pcf%froot_mr
	froot_r					=> cmic%froot_r
	
	annsum_npp            			=> pepv%annsum_npp
	annavg_agnpp          			=> pcf%annavg_agnpp
	annavg_bgnpp          			=> pcf%annavg_bgnpp
	col_npp            				=> pcf_a%npp
	col_rr            				=> pcf_a%rr
	cannsum_npp 				=> cps%cannsum_npp

	decomp_cpools_vr              		=> ccs%decomp_cpools_vr
	decomp_npools_vr              		=> cns%decomp_npools_vr
   
	vwc           				=> cws%h2osoi_vol
	vwcsat           				=> cps%watsat

	sucsat                			=> cps%sucsat
	soilpsi               				=> cps%soilpsi
   
#if (defined HUM_HOL)
	qflx_lat_aqu   				=> cwf%qflx_lat_aqu
	wa   						=> cws%wa
	zwt            				=> cws%zwt
#endif

	caces_unsat_temp 			= 0._r8
	caces_sat_temp 				= 0._r8
	cdocs_unsat_temp 			= 0._r8
	cdocs_sat_temp 				= 0._r8
	dt = real( get_step_size(), r8 )
	
	rgasm					= rgas / 1000.
	do g =lbg, ubg
	forc_pch4(g) 				= atmch4 * forc_pbot(g)
	forc_po2(g) 				= atmo2 * forc_pbot(g)
	forc_pco2(g) 				= atmco2 * forc_pbot(g)
	forc_ph2(g) 				= atmh2 * forc_pbot(g)
	c_atm(g,1) 				= forc_pch4(g) / rgasm / forc_t(g)
	c_atm(g,2) 				= forc_po2(g) / rgasm / forc_t(g)
	c_atm(g,3) 				= forc_pco2(g) / rgasm / forc_t(g)
	c_atm(g,4) 				= forc_ph2(g) / rgasm / forc_t(g)
!write(iulog,*) " c_atm(g,2): ", c_atm(g,1), c_atm(g,2), c_atm(g,3), c_atm(g,4)
	end do

#if (defined HUM_HOL)
	finundated(c) = 0.99
#endif

!#ifdef MODELTEST
       do j = 1,nlevdecomp
            do fc = 1,num_soilc
               c = filter_soilc(fc)
	       
	if(j > jwaterhead_unsat(c)) then
		micfinundated = 0.99
	else
		micfinundated = finundated(c)
	end if  
!	      write(iulog,*) "microbial before: ",decomp_cpools_vr(c,j,i_dom), cdocs(c,j), decomp_npools_vr(c,j,i_dom), cdons(c,j)
      decomp_cpools_vr(c,j,i_dom) = max(0._r8, decomp_cpools_vr(c,j,i_dom))
      cdocs(c,j)			= decomp_cpools_vr(c,j,i_dom) 
      cdocs_unsat(c,j) 		= cdocs(c,j) !* (1. - micfinundated) ! concentration
      cdocs_sat(c,j) 		= cdocs(c,j) !* micfinundated ! concentration
      cdons(c,j) 			= decomp_npools_vr(c,j,i_dom)
      cdons_unsat(c,j) 	= cdons(c,j) !* (1. - micfinundated)  ! concentration
      cdons_sat(c,j) 		= cdons(c,j) !* micfinundated ! concentration
      cdons_min(c,j) 		= 0_r8
!      write(iulog,*) " before: ", decomp_cpools_vr(c,j,i_dom)
           end do
      end do
      
	do fc = 1,num_soilc
        c = filter_soilc(fc)
!write(iulog,*) " c: ", c
	l = clandunit(c)     
	if (ltype(l) == istsoil .or. ltype(l) == istcrop) then 
	do j = 1,nlevsoi
	cdocs_pre(c,j)				= cdocs(c,j)
	
	cdocs_unsat(c,j) 			= cdocs_unsat(c,j) / 12. 			!/ dz(c,j) gC/m3 -> mmmol/L or mol/m3
	caces_unsat(c,j) 			= caces_unsat(c,j) / 12. 			!/ dz(c,j)
	cacebios_unsat(c,j) 			= cacebios_unsat(c,j) / 12. 			!/ dz(c,j)
	cco2bios_unsat(c,j) 			= cco2bios_unsat(c,j) / 12. 			!/ dz(c,j)
	caerch4bios_unsat(c,j) 		= caerch4bios_unsat(c,j) / 12. 		!/ dz(c,j)
	canaerch4bios_unsat(c,j) 		= canaerch4bios_unsat(c,j) / 12. 		!/ dz(c,j)
	ccon_o2s_unsat(c,j) 			= ccon_o2s_unsat(c,j) / 32. 		!/ dz(c,j)
	ccon_ch4s_unsat(c,j) 		= ccon_ch4s_unsat(c,j) / 12. 		!/ dz(c,j)
	ccon_h2s_unsat(c,j) 			= ccon_h2s_unsat(c,j) / 2. 			!/ dz(c,j)
	ccon_co2s_unsat(c,j) 		= ccon_co2s_unsat(c,j) / 12. 		!/ dz(c,j)
!	write(iulog,*)"c ", c, " j ",j," ", ccon_o2s_unsat(c,j), " ",ccon_o2s_sat(c,j)
	
	cdocs_sat(c,j) 				= cdocs_sat(c,j) / 12. 			!/ dz(c,j)
	caces_sat(c,j) 				= caces_sat(c,j) / 12. 			!/ dz(c,j)
	cacebios_sat(c,j) 			= cacebios_sat(c,j) / 12. 			!/ dz(c,j)
	cco2bios_sat(c,j) 			= cco2bios_sat(c,j) / 12. 			!/ dz(c,j)
	caerch4bios_sat(c,j) 			= caerch4bios_sat(c,j) / 12. 		!/ dz(c,j)
	canaerch4bios_sat(c,j) 		= canaerch4bios_sat(c,j) / 12. 		!/ dz(c,j)
	ccon_o2s_sat(c,j) 			= ccon_o2s_sat(c,j) / 32. 			!/ dz(c,j)
	ccon_ch4s_sat(c,j) 			= ccon_ch4s_sat(c,j) / 12. 			!/ dz(c,j)
	ccon_h2s_sat(c,j) 			= ccon_h2s_sat(c,j) / 2. 			!/ dz(c,j)
	ccon_co2s_sat(c,j) 			= ccon_co2s_sat(c,j) / 12. 			!/ dz(c,j)	
	end do
	end if
	end do

	
! vertically distributed root respiration
	roothr_vr(:,:) = 0._r8
	rootfraction(:,:)=0.0_r8
	do fp = 1, num_soilc
		c = filter_soilc(fp)
		roothr_vr(c,:) = 0.0_r8
		rootfraction(c,:)=0.0_r8
		nppratio(c) = max(1e-9,col_rr(c)) * 1e+6 / max(cannsum_npp(c), 0.01)
!write(iulog, *) nppratio(c)," nppratio", col_npp(c), cannsum_npp(c)
		nppratio(c) = min(1.0, nppratio(c))
	end do
	do j=1,nlevsoi
		do fp = 1, num_soilp
		p = filter_soilp(fp)
		c = pcolumn(p)
		if (wtcol(p) > 0._r8 .and. ptype(p) /= noveg) then
		roothr_vr(c,j) = roothr_vr(c,j) + roothr(p)*rootfr_vr(p,j)*wtcol(p)
		froot_r(c,j) = froot_r(c,j) + froot_mr(p)*rootfr_vr(p,j)*wtcol(p)
		rootfraction(c,j) = rootfraction(c,j) + rootfr_vr(p,j)*wtcol(p)
		!write(iulog,*) "roothr_vr(c,j)",j, " j ", roothr_vr(c,j), "roothr(p)", roothr(p), "rootfr_vr(p,j)", rootfr_vr(p,j), "wtcol(p)", wtcol(p), "rootfraction(c,j)", rootfraction(c,j)
		ccon_co2s_sat(c,j) = ccon_co2s_sat(c,j) + froot_r(c,j) / 12.0 + hr_vr(c,j) / 12.0
		
               !~ anpp = annsum_npp(p) ! g C / m^2/yr
               !~ anpp = max(anpp, 0._r8) ! NPP can be negative b/c of consumption of storage pools
                  !~ if (annavg_agnpp(p) /= spval .and. annavg_bgnpp(p) /= spval .and. &
                    !~ annavg_agnpp(p) > 0._r8 .and. annavg_bgnpp(p) > 0._r8) then
                    !~ nppratio = annavg_bgnpp(p) / (annavg_agnpp(p) + annavg_bgnpp(p))
                  !~ else
                    !~ nppratio = 0.01_r8
                  !~ end if
                    !~ nppratio = max(0.01_r8, nppratio)		
	        end if
		end do
	end do
		
! end of root respiration
!	hr_vr(:,:)  some and litter decomposing c release
		
!~ do fc = 1,num_soilc
!~ c = filter_soilc(fc)
	!~ do j = 1,nlevsoi
!~ write(iulog,*)"microbe before diffusion: o2,ch4,h2,co2 ",ccon_o2s_sat(c,j),ccon_ch4s_sat(c,j),ccon_h2s_sat(c,j),ccon_co2s_sat(c,j)
	!~ end do
!~ end do


!	call seasonality(lbc, ubc, lbp, ubp, num_soilc, filter_soilc, num_soilp, filter_soilp)
	call get_waterhead(lbc, ubc, num_soilc, filter_soilc,jwaterhead_unsat)
	call gas_diffusion(lbc, ubc, num_soilc, filter_soilc)
#if (defined HUM_HOL)
	call lateral_bgc(lbc, ubc, num_soilc, filter_soilc)
#endif

!~ do fc = 1,num_soilc
!~ c = filter_soilc(fc)
	!~ do j = 1,nlevsoi
!~ write(iulog,*)"microbe after diffusion: o2,ch4,h2,co2 ",ccon_o2s_sat(c,j),ccon_ch4s_sat(c,j),ccon_h2s_sat(c,j),ccon_co2s_sat(c,j)
	!~ end do
!~ end do

	do fc = 1,num_soilc
        c = filter_soilc(fc)
	waterhead_unsat(c) = zi(c, jwaterhead_unsat(c))
!write(iulog,*)"water head: ",jwaterhead_unsat(c)
	end do

!	separate carbon respiration to ten soil layers in saturated and unsaturated fraction
	do fc = 1,num_soilc
        c = filter_soilc(fc)
		do j = 1, nlevsoi
		if(j > jwaterhead_unsat(c)) then
		micfinundated = 0.99
		co2_decomp_depth_sat(c,j) = (roothr_vr(c,j) + hr_vr(c,j)) !* micfinundated ! concentration
		co2_decomp_depth_unsat(c,j) = (roothr_vr(c,j) + hr_vr(c,j)) !* (1. - micfinundated) ! concentration
		
		o2_decomp_depth_unsat(c,j) = co2_decomp_depth_unsat(c,j)
		o2_decomp_depth_sat(c,j) = co2_decomp_depth_sat(c,j)
		else
		micfinundated = finundated(c)
		co2_decomp_depth_sat(c,j) = (roothr_vr(c,j) + hr_vr(c,j)) !* micfinundated ! concentration
		co2_decomp_depth_unsat(c,j) = (roothr_vr(c,j) + hr_vr(c,j)) !* (1. - micfinundated) ! concentration
		
		o2_decomp_depth_unsat(c,j) = co2_decomp_depth_unsat(c,j)
		o2_decomp_depth_sat(c,j) = co2_decomp_depth_sat(c,j)		
		end if
		end do
	end do


      ! unsaturated fraction
	! begins of oxygen transport from atmosphere to the soil
	! end of oxygen transport from atmosphere to the soil

	do fc = 1,num_soilc
        c = filter_soilc(fc)
	g = cgridcell(c)
	do j = 1,nlevsoi
	AceProd = 0.
	ACH2Prod = 0.
	ACCO2Prod = 0.
	H2CH4Prod = 0.
	H2AceProd = 0.
	H2Cons = 0.
	CH4PlantFlux = 0.
	CH4Ebull = 0.
	PlantO2Cons = 0.
	AerO2Cons = 0.
	CH4O2Cons = 0.
	dtempratio = 0.
	dtemph2 = 0.
	H2PlantFlux = 0.
	CH4Prod = 0.
	CH4Oxid = 0.
	O2PlantFlux = 0.
	CO2Prod = 0.
	CO2PlantFlux = 0.
	AOMCH4Oxid = 0.
	AceCons = 0.
	H2CO2Cons = 0.
	AceMethanogenGrowth = 0.
	AceMethanogenDying = 0.
	H2MethanogenGrowth = 0.
	H2MethanogenDying = 0.
	MethanotrophGrowth = 0.
	MethanotrophDying = 0.
	AOMMethanotrophGrowth = 0.
	AOMMethanotrophDying = 0.
	
	l=clandunit(c)
	if(ltype(l) == istwet) then
	ccon_o2s_unsat(c,j) = 0._r8
	endif

! below water table, larger than the number is below water table, all saturated
if(j >= jwaterhead_unsat(c)) then
!write(iulog,*) "c: ", c, " original soil ph: ", origionalsoilph(c), " ",  j, " j and water table ",  jwaterhead_unsat(c)
	if(origionalsoilph(c) > 5.5 .and. caces_unsat(c,j) > 0) then
	soilpH_unsat(c,j) = -1 * log10((10**(-origionalsoilph(c)) + 0.0042 * 1e-6 * caces_unsat(c,j)))
	else
	soilpH_unsat(c,j) = origionalsoilph(c)
	end if

	pHeffect = (soilpH_unsat(c,j) - pHmin) * (soilpH_unsat(c,j) - pHmax) / ((soilpH_unsat(c,j) - pHmin) &
		* (soilpH_unsat(c,j) - pHmax) - (soilpH_unsat(c,j) - pHopt) * (soilpH_unsat(c,j) - pHopt))
	pHeffect = min(1.0_r8, pHeffect)

	ACConcentration = 0.0
	ACConcentration = ACConcentration + cdocs_unsat(c,j) ! m mol C / m3
 !write(iulog,*) dt, "ACConcentration: ", ACConcentration, " cdocs(c,j): ", cdocs(c,j)

! Xiaofeng replaced the following conditional code with control by oxygen with the control with soil moiture
	!~ if(ccon_o2s_unsat(c,j) <= 1) then
	!~ AceProd = 2.0 / 3.0 * m_dAceProdACmax * ACConcentration / (ACConcentration + m_dKAce) &
		!~ * (m_dACMinQ10 ** ((soiltemp(c,j) - 286.65) / 10.)) * pHeffect &
		!~ * (1. - (caces_unsat(c,j) / (caces_unsat(c,j) + 100.)))
	!~ ACH2Prod = AceProd / 6.0
	!~ ACCO2Prod = 0.5 * AceProd
	!~ IsH2Production = 1
	!~ else
	!~ AceProd = 2.0 / 3.0 * m_dAceProdACmax * (1 - ccon_o2s_unsat(c,j) / (ccon_o2s_unsat(c,j) + m_dKAceProdO2)) &
		!~ * (m_dAceProdQ10 ** (soiltemp(c,j) - 286.65) / 10.) * pHeffect &
		!~ * (1. - (caces_unsat(c,j) / (caces_unsat(c,j) + 100.)))
	!~ ACCO2Prod = 0.5 * AceProd
	!~ ACH2Prod = 0
	!~ IsH2Production = 0
	!~ end if
	
	IsH2Production = 1
	minpsi = -10.0_r8;
	maxpsi = sucsat(c,j) * (-9.8e-6_r8)
	psi = min(soilpsi(c,j),maxpsi)
	    
	AceProd = 2.0 / 3.0 * m_dAceProdACmax * ACConcentration / (ACConcentration + m_dKAce) &
		* (m_dACMinQ10 ** ((soiltemp(c,j) - 286.65) / 10.)) * pHeffect &
		* (1. - (caces_unsat(c,j) / (caces_unsat(c,j) + 0.1))) * (log(minpsi/psi)/log(minpsi/maxpsi))
	ACH2Prod = AceProd / 6.0
	ACCO2Prod = 0.5 * AceProd
! Xiaofeng replaced the above conditional code with control by oxygen with the control with soil moiture

	if(AceProd < 0) then
	AceProd = 0
	end if

	if(ACConcentration > (3.0 / 2.0 * AceProd * dt)) then
	AceProd = AceProd
	else	
	AceProd = ACConcentration * 2.0 / 3.0 / dt
	ACCO2Prod = 0.5 * AceProd
	end if
	
!	caces_unsat_prod(c,j) = AceProd

	if(IsH2Production == 1.0)  then
	ACH2Prod = (AceProd / 6.0)
	ACCO2Prod = 0.5 * AceProd
	else
	ACCO2Prod = 0.5 * AceProd
	ACH2Prod = 0.0
	end if

        ccon_h2s_unsat(c,j) = ccon_h2s_unsat(c,j) + ACH2Prod * dt !
!write(iulog,*) "ccon_h2s(c,j): ",	ccon_h2s(c,j)	 
	ACConcentration = ACConcentration - (3.0 / 2.0 * AceProd * dt)
	
	if(ACConcentration < 0.0) then
	ACConcentration = 0.0
	end if
	
!	if(ccon_h2s_unsat(c,j) <= m_dCH4H2min) then
	! Xiaofeng replaced the following two lines code with new mechanism of CH4 production from CO2 
	!back to orgional on 7/11/2013
	H2AceProd = 0
	H2CH4Prod = 0
	! end
	!HCH4Prod = m_dGrowRH2Methanogens / m_dYH2Methanogens * cco2bios_unsat(c,j) &
	!	* ccon_co2s_unsat(c,j) / (ccon_co2s_unsat(c,j)  + m_dKCO2ProdCH4) &
	!	* (m_dH2CH4ProdQ10 ** ((soiltemp(c,j) - 286.65) / 10.)) * pHeffect / 10.
	!H2AceProd = 0	
	!print *, "HCH4Prod: ", HCH4Prod	
!	else
!	if(ccon_h2s_unsat(c,j) <= m_dAceH2min) then
	H2CH4Prod = m_dGrowRH2Methanogens / m_dYH2Methanogens * cco2bios_unsat(c,j) &
	* ccon_h2s_unsat(c,j) / ( ccon_h2s_unsat(c,j) + m_dKH2ProdCH4) &
		* ccon_co2s_unsat(c,j) / (ccon_co2s_unsat(c,j) + m_dKCO2ProdCH4) &
		* (m_dH2CH4ProdQ10 ** ((soiltemp(c,j) - 286.65) / 10.)) * pHeffect &
		* (1.0 - min(1.0, ccon_co2s_unsat(c,j) / 9.2))  !9.375 is the 21% oxgen
!	H2AceProd = 0
!write(iulog,*) "H2CH4Prod: ", H2CH4Prod
!	else
	H2AceProd = m_dH2ProdAcemax * ccon_h2s_unsat(c,j) / (ccon_h2s_unsat(c,j) + m_dKH2ProdAce) &
	* cco2bios_unsat(c,j) / (cco2bios_unsat(c,j) + m_dKCO2ProdAce) &
		* (m_dH2AceProdQ10 ** ((soiltemp(c,j) - 286.65) / 10.)) * pHeffect
!	H2CH4Prod = 0
!	end if
!	end if
			
	H2Cons = 4. * (H2AceProd + H2CH4Prod)

	H2Cons = max(0._r8, H2Cons)
	ccon_h2s_unsat(c,j) = max(0._r8, ccon_h2s_unsat(c,j))	
	
	if(ccon_h2s_unsat(c,j) >= (H2Cons * dt)) then
	H2Cons = H2Cons
	else
	dTempratio = ccon_h2s_unsat(c,j) / (H2Cons * dt)
	dTempH2 = ccon_h2s_unsat(c,j)
	H2Cons = ccon_h2s_unsat(c,j) / dt
!	ccon_h2s_unsat(c,j) = 0
	H2AceProd = H2AceProd * dTempratio
	H2CH4Prod = H2CH4Prod * dTempratio ! dTempH2 - H2AceProd * dt
	end if
		
	ccon_h2s_unsat(c,j) = ccon_h2s_unsat(c,j) - H2Cons * dt
	ccon_h2s_unsat(c,j) = max(0._r8, ccon_h2s_unsat(c,j))
	H2AceProd = max(0._r8, H2AceProd)
	H2CH4Prod = max(0._r8, H2CH4Prod)
	
	if(ccon_h2s_unsat(c,j)>g_dMaxH2inWater) then
	H2PlantFlux = m_dPlantTrans *  rootfraction(c,j) * (ccon_h2s_unsat(c,j) - g_dMaxH2inWater) * nppratio(c)*exp(-z(c,j)/0.25)/z(c,j) !* bgnpp_timestep(c) / bgnpp_avg(c)
	else
	H2PlantFlux = 0._r8
	end if
	
	ccon_h2s_unsat(c,j) = ccon_h2s_unsat(c,j) - H2PlantFlux * dt
	ccon_h2s_unsat(c,j) = max(0._r8, ccon_h2s_unsat(c,j))
	
	caces_unsat(c,j) = caces_unsat(c,j) + (AceProd + H2AceProd) * dt
	
!write(iulog,*) "here", cacebios_unsat(c,j)
	!	// For acetate dyndamics
	AceCons = m_dGrowRAceMethanogens / m_dYAceMethanogens * cacebios_unsat(c,j) &
	* caces_unsat(c,j)  / (caces_unsat(c,j)  + m_dKCH4ProdAce) &
		* (m_dCH4ProdQ10 ** ((soiltemp(c,j) - 286.65) / 10.)) * pHeffect !* 1. / (1. + ccon_co2s_unsat(c,j))
	!else
!write(iulog,*)"acecons: ", m_dGrowRAceMethanogens, m_dYAceMethanogens, cacebios(c,j), caces(c,j), m_dKCH4ProdAce, m_dCH4ProdQ10, AceCons, ccon_co2s(c,j) 
	!endif
	
	if((caces_unsat(c,j) - (AceCons * dt))>0) then
	AceCons = AceCons
	else
	AceCons = caces_unsat(c,j) / dt
	end if

	caces_unsat(c,j) = caces_unsat(c,j) - AceCons * dt

	caces_unsat(c,j) = max(0._r8,caces_unsat(c,j))

	!// For CH4 dyndamics
!write(iulog,*) "m_drCH4Prod * (1 - m_dYAceMethanogens): ", m_drCH4Prod, " ", (1 - m_dYAceMethanogens) * AceCons, "H2CH4Prod: ", H2CH4Prod
	CH4Prod = m_drCH4Prod * (1 - m_dYAceMethanogens) * AceCons + H2CH4Prod + HCH4Prod
	
	CH4Prod = max(0._r8, CH4Prod)
	
	ccon_ch4s_unsat(c,j) = ccon_ch4s_unsat(c,j) + CH4Prod * dt

	CH4Oxid = m_dGrowRMethanotrophs / m_dYMethanotrophs * caerch4bios_unsat(c,j) &
	* ccon_ch4s_unsat(c,j) / (ccon_ch4s_unsat(c,j) + m_dKCH4OxidCH4) &
		* ccon_o2s_unsat(c,j) / (ccon_o2s_unsat(c,j) + m_dKCH4OxidO2) * (m_dCH4OxidQ10 ** ((soiltemp(c,j) - 286.65) / 10.)) * pHeffect
		
	if(CH4Oxid > 0) then
	CH4Oxid = CH4Oxid
	else
	CH4Oxid = 0
	end if

	if(ccon_ch4s_unsat(c,j) > (CH4Oxid * dt) .and. ccon_o2s_unsat(c,j)>(2.0*CH4Oxid*dt)) then
	CH4Oxid = CH4Oxid
		else
		if(ccon_ch4s_unsat(c,j)<=0 .or. ccon_o2s_unsat(c,j)<=0) then
		CH4Oxid = 0
		else
		CH4Oxid = min(ccon_ch4s_unsat(c,j) / dt, ccon_o2s_unsat(c,j) / 2.0 / dt)
		end if
	end if

	if(CH4Oxid < 0) then
	CH4Oxid = 0
	end if

	ccon_ch4s_unsat(c,j) = ccon_ch4s_unsat(c,j) - CH4Oxid * dt

	ccon_ch4s_unsat(c,j) = max(0._r8, ccon_ch4s_unsat(c,j))

	if(ccon_ch4s_unsat(c,j) < 0._r8) then
	ccon_ch4s_unsat(c,j) = 0._r8
	end if
	
!write(iulog,*) "CH4PlantFlux: ", m_dPlantTrans, " ", rootfr_vr(c,j), ccon_ch4s_unsat(c,j), m_dCH4min !, tempavg_bgnpp(c), annavg_bgnpp(c)
	if(soiltemp(c,jwaterhead_unsat(c)) > -0.1 .and. ccon_ch4s_unsat(c,j) > m_dCH4min/1.0) then
	CH4PlantFlux = m_dPlantTrans *  rootfraction(c,j) * (ccon_ch4s_unsat(c,j) - m_dCH4min / 1.0) * nppratio(c)*exp(-z(c,j)/0.25)/z(c,j)  !* bgnpp_timestep(c) / bgnpp_avg(c)		
!	ccon_ch4s_unsat(c,j) = ccon_ch4s_unsat(c,j) - CH4PlantFlux
	CH4Ebull = max((ccon_ch4s_unsat(c,j) - m_dCH4min), 0._r8) * nppratio(c) * exp(-z(c,j)/0.35) !* (2.0 ** ((soiltemp(c,j) - 286.65) / 10.)) ! (exp(0.5 * (j - 0.5))-1) !* exp(-z(c,j)/1.25)   ! mmol/L     ! current these two equation have same threshold, will need to be corrected later
	else
	CH4PlantFlux = 0.0 * m_dPlantTrans *  rootfraction(c,j)*(ccon_ch4s_unsat(c,j) - m_dCH4min / 1.0) * nppratio(c)*exp(-z(c,j)/0.25)/z(c,j)
	!CH4Ebull = 0._r8
	endif

	CH4PlantFlux = max(CH4PlantFlux, 0._r8)
	CH4Ebull = max(CH4Ebull, 0._r8)
	CH4Ebull = min(CH4Ebull, ccon_ch4s_unsat(c,j) / 2.0)
	
	ccon_ch4s_unsat(c,j) = ccon_ch4s_unsat(c,j) - CH4PlantFlux * dt - CH4Ebull * dt
		
	!	// For O2 dyndamics
	AerO2Cons = m_drAer * ACCO2Prod
	CH4O2Cons = m_drCH4Oxid * CH4Oxid
	ccon_o2s_unsat(c,j) = ccon_o2s_unsat(c,j) - (AerO2Cons + CH4O2Cons) * dt
	
	if(ccon_o2s_unsat(c,j) < 0._r8) then
	ccon_o2s_unsat(c,j) = 0._r8
	end if

	O2PlantFlux = m_dPlantTrans*rootfraction(c,j)*(ccon_o2s_unsat(c,j) - c_atm(g,2))*nppratio(c)*exp(-z(c,j)/0.15)/z(c,j)  !* bgnpp_timestep(c) / bgnpp_avg(c) * exp(-z(c,j) / 0.1) 
	if(ccon_o2s_unsat(c,j) > c_atm(g,2) .or. soiltemp(c,jwaterhead_unsat(c)) < -0.1) then
	O2PlantFlux = 0._r8
	else
	O2PlantFlux = O2PlantFlux !* dz(c,j)
	!O2PlantFlux = max(O2PlantFlux, (ccon_o2s_unsat(c,j) - c_atm(g,2))) * dz(c,j)
	end if
	
	ccon_o2s_unsat(c,j) = ccon_o2s_unsat(c,j) - O2PlantFlux * dt !/ dz(c,j)
	
	ccon_o2s_unsat(c,j) = min(ccon_o2s_unsat(c,j), c_atm(g,2))

	PlantO2Cons = O2PlantFlux * 0.001
	
	if(ccon_o2s_unsat(c,j) >= (PlantO2Cons*dt)) then
	PlantO2Cons = PlantO2Cons
	else
	PlantO2Cons = ccon_o2s_unsat(c,j)/dt
	end if
		
	ccon_o2s_unsat(c,j) = ccon_o2s_unsat(c,j) - PlantO2Cons * dt
	
	ccon_o2s_unsat(c,j) = max(0.0, ccon_o2s_unsat(c,j) - o2_decomp_depth_unsat(c,j) * dt)
		
	if(ccon_o2s_unsat(c,j) < 0._r8) then
	ccon_o2s_unsat(c,j) = 0._r8
	end if
	
	!// For CO2 dyndamics
!write(iulog,*) "ACCO2Prod: ", ACCO2Prod, "CH4Prod: ", m_drCH4Prod * (1 - m_dYAceMethanogens) * AceCons, "CH4Oxid: ", CH4Oxid
	CO2Prod = ACCO2Prod + m_drCH4Prod * (1 - m_dYAceMethanogens) * AceCons + CH4Oxid + PlantO2Cons
	H2CO2Cons = 2 * H2AceProd + H2CH4Prod + HCH4Prod
	
	ccon_co2s_unsat(c,j) = ccon_co2s_unsat(c,j) + (CO2Prod - H2CO2Cons) * dt
	
	CO2PlantFlux = ccon_co2s_unsat(c,j) * 0.001
		
	ccon_co2s_unsat(c,j) = ccon_co2s_unsat(c,j) - CO2PlantFlux * dt

	if(ccon_co2s_unsat(c,j) < 0) then
	ccon_co2s_unsat(c,j) = 0
	end if
	
	!print *, "h2: ", m_dH2,  "O2: ", m_dO2,  "CO2: ",m_dCO2
	!if((ccon_h2s(c,j) * 0.5e-7) < 0.001 .and. ccon_co2s(c,j) < 0.00001 .and. ccon_co2s(c,j) < 0.001) then
	!AOM = ccon_ch4s(c,j) * 0.001 * (hco3 + 1.) / (hco3 + 0.1)
	!AOM = 0.085 * ccon_ch4s(c,j) / (ccon_ch4s(c,j) + 37.)
	!ccon_ch4s(c,j) = ccon_ch4s(c,j) - AOM
	!hco3 = hco3 + AOM
	!end if
	
!	if(ccon_o2s(c,j) < 0.0001) then
	AOMCH4Oxid = m_dGrowRAOMMethanotrophs / m_dYAOMMethanotrophs * canaerch4bios_unsat(c,j) * &
		ccon_ch4s_unsat(c,j)  / (ccon_ch4s_unsat(c,j)  + m_dKAOMCH4OxidCH4) * (m_dAOMCH4OxidQ10 ** ((soiltemp(c,j) - 13.5) / 10.)) &
		* (1.0 - min(1.0, ccon_o2s_unsat(c,j) / 4.6)) * pHeffect
!	endif
!	
	if((AOMCH4Oxid *dt)< ccon_ch4s_unsat(c,j)) then
	AOMCH4Oxid = AOMCH4Oxid
	else
	AOMCH4Oxid = ccon_ch4s_unsat(c,j)/dt
	end if	
		
	ccon_ch4s_unsat(c,j) = ccon_ch4s_unsat(c,j) - AOMCH4Oxid*dt
	ccon_ch4s_unsat(c,j) = max(0._r8, ccon_ch4s_unsat(c,j))
	
!	// For Microbe dyndamics
	AceMethanogenGrowth = m_dYAceMethanogens * AceCons * pHeffect
	AceMethanogenDying = m_dDeadRAceMethanogens * cacebios_unsat(c,j) * pHeffect
	AceMethanogenGrowth = max(0._r8, AceMethanogenGrowth)
	AceMethanogenDying = max(0._r8, AceMethanogenDying)

	!print *, AceMethanogenGrowth," ",AceMethanogenDying
	H2MethanogenGrowth = m_dYH2Methanogens * 4. * H2CH4Prod * pHeffect
	H2MethanogenDying = m_dDeadRH2Methanogens * cco2bios_unsat(c,j) * pHeffect
	H2MethanogenGrowth = max(0._r8, H2MethanogenGrowth)
	H2MethanogenDying = max(0._r8, H2MethanogenDying)

	MethanotrophGrowth = m_dYMethanotrophs * CH4Oxid * pHeffect
	MethanotrophDying = m_dDeadRMethanotrophs * caerch4bios_unsat(c,j) * pHeffect
	MethanotrophGrowth = max(0._r8, MethanotrophGrowth)
	MethanotrophDying = max(0._r8, MethanotrophDying)

	AOMMethanotrophGrowth = m_dYAOMMethanotrophs * AOMCH4Oxid * pHeffect
	AOMMethanotrophDying = m_dDeadRAOMMethanotrophs * canaerch4bios_unsat(c,j) * pHeffect
	AOMMethanotrophGrowth = max(0._r8, AOMMethanotrophGrowth)
	AOMMethanotrophDying = max(0._r8, AOMMethanotrophDying)

	!~ cacebios_unsat(c,j) = cacebios_unsat(c,j) + min((AceMethanogenGrowth - AceMethanogenDying), 0._r8)
	!~ cco2bios_unsat(c,j) = cco2bios_unsat(c,j) + min((H2MethanogenGrowth - H2MethanogenDying), 0._r8)
	!~ caerch4bios_unsat(c,j) = caerch4bios_unsat(c,j) + min((MethanotrophGrowth - MethanotrophDying), 0._r8)
	!~ canaerch4bios_unsat(c,j) = canaerch4bios_unsat(c,j) + min((AOMMethanotrophGrowth - AOMMethanotrophDying), 0._r8)
	
	cacebios_unsat(c,j) = cacebios_unsat(c,j) + (AceMethanogenGrowth - AceMethanogenDying) * dt
	cco2bios_unsat(c,j) = cco2bios_unsat(c,j) + (H2MethanogenGrowth - H2MethanogenDying) * dt
	caerch4bios_unsat(c,j) = caerch4bios_unsat(c,j) + (MethanotrophGrowth - MethanotrophDying) * dt
	canaerch4bios_unsat(c,j) = canaerch4bios_unsat(c,j) + (AOMMethanotrophGrowth - AOMMethanotrophDying) * dt
	
	cacebios_unsat(c,j) = max(cacebios_unsat(c,j),MFGbiomin)
	cco2bios_unsat(c,j) = max(cco2bios_unsat(c,j),MFGbiomin)
	caerch4bios_unsat(c,j) = max(caerch4bios_unsat(c,j),MFGbiomin)
	canaerch4bios_unsat(c,j) = max(canaerch4bios_unsat(c,j),MFGbiomin)
		
!write(iulog,*) "xiaofeng here ", c,j,cacebios_unsat(c,j),cco2bios_unsat(c,j),caerch4bios_unsat(c,j),canaerch4bios_unsat(c,j)
	cdocs_unsat(c,j) = ACConcentration
	caces_unsat_prod(c,j)					= AceProd
	caces_unsat_prod_h2(c,j)				= H2AceProd

	ch4_prod_ace_depth_unsat(c,j) 			= CH4Prod 
	ch4_prod_co2_depth_unsat(c,j) 			= H2CH4Prod
	ch4_oxid_o2_depth_unsat(c,j) 			= CH4Oxid
	ch4_oxid_aom_depth_unsat(c,j) 			= AOMCH4Oxid
	ch4_aere_depth_unsat(c,j) 				= CH4PlantFlux
	ch4_dif_depth_unsat(c,j) 				= 0._r8
	ch4_ebul_depth_unsat(c,j) 				= CH4Ebull

	co2_prod_ace_depth_unsat(c,j) 			= ACCO2Prod
!	co2_decomp_depth_unsat(c,j) 		! this has been calclated in previous part of this subroutine
	co2_cons_depth_unsat(c,j) 				= 0._r8
	co2_ebul_depth_unsat(c,j) 				= 0._r8
	co2_aere_depth_unsat(c,j) 				= 0._r8
	co2_dif_depth_unsat(c,j) 				= 0._r8

	o2_cons_depth_unsat(c,j) 				= AerO2Cons + CH4O2Cons
	o2_aere_depth_unsat(c,j) 				= O2PlantFlux
	o2_aere_oxid_depth_unsat(c,j) 			= PlantO2Cons
	o2_decomp_depth_unsat(c,j) 			= co2_decomp_depth_unsat(c,j)
	o2_dif_depth_unsat(c,j) 				= 0._r8

	h2_prod_depth_unsat(c,j) 				= ACH2Prod
	h2_cons_depth_unsat(c,j) 				= H2Cons
	h2_aere_depth_unsat(c,j) 				= H2PlantFlux
	h2_diff_depth_unsat(c,j) 				= 0._r8
	h2_ebul_depth_unsat(c,j) 				= 0._r8
else
! above water table in unsaturated fraction of soil column
	if(origionalsoilph(c) > 5.5 .and. caces_unsat(c,j) > 0) then
	soilpH_unsat(c,j) = -1 * log10((10**(-origionalsoilph(c)) + 0.0042 * 1e-6 * caces_unsat(c,j)))
	else
	soilpH_unsat(c,j) = origionalsoilph(c)
	end if

	pHeffect = (soilpH_unsat(c,j) - pHmin) * (soilpH_unsat(c,j) - pHmax) / ((soilpH_unsat(c,j) - pHmin) &
		* (soilpH_unsat(c,j) - pHmax) - (soilpH_unsat(c,j) - pHopt) * (soilpH_unsat(c,j) - pHopt))
	pHeffect = min(1.0_r8, pHeffect)

	ACConcentration = 0.0
	ACConcentration = ACConcentration + cdocs_unsat(c,j) ! m mol C / m3
	
!	ccon_ch4s_unsat(c,j)	= c_atm(g,1)
!	ccon_o2s_unsat(c,j)	= c_atm(g,2)
!	ccon_co2s_unsat(c,j)	= c_atm(g,3)
!	ccon_h2s_unsat(c,j)	= c_atm(g,4)
	
 !~ write(iulog,*) "j: ", j, " o2: ", ccon_o2s_unsat(c,j)

! Xiaofeng replaced the following conditional code with control by oxygen with the control with soil moiture
	AceProd = 0.
	ACCO2Prod = 0.5 * AceProd
	ACH2Prod = 0
	IsH2Production = 0
	if(AceProd < 0) then
	AceProd = 0
	end if
	
	!~ IsH2Production = 1
	!~ minpsi = -10.0_r8;
	!~ maxpsi = sucsat(c,j) * (-9.8e-6_r8)
	!~ psi = min(soilpsi(c,j),maxpsi)
	    
	!~ AceProd = 2.0 / 3.0 * m_dAceProdACmax * ACConcentration / (ACConcentration + m_dKAce) &
		!~ * (m_dACMinQ10 ** ((soiltemp(c,j) - 286.65) / 10.)) * pHeffect &
		!~ * (1. - (caces_sat(c,j) / (caces_sat(c,j) + 100.))) * (log(minpsi/psi)/log(minpsi/maxpsi))
	!~ ACH2Prod = AceProd / 6.0
	!~ ACCO2Prod = 0.5 * AceProd
! Xiaofeng replaced the above conditional code with control by oxygen with the control with soil moiture

	if(ACConcentration > (3.0 / 2.0 * AceProd * dt)) then
	AceProd = AceProd
	else	
	AceProd = ACConcentration * 2.0 / 3.0 / dt
	end if
	
!	caces_unsat_prod(c,j) = AceProd

	if(IsH2Production == 1.0)  then
	ACH2Prod = (AceProd / 6.0)
	ACCO2Prod = 0.5 * AceProd
	else
	ACCO2Prod = 0.5 * AceProd
	ACH2Prod = 0.0
	end if

    ccon_h2s_unsat(c,j) = ccon_h2s_unsat(c,j) + ACH2Prod * dt !
!write(iulog,*) "ccon_h2s(c,j): ",	ccon_h2s(c,j)	 
	ACConcentration = ACConcentration - (3.0 / 2.0 * AceProd) * dt
	
	if(ACConcentration < 0.0) then
	ACConcentration = 0.0
	end if
	
	!~ if(ccon_h2s_unsat(c,j) <= m_dCH4H2min) then
	!~ H2AceProd = 0
	!~ H2CH4Prod = 0
	!~ else
	!~ if(ccon_h2s_unsat(c,j) <= m_dAceH2min) then
	!~ H2CH4Prod = m_dGrowRH2Methanogens / m_dYH2Methanogens * cco2bios_unsat(c,j) &
	!~ * ccon_h2s_unsat(c,j) / ( ccon_h2s_unsat(c,j) + m_dKH2ProdCH4) &
		!~ * ccon_co2s_unsat(c,j) / (ccon_co2s_unsat(c,j) + m_dKCO2ProdCH4) &
		!~ * (m_dH2CH4ProdQ10 ** ((soiltemp(c,j) - 286.65) / 10.)) * pHeffect
	!~ H2AceProd = 0
	!~ !print *, "H2CH4Prod: ", H2CH4Prod
	!~ else
	!~ H2AceProd = m_dH2ProdAcemax * ccon_h2s_unsat(c,j) / (ccon_h2s_unsat(c,j) + m_dKH2ProdAce) &
	!~ * cco2bios_unsat(c,j) / (cco2bios_unsat(c,j) + m_dKCO2ProdAce) &
		!~ * (m_dH2AceProdQ10 ** ((soiltemp(c,j) - 286.65) / 10.)) * pHeffect
	!~ H2CH4Prod = 0

	!~ end if
	!~ end if
		
	H2AceProd = 0
	H2CH4Prod = 0
	
	H2Cons = 4. * (H2AceProd + H2CH4Prod)

	H2Cons = max(0._r8, H2Cons)
	ccon_h2s_unsat(c,j) = max(0._r8, ccon_h2s_unsat(c,j))	
	
	if(ccon_h2s_unsat(c,j) >= (H2Cons * dt)) then
	H2Cons = H2Cons
	ccon_h2s_unsat(c,j) = ccon_h2s_unsat(c,j) - H2Cons * dt
	Else
	dTempratio = ccon_h2s_unsat(c,j) / (H2Cons * dt)
	dTempH2 = ccon_h2s_unsat(c,j)
!	H2Cons = ccon_h2s_unsat(c,j)
	ccon_h2s_unsat(c,j) = 0
	H2AceProd = H2AceProd * dTempratio
	H2CH4Prod = dTempH2 - H2AceProd * dt
	end if
		
	ccon_h2s_unsat(c,j) = max(0._r8, ccon_h2s_unsat(c,j))
	H2AceProd = max(0._r8, H2AceProd)
	H2CH4Prod = max(0._r8, H2CH4Prod)
	
	H2PlantFlux = 0._r8

	!	// For acetate dyndamics
	!~ AceCons = m_dGrowRAceMethanogens / m_dYAceMethanogens * cacebios_unsat(c,j) &
	!~ * caces_unsat(c,j)  / (caces_unsat(c,j)  + m_dKCH4ProdAce) &
		!~ * (m_dCH4ProdQ10 ** ((soiltemp(c,j) - 286.65) / 10.)) * pHeffect * 1. / (1. + ccon_co2s_unsat(c,j))
		
	AceCons = 0.
        OxidAce2CO2 = 0.05 * min(caces_unsat(c,j), ccon_o2s_unsat(c,j)) * ccon_o2s_unsat(c,j) / (ccon_o2s_unsat(c,j) + m_dKAceProdO2)
	! write(*,*) OxidAce2CO2, " c ", c, " j ",  j, " ace ", caces_unsat(c,j), " o2 ",ccon_o2s_unsat(c,j)
	if(caces_unsat(c,j) > ((AceCons + OxidAce2CO2)*dt)) then
	AceCons = AceCons
	else
	AceCons = (caces_unsat(c,j) - OxidAce2CO2 * dt) / dt
	end if

	caces_unsat(c,j) = caces_unsat(c,j) - (AceCons + OxidAce2CO2) * dt
	caces_unsat(c,j) = max(0._r8 ,caces_unsat(c,j))

	!// For CH4 dyndamics
!	write(iulog,*) "m_drCH4Prod * (1 - m_dYAceMethanogens): ", m_drCH4Prod, " ", (1 - m_dYAceMethanogens) * AceCons, "H2CH4Prod: ", H2CH4Prod
	CH4Prod = m_drCH4Prod * (1 - m_dYAceMethanogens) * AceCons + H2CH4Prod + HCH4Prod
	
	CH4Prod = max(0._r8, CH4Prod)
	
	ccon_ch4s_unsat(c,j) = ccon_ch4s_unsat(c,j) + CH4Prod * dt

	CH4Oxid = m_dGrowRMethanotrophs / m_dYMethanotrophs * caerch4bios_unsat(c,j) &
	* ccon_ch4s_unsat(c,j) / (ccon_ch4s_unsat(c,j) + m_dKCH4OxidCH4) &
		* ccon_o2s_unsat(c,j) / (ccon_o2s_unsat(c,j) + m_dKCH4OxidO2) * (m_dCH4OxidQ10 ** ((soiltemp(c,j) - 286.65) / 10.)) * pHeffect

	if(ccon_ch4s_unsat(c,j) > (CH4Oxid*dt) .and. ccon_o2s_unsat(c,j)>(2.0*CH4Oxid*dt)) then
	CH4Oxid = CH4Oxid
		else
		if(ccon_ch4s_unsat(c,j)<=0 .or. ccon_o2s_unsat(c,j)<=0) then
		CH4Oxid = 0
		else
		CH4Oxid = min(0.8 * ccon_ch4s_unsat(c,j)/dt, ccon_o2s_unsat(c,j) / 2.0/dt)
		end if
	end if

	CH4Oxid = max(0._r8 ,CH4Oxid)

	ccon_ch4s_unsat(c,j) = ccon_ch4s_unsat(c,j) - CH4Oxid * dt

	ccon_ch4s_unsat(c,j) = max(0._r8, ccon_ch4s_unsat(c,j))

	if(ccon_ch4s_unsat(c,j) < 0) then
	ccon_ch4s_unsat(c,j) = 0
	end if
	
	CH4PlantFlux = 0._r8
	CH4Ebull = 0._r8
		
	CH4PlantFlux = max(CH4PlantFlux, 0._r8)
	CH4Ebull = max(CH4Ebull, 0._r8)
		
	!	// For O2 dyndamics
	AerO2Cons = m_drAer * ACCO2Prod
	CH4O2Cons = m_drCH4Oxid * CH4Oxid
	ccon_o2s_unsat(c,j) = ccon_o2s_unsat(c,j) - (AerO2Cons + CH4O2Cons - OxidAce2CO2) * dt
	
	ccon_o2s_unsat(c,j) = max(0.0, ccon_o2s_unsat(c,j) - o2_decomp_depth_unsat(c,j))
	
	if(ccon_o2s_unsat(c,j) < 0) then
	ccon_o2s_unsat(c,j) = 0
	end if

	!// For CO2 dyndamics
!print *, "ACCO2Prod: ", ACCO2Prod, "CH4Prod: ", m_drCH4Prod * (1 - m_dYAceMethanogens) * AceCons, "CH4Oxid: ", CH4Oxid
	CO2Prod = ACCO2Prod + m_drCH4Prod * (1 - m_dYAceMethanogens) * AceCons + CH4Oxid + PlantO2Cons + OxidAce2CO2
	H2CO2Cons = 2 * H2AceProd + H2CH4Prod + HCH4Prod
	
	ccon_co2s_unsat(c,j) = ccon_co2s_unsat(c,j) + (CO2Prod - H2CO2Cons) * dt
	
	CO2PlantFlux = 0_r8

	if(ccon_co2s_unsat(c,j) < 0) then
	ccon_co2s_unsat(c,j) = 0
	end if
	
	!print *, "h2: ", m_dH2,  "O2: ", m_dO2,  "CO2: ",m_dCO2
	!if((ccon_h2s(c,j) * 0.5e-7) < 0.001 .and. ccon_co2s(c,j) < 0.00001 .and. ccon_co2s(c,j) < 0.001) then
	!AOM = ccon_ch4s(c,j) * 0.001 * (hco3 + 1.) / (hco3 + 0.1)
	!AOM = 0.085 * ccon_ch4s(c,j) / (ccon_ch4s(c,j) + 37.)
	!ccon_ch4s(c,j) = ccon_ch4s(c,j) - AOM
	!hco3 = hco3 + AOM
	!end if
	
!	if(ccon_o2s(c,j) < 0.0001) then
	AOMCH4Oxid = m_dGrowRAOMMethanotrophs / m_dYAOMMethanotrophs * canaerch4bios_unsat(c,j) * &
		ccon_ch4s_unsat(c,j)  / (ccon_ch4s_unsat(c,j)  + m_dKAOMCH4OxidCH4) * (m_dAOMCH4OxidQ10 ** ((soiltemp(c,j) - 13.5) / 10.)) &
		* (1. - min(1.0, ccon_o2s_unsat(c,j) / 4.6)) * pHeffect
!	endif
!	
	if((AOMCH4Oxid*dt) < ccon_ch4s_unsat(c,j)) then
	AOMCH4Oxid = AOMCH4Oxid
	else
	AOMCH4Oxid = ccon_ch4s_unsat(c,j)/dt
	end if	
		
	ccon_ch4s_unsat(c,j) = ccon_ch4s_unsat(c,j) - AOMCH4Oxid*dt
	ccon_ch4s_unsat(c,j) = max(0._r8, ccon_ch4s_unsat(c,j))
	
!	// For Microbe dyndamics
	AceMethanogenGrowth = m_dYAceMethanogens * AceCons * pHeffect
	AceMethanogenDying = m_dDeadRAceMethanogens * cacebios_unsat(c,j) * pHeffect
	AceMethanogenGrowth = max(0._r8, AceMethanogenGrowth)
	AceMethanogenDying = max(0._r8, AceMethanogenDying)

	!print *, AceMethanogenGrowth," ",AceMethanogenDying
	H2MethanogenGrowth = m_dYH2Methanogens * 4. * H2CH4Prod * pHeffect
	H2MethanogenDying = m_dDeadRH2Methanogens * cco2bios_unsat(c,j) * pHeffect
	H2MethanogenGrowth = max(0._r8, H2MethanogenGrowth)
	H2MethanogenDying = max(0._r8, H2MethanogenDying)

	MethanotrophGrowth = m_dYMethanotrophs * CH4Oxid * pHeffect
	MethanotrophDying = m_dDeadRMethanotrophs * caerch4bios_unsat(c,j) * pHeffect
	MethanotrophGrowth = max(0._r8, MethanotrophGrowth)
	MethanotrophDying = max(0._r8, MethanotrophDying)

	AOMMethanotrophGrowth = m_dYAOMMethanotrophs * AOMCH4Oxid * pHeffect
	AOMMethanotrophDying = m_dDeadRAOMMethanotrophs * canaerch4bios_unsat(c,j) * pHeffect
	AOMMethanotrophGrowth = max(0._r8, AOMMethanotrophGrowth)
	AOMMethanotrophDying = max(0._r8, AOMMethanotrophDying)

	!~ cacebios_unsat(c,j) = cacebios_unsat(c,j) + min((AceMethanogenGrowth - AceMethanogenDying), 0._r8)
	!~ cco2bios_unsat(c,j) = cco2bios_unsat(c,j) + min((H2MethanogenGrowth - H2MethanogenDying), 0._r8)
	!~ caerch4bios_unsat(c,j) = caerch4bios_unsat(c,j) + min((MethanotrophGrowth - MethanotrophDying), 0._r8)
	!~ canaerch4bios_unsat(c,j) = canaerch4bios_unsat(c,j) + min((AOMMethanotrophGrowth - AOMMethanotrophDying), 0._r8)

	cacebios_unsat(c,j) = cacebios_unsat(c,j) + (AceMethanogenGrowth - AceMethanogenDying)*dt
	cco2bios_unsat(c,j) = cco2bios_unsat(c,j) + (H2MethanogenGrowth - H2MethanogenDying)*dt
	caerch4bios_unsat(c,j) = caerch4bios_unsat(c,j) + (MethanotrophGrowth - MethanotrophDying)*dt
	canaerch4bios_unsat(c,j) = canaerch4bios_unsat(c,j) + (AOMMethanotrophGrowth - AOMMethanotrophDying)*dt
	
	cacebios_unsat(c,j) = max(cacebios_unsat(c,j),MFGbiomin)
	cco2bios_unsat(c,j) = max(cco2bios_unsat(c,j),MFGbiomin)
	caerch4bios_unsat(c,j) = max(caerch4bios_unsat(c,j),MFGbiomin)
	canaerch4bios_unsat(c,j) = max(canaerch4bios_unsat(c,j),MFGbiomin)
	
!write(iulog,*) "xiaofeng here2 ", c,j,cacebios_unsat(c,j),cco2bios_unsat(c,j),caerch4bios_unsat(c,j),canaerch4bios_unsat(c,j)	
	cdocs_unsat(c,j) = ACConcentration
	caces_unsat_prod(c,j)					= AceProd
	caces_unsat_prod_h2(c,j)				= H2AceProd
	
	ch4_prod_ace_depth_unsat(c,j) 			= CH4Prod 
	ch4_prod_co2_depth_unsat(c,j) 			= H2CH4Prod
	ch4_oxid_o2_depth_unsat(c,j) 			= CH4Oxid
	ch4_oxid_aom_depth_unsat(c,j) 			= AOMCH4Oxid
	ch4_aere_depth_unsat(c,j) 				= CH4PlantFlux
	ch4_dif_depth_unsat(c,j) 				= 0._r8
	ch4_ebul_depth_unsat(c,j) 				= CH4Ebull

	co2_prod_ace_depth_unsat(c,j) 			= ACCO2Prod
!	co2_decomp_depth_unsat(c,j) 		! this has been calclated in previous part of this subroutine
	co2_cons_depth_unsat(c,j) 				= 0._r8
	co2_ebul_depth_unsat(c,j) 				= 0._r8
	co2_aere_depth_unsat(c,j) 				= 0._r8
	co2_dif_depth_unsat(c,j) 				= 0._r8

	o2_cons_depth_unsat(c,j) 				= AerO2Cons + CH4O2Cons
	o2_aere_depth_unsat(c,j) 				= O2PlantFlux
	o2_aere_oxid_depth_unsat(c,j) 			= PlantO2Cons
	o2_decomp_depth_unsat(c,j) 				= co2_decomp_depth_unsat(c,j)
	o2_dif_depth_unsat(c,j) 				= 0._r8

	h2_prod_depth_unsat(c,j) 				= ACH2Prod
	h2_cons_depth_unsat(c,j) 				= H2Cons
	h2_aere_depth_unsat(c,j) 				= H2PlantFlux
	h2_diff_depth_unsat(c,j) 				= 0._r8
	h2_ebul_depth_unsat(c,j) 				= 0._r8
end if
! above water table in unsaturated fraction of soil column
       end do
! average to grid level fluxes or state variables    

	ch4_surf_aere_unsat(c)				= 0._r8
	ch4_surf_ebul_unsat(c) 				= 0._r8
	ch4_surf_dif_unsat(c) 					= 0._r8
	ch4_surf_netflux_unsat(c) 				= 0._r8

	co2_surf_aere_unsat(c)				= 0._r8
	co2_surf_ebul_unsat(c) 				= 0._r8
	co2_surf_dif_unsat(c) 					= 0._r8
	co2_surf_netflux_unsat(c) 				= 0._r8

	o2_surf_aere_unsat(c)					= 0._r8
	o2_surf_dif_unsat(c) 					= 0._r8
	o2_surf_netflux_unsat(c) 				= 0._r8

	h2_surf_aere_unsat(c)					= 0._r8
	h2_surf_ebul_unsat(c) 					= 0._r8
	h2_surf_dif_unsat(c) 					= 0._r8
	h2_surf_netflux_unsat(c) 				= 0._r8
!ch4, o2, co2, h2,
	tem1 = Henry_kHpc_w(1) * exp (Henry_C_w(1) * (1/soiltemp(c,1) - 1/kh_tbase))
	tem2 = Henry_kHpc_w(2) * exp (Henry_C_w(2) * (1/soiltemp(c,1) - 1/kh_tbase))
	tem3 = Henry_kHpc_w(3) * exp (Henry_C_w(3) * (1/soiltemp(c,1) - 1/kh_tbase))
	tem4 = Henry_kHpc_w(4) * exp (Henry_C_w(4) * (1/soiltemp(c,1) - 1/kh_tbase))
	
	tem1 = 1. / tem1 * atmch4 * 1000. / 16.
	tem2 = 1. / tem2 * atmo2 * 1000. / 32.
	tem3 = 1. / tem3 * atmco2 * 1000. / 44.
	tem4 = 1. / tem4 * atmh2 * 1000. / 2.
!write(iulog,*) "xiaofeng tem1: ", tem1, tem2, tem3, tem4
!write(iulog,*) "xiaofeng ccon_ch4s_unsat(c,1): ", ccon_ch4s_unsat(c,1), ccon_co2s_unsat(c,1),ccon_h2s_unsat(c,1),ccon_o2s_unsat(c,1)

if((soiltemp(c,jwaterhead_unsat(c)) < SHR_CONST_TKFRZ) .or. (jwaterhead_unsat(c) .lt. 9)) then
	ch4_surf_dif_unsat(c) 		= 0._r8
!	ch4_surf_aere_unsat(c) 	= 0._r8
	ch4_surf_ebul_unsat(c) 	= 0._r8
	
	h2_surf_dif_unsat(c) 		= 0._r8
!	h2_surf_aere_unsat(c) 	= 0._r8
	h2_surf_ebul_unsat(c) 		= 0._r8
	
	co2_surf_dif_unsat(c) 		= 0._r8
!	co2_surf_aere_unsat(c)	 = 0._r8
	co2_surf_ebul_unsat(c) 	= 0._r8
	
	o2_surf_dif_unsat(c) 		= 0._r8
!	o2_surf_aere_unsat(c) 	= 0._r8	
else
	ch4_surf_dif_unsat(c) 	= (ccon_ch4s_unsat(c,j) - tem1)/dt	! max(0._r8, (ccon_ch4s_unsat(c,1) - tem1))  ! solubility of ch4 is 0.0000227g/L
	ccon_ch4s_unsat(c,jwaterhead_unsat(c)) = ccon_ch4s_unsat(c,jwaterhead_unsat(c)) - ch4_surf_dif_unsat(c) * dt
	ch4_surf_dif_unsat(c) 	= ch4_surf_dif_unsat(c) * dz(c,jwaterhead_unsat(c))

	o2_surf_dif_unsat(c) 	= (ccon_o2s_unsat(c,jwaterhead_unsat(c)) - tem2)/dt		! max(0._r8, (ccon_o2s_unsat(c,1) - tem4))  ! solubility of ch4 is 0.00004/L
	ccon_o2s_unsat(c,jwaterhead_unsat(c)) = ccon_o2s_unsat(c,jwaterhead_unsat(c)) - o2_surf_dif_unsat(c) * dt
	o2_surf_dif_unsat(c) 	= o2_surf_dif_unsat(c) * dz(c,jwaterhead_unsat(c))

	co2_surf_dif_unsat(c) 	= (ccon_co2s_unsat(c,jwaterhead_unsat(c)) - tem3)/dt		! max(0._r8, (ccon_co2s_unsat(c,1) - tem3))  ! solubility of ch4 is 0.002/L
	ccon_co2s_unsat(c,jwaterhead_unsat(c)) = ccon_co2s_unsat(c,jwaterhead_unsat(c)) - co2_surf_dif_unsat(c) * dt
	co2_surf_dif_unsat(c) 	= co2_surf_dif_unsat(c) * dz(c,jwaterhead_unsat(c))

	h2_surf_dif_unsat(c) 	= (ccon_h2s_unsat(c,jwaterhead_unsat(c)) - tem4)/dt 		! max(0._r8, (ccon_h2s_unsat(c,1) - tem2))  ! solubility of ch4 is 0.0000015/L
	ccon_h2s_unsat(c,jwaterhead_unsat(c)) = ccon_h2s_unsat(c,jwaterhead_unsat(c)) - h2_surf_dif_unsat(c) * dt
	h2_surf_dif_unsat(c) 	= h2_surf_dif_unsat(c) * dz(c,jwaterhead_unsat(c))

	do j = 1,nlevsoi
	ch4_surf_aere_unsat(c)			= ch4_surf_aere_unsat(c) + ch4_aere_depth_unsat(c,j) * dz(c,j)
!	ch4_surf_ebul_unsat(c) 			= ch4_surf_ebul_unsat(c) + ch4_ebul_depth_unsat(c,j) * dz(c,j)
	
	co2_surf_aere_unsat(c)			= co2_surf_aere_unsat(c) + co2_aere_depth_unsat(c,j) * dz(c,j)
!	co2_surf_ebul_unsat(c) 			= co2_surf_ebul_unsat(c) + co2_ebul_depth_unsat(c,j) * dz(c,j)
	
	o2_surf_aere_unsat(c)				= o2_surf_aere_unsat(c) + o2_aere_depth_unsat(c,j) * dz(c,j)
	
	h2_surf_aere_unsat(c)				= h2_surf_aere_unsat(c) + h2_aere_depth_unsat(c,j) * dz(c,j)
!	h2_surf_ebul_unsat(c) 				= h2_surf_ebul_unsat(c) + h2_ebul_depth_unsat(c,j) * dz(c,j)
	end do
!end if
end if  ! end if of the frozen mechanism in trapping gases in soil

	ch4_surf_netflux_unsat(c) 			= ch4_surf_netflux_unsat(c) + ch4_surf_dif_unsat(c) + ch4_surf_aere_unsat(c) + ch4_surf_ebul_unsat(c)
	co2_surf_netflux_unsat(c) 			= co2_surf_netflux_unsat(c) + co2_surf_dif_unsat(c) + co2_surf_aere_unsat(c) + co2_surf_ebul_unsat(c)
	o2_surf_netflux_unsat(c) 			= o2_surf_netflux_unsat(c) + o2_surf_dif_unsat(c) + o2_surf_aere_unsat(c)
	h2_surf_netflux_unsat(c) 			= h2_surf_netflux_unsat(c) + h2_surf_dif_unsat(c) + h2_surf_aere_unsat(c) + h2_surf_ebul_unsat(c)
	
!write(iulog,*) "xiaofeng here! ", ch4_surf_netflux_unsat(c), ch4_surf_dif_unsat(c), ch4_surf_aere_unsat(c), ch4_surf_ebul_unsat(c)
	end do

	do fc = 1,num_soilc
        c = filter_soilc(fc)
	l = clandunit(c)     
	if (ltype(l) == istsoil .or. ltype(l) == istcrop) then 
	do j = 1,nlevsoi
	cdocs_unsat(c,j) 				= cdocs_unsat(c,j) * 12.			! convert from mmol/m3 to gC/m3
	caces_unsat(c,j) 				= caces_unsat(c,j) * 12.			! convert from mmol/m3 to gC/m3
	cacebios_unsat(c,j) 				= cacebios_unsat(c,j) * 12.		! convert from mmol/m3 to gC/m3
	cco2bios_unsat(c,j) 				= cco2bios_unsat(c,j) * 12. 		! convert from mmol/m3 to gC/m3
	caerch4bios_unsat(c,j) 			= caerch4bios_unsat(c,j) * 12.		! convert from mmol/m3 to gC/m3
	canaerch4bios_unsat(c,j) 			= canaerch4bios_unsat(c,j) * 12.	! convert from mmol/m3 to gC/m3
	ccon_o2s_unsat(c,j) 				= ccon_o2s_unsat(c,j) * 32.		! convert from mmol/m3 to gC/m3
	ccon_ch4s_unsat(c,j) 			= ccon_ch4s_unsat(c,j) * 12.		! convert from mmol/m3 to gC/m3
	ccon_h2s_unsat(c,j) 				= ccon_h2s_unsat(c,j) * 2.		! convert from mmol/m3 to gC/m3
	ccon_co2s_unsat(c,j) 			= ccon_co2s_unsat(c,j) * 12.		! convert from mmol/m3 to gC/m3
	end do
	end if
	end do
! end of the unsaturated simulation	


! Saturated fraction
      
	do fc = 1,num_soilc
        c = filter_soilc(fc)
	g = cgridcell(c)
	do j = 1,nlevsoi
	AceProd = 0._r8
	ACH2Prod = 0._r8
	ACCO2Prod = 0._r8
	H2CH4Prod = 0._r8
	H2AceProd = 0._r8
	H2Cons = 0._r8
	CH4PlantFlux = 0._r8
	CH4Ebull = 0._r8
	PlantO2Cons = 0._r8
	AerO2Cons = 0._r8
	CH4O2Cons = 0._r8
	dtempratio = 0._r8
	dtemph2 = 0._r8
	H2PlantFlux = 0._r8
	CH4Prod = 0._r8
	CH4Oxid = 0._r8
	O2PlantFlux = 0._r8
	CO2Prod = 0._r8
	CO2PlantFlux = 0._r8
	AOMCH4Oxid = 0._r8
	AceCons = 0._r8
	H2CO2Cons = 0._r8
	AceMethanogenGrowth = 0._r8
	AceMethanogenDying = 0._r8
	H2MethanogenGrowth = 0._r8
	H2MethanogenDying = 0._r8
	MethanotrophGrowth = 0._r8
	MethanotrophDying = 0._r8
	AOMMethanotrophGrowth = 0._r8
	AOMMethanotrophDying = 0._r8
	
	l=clandunit(c)
	if(ltype(l) == istwet) then
	ccon_o2s_unsat(c,j) = 0._r8
	endif

!write(iulog,*) "c: ", c, " original soil ph: ", origionalsoilph(c)
	if(origionalsoilph(c) > 5.5 .and. caces_sat(c,j) > 0) then
	soilpH_sat(c,j) = -1 * log10((10**(-origionalsoilph(c)) + 0.0042 * 1e-6 * caces_sat(c,j)))
	else
	soilpH_sat(c,j) = origionalsoilph(c)
	end if

	pHeffect = (soilpH_sat(c,j) - pHmin) * (soilpH_sat(c,j) - pHmax) / ((soilpH_sat(c,j) - pHmin) &
		* (soilpH_sat(c,j) - pHmax) - (soilpH_sat(c,j) - pHopt) * (soilpH_sat(c,j) - pHopt))
	pHeffect = min(1.0_r8, pHeffect)

	ACConcentration = 0.0
	ACConcentration = ACConcentration + cdocs_sat(c,j)

! Xiaofeng replaced the following conditional code with control by oxygen with the control with soil moiture
	!~ if(ccon_o2s_sat(c,j) <= 1) then
	!~ AceProd = 2.0 / 3.0 * m_dAceProdACmax * ACConcentration / (ACConcentration + m_dKAce) &
		!~ * (m_dACMinQ10 ** ((soiltemp(c,j) - 286.65) / 10.)) * pHeffect &
		!~ * (1. - (caces_sat(c,j) / (caces_sat(c,j) + 100.)))
	!~ ACH2Prod = AceProd / 6.0
	!~ ACCO2Prod = 0.5 * AceProd
	!~ IsH2Production = 1
	!~ else
	!~ AceProd = 2.0 / 3.0 * m_dAceProdACmax * (1 - ccon_o2s_sat(c,j) / (ccon_o2s_sat(c,j) + m_dKAceProdO2)) &
		!~ * (m_dAceProdQ10 ** (soiltemp(c,j) - 286.65) / 10.) * pHeffect &
		!~ * (1. - (caces_sat(c,j) / (caces_sat(c,j) + 100.)))
	!~ ACCO2Prod = 0.5 * AceProd
	!~ ACH2Prod = 0
	!~ IsH2Production = 0
	!~ end if
	
	IsH2Production = 1
	minpsi = -10.0_r8;
	maxpsi = sucsat(c,j) * (-9.8e-6_r8)
	psi = min(soilpsi(c,j),maxpsi)
	    
	AceProd = 2.0 / 3.0 * m_dAceProdACmax * ACConcentration / (ACConcentration + m_dKAce) &
		* (m_dACMinQ10 ** ((soiltemp(c,j) - 286.65) / 10.)) * pHeffect &
		* (1. - (caces_sat(c,j) / (caces_sat(c,j) + 0.1))) * (log(minpsi/psi)/log(minpsi/maxpsi))
	ACH2Prod = AceProd / 6.0
	ACCO2Prod = 0.5 * AceProd
! Xiaofeng replaced the above conditional code with control by oxygen with the control with soil moiture

	if(AceProd < 0) then
	AceProd = 0
	end if

	if(ACConcentration > (3.0 / 2.0 * AceProd*dt)) then
	AceProd = AceProd
	else	
	AceProd = ACConcentration * 2.0 / 3.0 /dt
	end if
	
!	caces_sat_prod(c,j) = AceProd

	if(IsH2Production == 1.0)  then
	ACH2Prod = (AceProd / 6.0)
	ACCO2Prod = 0.5 * AceProd
	else
	ACCO2Prod = 0.5 * AceProd
	ACH2Prod = 0.0
	end if

        ccon_h2s_sat(c,j) = ccon_h2s_sat(c,j) + ACH2Prod * dt!
!write(iulog,*) "ccon_h2s(c,j): ",	ccon_h2s(c,j)	 
	ACConcentration = ACConcentration - (3.0 / 2.0 * AceProd) * dt
	
	if(ACConcentration < 0.0) then
	ACConcentration = 0.0
	end if
	
	!print *, "H2: ", m_dH2, " ch4h2min: ", m_dCH4H2min, "CO2: ", m_dCO2
!	if(ccon_h2s_sat(c,j) <= m_dCH4H2min) then
	! Xiaofeng replaced the following two lines code with new mechanism of CH4 production from CO2 
	!back to orgional on 7/11/2013
	H2AceProd = 0
	H2CH4Prod = 0
	! end
	!HCH4Prod = m_dGrowRH2Methanogens / m_dYH2Methanogens * cco2bios_unsat(c,j) &
	!	* ccon_co2s_unsat(c,j) / (ccon_co2s_unsat(c,j)  + m_dKCO2ProdCH4) &
	!	* (m_dH2CH4ProdQ10 ** ((soiltemp(c,j) - 286.65) / 10.)) * pHeffect / 10.
	!H2AceProd = 0	
	!print *, "HCH4Prod: ", HCH4Prod	
!	else
!	if(ccon_h2s_sat(c,j) <= m_dAceH2min) then
	H2CH4Prod = m_dGrowRH2Methanogens / m_dYH2Methanogens * cco2bios_sat(c,j) &
	* ccon_h2s_sat(c,j) / ( ccon_h2s_sat(c,j) + m_dKH2ProdCH4) &
		* ccon_co2s_sat(c,j) / (ccon_co2s_sat(c,j) + m_dKCO2ProdCH4) &
		* (m_dH2CH4ProdQ10 ** ((soiltemp(c,j) - 286.65) / 10.)) * pHeffect&
		* (1.0 - min(1.0, ccon_co2s_unsat(c,j) / 9.2))
!	H2AceProd = 0
	!print *, "H2CH4Prod: ", H2CH4Prod
!	else
	H2AceProd = m_dH2ProdAcemax * ccon_h2s_sat(c,j) / (ccon_h2s_sat(c,j) + m_dKH2ProdAce) &
	* cco2bios_sat(c,j) / (cco2bios_sat(c,j) + m_dKCO2ProdAce) &
		* (m_dH2AceProdQ10 ** ((soiltemp(c,j) - 286.65) / 10.)) * pHeffect
!	H2CH4Prod = 0

!	end if
!	end if
			
	H2Cons = 4. * (H2AceProd + H2CH4Prod)

	H2Cons = max(0._r8, H2Cons)
	ccon_h2s_sat(c,j) = max(0._r8, ccon_h2s_sat(c,j))	
	
	if(ccon_h2s_sat(c,j) >= (H2Cons * dt)) then
	H2Cons = H2Cons
	H2AceProd = H2AceProd
	H2CH4Prod = H2CH4Prod
	else
	dTempratio = ccon_h2s_sat(c,j) / (H2Cons * dt)
	dTempH2 = ccon_h2s_sat(c,j)
	H2Cons = ccon_h2s_sat(c,j) / dt
	H2AceProd = H2AceProd * dTempratio
	H2CH4Prod = H2CH4Prod * dTempratio
	end if

	ccon_h2s_sat(c,j) = ccon_h2s_sat(c,j) - H2Cons * dt
		
	ccon_h2s_sat(c,j) = max(0._r8, ccon_h2s_sat(c,j))
	H2AceProd = max(0._r8, H2AceProd)
	H2CH4Prod = max(0._r8, H2CH4Prod)
	
	!print *, "H2cons: ", m_dH2
	!// For H2 dyndamics
	!if(m_dH2>g_dMaxH2inWater) then
	!H2PlantFlux = m_dPlantTrans * pow(RootDistribution, 0.5) * RootFactor * (m_dH2 - g_dMaxH2inWater)		
	!if((m_dH2 - H2PlantFlux)>0) then
	!H2PlantFlux = H2PlantFlux
	!else
	!H2PlantFlux = (m_dH2 - g_dMaxH2inWater)
	!end if
	!else
	if(ccon_h2s_sat(c,j)>g_dMaxH2inWater .and. soiltemp(c,jwaterhead_unsat(c)) > -0.1) then
	H2PlantFlux = m_dPlantTrans *  rootfraction(c,j) * (ccon_h2s_sat(c,j) - g_dMaxH2inWater) * nppratio(c)*exp(-z(c,j)/0.25)/z(c,j)  !* bgnpp_timestep(c) / bgnpp_avg(c)
	else
	H2PlantFlux = 0._r8
	end if
	
	ccon_h2s_sat(c,j) = ccon_h2s_sat(c,j) - H2PlantFlux * dt
	ccon_h2s_sat(c,j) = max(0._r8, ccon_h2s_sat(c,j))
	
	caces_sat(c,j) = caces_sat(c,j) + (AceProd + H2AceProd) * dt

	if(cacebios_sat(c,j)<1E-5) then 
		cacebios_sat(c,j) = 1E-5
	end if

	!	// For acetate dyndamics
	AceCons = m_dGrowRAceMethanogens / m_dYAceMethanogens * cacebios_sat(c,j) &
	* caces_sat(c,j)  / (caces_sat(c,j)  + m_dKCH4ProdAce) &
		* (m_dCH4ProdQ10 ** ((soiltemp(c,j) - 286.65) / 10.)) * pHeffect !* 1. / (1. + ccon_co2s_sat(c,j))
	!else
!	write(iulog,*)"acecons: ", m_dGrowRAceMethanogens, m_dYAceMethanogens, cacebios(c,j), caces(c,j), m_dKCH4ProdAce, m_dCH4ProdQ10, AceCons, ccon_co2s(c,j) 
	!endif
	
	if(caces_sat(c,j) > (AceCons * dt)) then
	AceCons = AceCons
	else
	AceCons = caces_sat(c,j) / dt
	end if

	caces_sat(c,j) = caces_sat(c,j) - AceCons * dt

	if(caces_sat(c,j) < 0) then
	caces_sat(c,j) = 0
	end if

	!// For CH4 dyndamics
!	write(iulog,*) "m_drCH4Prod * (1 - m_dYAceMethanogens): ", m_drCH4Prod, " ", (1 - m_dYAceMethanogens) * AceCons, "H2CH4Prod: ", H2CH4Prod
	CH4Prod = m_drCH4Prod * (1 - m_dYAceMethanogens) * AceCons + H2CH4Prod + HCH4Prod
	
	CH4Prod = max(0._r8, CH4Prod)
	
	ccon_ch4s_sat(c,j) = ccon_ch4s_sat(c,j) + CH4Prod * dt

	CH4Oxid = m_dGrowRMethanotrophs / m_dYMethanotrophs * caerch4bios_sat(c,j) &
	* ccon_ch4s_sat(c,j) / (ccon_ch4s_sat(c,j) + m_dKCH4OxidCH4) &
		* ccon_o2s_sat(c,j) / (ccon_o2s_sat(c,j) + m_dKCH4OxidO2) * (m_dCH4OxidQ10 ** ((soiltemp(c,j) - 286.65) / 10.)) * pHeffect
		
	if((ccon_ch4s_sat(c,j) > CH4Oxid * dt) .and. ccon_o2s_sat(c,j)>(2.0*CH4Oxid*dt)) then
	CH4Oxid = CH4Oxid
		else
		if(ccon_ch4s_sat(c,j)<=0 .or. ccon_o2s_sat(c,j)<=0) then
		CH4Oxid = 0
		else
		CH4Oxid = min(ccon_ch4s_sat(c,j)/dt, ccon_o2s_sat(c,j)/2.0/dt)
		end if
	end if

	CH4Oxid = max(0._r8, CH4Oxid)
	
	ccon_ch4s_sat(c,j) = ccon_ch4s_sat(c,j) - CH4Oxid * dt

	ccon_ch4s_sat(c,j) = max(0._r8, ccon_ch4s_sat(c,j))

	if(ccon_ch4s_sat(c,j) < 0) then
	ccon_ch4s_sat(c,j) = 0
	end if
	
	if(soiltemp(c,jwaterhead_unsat(c)) > -0.1 .and. ccon_ch4s_sat(c,j) > m_dCH4min/1.0) then
	CH4PlantFlux = m_dPlantTrans *  rootfraction(c,j) * (ccon_ch4s_sat(c,j) - m_dCH4min/1.0) * nppratio(c)*exp(-z(c,j)/0.25)/z(c,j)  !* bgnpp_timestep(c) / bgnpp_avg(c)		
	CH4Ebull = max((ccon_ch4s_sat(c,j) - m_dCH4min) /dt, 0._r8) * nppratio(c) * exp(-z(c,j)/0.35) !* (2.0 ** ((soiltemp(c,j) - 286.65) / 10.)) !/ (exp(0.5 * (j - 0.5))-1) * 20.0 !* exp(-z(c,j)/1.25) 				! mmol/L
	else
	CH4PlantFlux = 0.0 * m_dPlantTrans * rootfraction(c,j) * (ccon_ch4s_sat(c,j) - m_dCH4min/1.0) * nppratio(c)*exp(-z(c,j)/0.25)/z(c,j)
	!CH4Ebull = 0._r8
	endif
	
	CH4PlantFlux = max(CH4PlantFlux, 0._r8)
	CH4Ebull = max(CH4Ebull, 0._r8)
!	CH4Ebull = min(CH4Ebull, ccon_ch4s_sat(c,j) / 2.0)
	
	if((ccon_ch4s_sat(c,j) - m_dCH4min) > (CH4PlantFlux + CH4Ebull) * dt) then
	ccon_ch4s_sat(c,j) = ccon_ch4s_sat(c,j) - CH4PlantFlux * dt - CH4Ebull * dt
	else
	tem4 = (ccon_ch4s_sat(c,j) - m_dCH4min) / ((CH4PlantFlux + CH4Ebull) * dt)
	tem4 = max(0._r8 , tem4)
	ccon_ch4s_sat(c,j) = ccon_ch4s_sat(c,j) - (CH4PlantFlux + CH4Ebull) * tem4 * dt
	CH4PlantFlux = CH4PlantFlux * tem4
	CH4Ebull = CH4Ebull * tem4
	endif
	
	ccon_ch4s_sat(c,j) = max(0._r8 , ccon_ch4s_sat(c,j))

	!	// For O2 dyndamics
	AerO2Cons = m_drAer * ACCO2Prod
	CH4O2Cons = m_drCH4Oxid * CH4Oxid
	ccon_o2s_sat(c,j) = ccon_o2s_sat(c,j) - (AerO2Cons + CH4O2Cons) * dt
	
	if(ccon_o2s_sat(c,j) < 0) then
	ccon_o2s_sat(c,j) = 0
	end if

	if(soiltemp(c,jwaterhead_unsat(c)) > -0.1 .and. ccon_o2s_sat(c,j) < c_atm(g,2)) then
	O2PlantFlux = m_dPlantTrans * rootfraction(c,j) * (ccon_o2s_sat(c,j) - c_atm(g,2)) * nppratio(c)*exp(-z(c,j)/0.15)/z(c,j) 
	else
	O2PlantFlux = 0._r8
	!O2PlantFlux = O2PlantFlux !* dz(c,j)
	!O2PlantFlux = max(O2PlantFlux, (ccon_o2s_sat(c,j) - c_atm(g,2))) * dz(c,j)
	end if
	
	ccon_o2s_sat(c,j) = ccon_o2s_sat(c,j) - O2PlantFlux * dt !/ dz(c,j)
	
	ccon_o2s_sat(c,j) = min(ccon_o2s_sat(c,j), c_atm(g,2))

	PlantO2Cons = O2PlantFlux * 0.0001
	
	if(ccon_o2s_sat(c,j) >= (PlantO2Cons*dt)) then
	PlantO2Cons = PlantO2Cons
	else
	PlantO2Cons = 0.9 * ccon_o2s_sat(c,j) / dt
	end if
	
	ccon_o2s_sat(c,j) = ccon_o2s_sat(c,j) - PlantO2Cons * dt
	
	ccon_o2s_sat(c,j) = max(0.0, ccon_o2s_sat(c,j) - o2_decomp_depth_sat(c,j))
	
	if(ccon_o2s_sat(c,j) < 0) then
	ccon_o2s_sat(c,j) = 0
	end if
	
	!// For CO2 dyndamics
	CO2Prod = ACCO2Prod + m_drCH4Prod * (1 - m_dYAceMethanogens) * AceCons + CH4Oxid + PlantO2Cons
	H2CO2Cons = 2 * H2AceProd + H2CH4Prod + HCH4Prod
	
	ccon_co2s_sat(c,j) = ccon_co2s_sat(c,j) + (CO2Prod - H2CO2Cons) * dt
	
	CO2PlantFlux = ccon_co2s_sat(c,j) * 0.00001
		
	ccon_co2s_sat(c,j) = ccon_co2s_sat(c,j) - CO2PlantFlux * dt

	if(ccon_co2s_sat(c,j) < 0) then
	ccon_co2s_sat(c,j) = 0
	end if
	
	!print *, "h2: ", m_dH2,  "O2: ", m_dO2,  "CO2: ",m_dCO2
	!if((ccon_h2s(c,j) * 0.5e-7) < 0.001 .and. ccon_co2s(c,j) < 0.00001 .and. ccon_co2s(c,j) < 0.001) then
	!AOM = ccon_ch4s(c,j) * 0.001 * (hco3 + 1.) / (hco3 + 0.1)
	!AOM = 0.085 * ccon_ch4s(c,j) / (ccon_ch4s(c,j) + 37.)
	!ccon_ch4s(c,j) = ccon_ch4s(c,j) - AOM
	!hco3 = hco3 + AOM
	!end if
	
!	if(ccon_o2s(c,j) < 0.0001) then
	AOMCH4Oxid = m_dGrowRAOMMethanotrophs / m_dYAOMMethanotrophs * canaerch4bios_sat(c,j) * &
		ccon_ch4s_sat(c,j)  / (ccon_ch4s_sat(c,j)  + m_dKAOMCH4OxidCH4) * (m_dAOMCH4OxidQ10 ** ((soiltemp(c,j) - 13.5) / 10.)) &
		* (1.0 - min(1.0, ccon_o2s_sat(c,j) / 4.6))* pHeffect
!	endif
!	
	if((AOMCH4Oxid * dt) < ccon_ch4s_sat(c,j)) then
	AOMCH4Oxid = AOMCH4Oxid
	else
	AOMCH4Oxid = ccon_ch4s_sat(c,j) / dt
	end if	
		
	ccon_ch4s_sat(c,j) = ccon_ch4s_sat(c,j) - AOMCH4Oxid * dt
	ccon_ch4s_sat(c,j) = max(0._r8, ccon_ch4s_sat(c,j))
	
!	// For Microbe dyndamics
	AceMethanogenGrowth = m_dYAceMethanogens * AceCons * pHeffect
	AceMethanogenDying = m_dDeadRAceMethanogens * cacebios_sat(c,j) * pHeffect
	AceMethanogenGrowth = max(0._r8, AceMethanogenGrowth)
	AceMethanogenDying = max(0._r8, AceMethanogenDying)

	!print *, AceMethanogenGrowth," ",AceMethanogenDying
	H2MethanogenGrowth 			= m_dYH2Methanogens * 4. * H2CH4Prod * pHeffect
	H2MethanogenDying 			= m_dDeadRH2Methanogens * cco2bios_sat(c,j) * pHeffect
	H2MethanogenGrowth 			= max(0._r8, H2MethanogenGrowth)
	H2MethanogenDying 			= max(0._r8, H2MethanogenDying)

	MethanotrophGrowth 			= m_dYMethanotrophs * CH4Oxid * pHeffect
	MethanotrophDying 			= m_dDeadRMethanotrophs * caerch4bios_sat(c,j) * pHeffect
	MethanotrophGrowth 			= max(0._r8, MethanotrophGrowth)
	MethanotrophDying 			= max(0._r8, MethanotrophDying)

	AOMMethanotrophGrowth 		= m_dYAOMMethanotrophs * AOMCH4Oxid * pHeffect
	AOMMethanotrophDying 		= m_dDeadRAOMMethanotrophs * canaerch4bios_sat(c,j) * pHeffect
	AOMMethanotrophGrowth 		= max(0._r8, AOMMethanotrophGrowth)
	AOMMethanotrophDying 		= max(0._r8, AOMMethanotrophDying)

	!~ cacebios_sat(c,j) 			= cacebios_sat(c,j) + min((AceMethanogenGrowth - AceMethanogenDying), 0._r8) 
	!~ cco2bios_sat(c,j) 			= cco2bios_sat(c,j) + min((H2MethanogenGrowth - H2MethanogenDying), 0._r8)
	!~ caerch4bios_sat(c,j) 		= caerch4bios_sat(c,j) + min((MethanotrophGrowth - MethanotrophDying), 0._r8)
	!~ canaerch4bios_sat(c,j) 		= canaerch4bios_sat(c,j) + min((AOMMethanotrophGrowth - AOMMethanotrophDying), 0._r8) 

	cacebios_sat(c,j) 			= cacebios_sat(c,j) + (AceMethanogenGrowth - AceMethanogenDying) * dt
	cco2bios_sat(c,j) 			= cco2bios_sat(c,j) + (H2MethanogenGrowth - H2MethanogenDying) * dt
	caerch4bios_sat(c,j) 			= caerch4bios_sat(c,j) + (MethanotrophGrowth - MethanotrophDying) * dt
	canaerch4bios_sat(c,j) 		= canaerch4bios_sat(c,j) + (AOMMethanotrophGrowth - AOMMethanotrophDying) * dt
	
	cacebios_sat(c,j) 			= max(cacebios_sat(c,j),MFGbiomin)
	cco2bios_sat(c,j) 			= max(cco2bios_sat(c,j),MFGbiomin)
	caerch4bios_sat(c,j) 			= max(caerch4bios_sat(c,j),MFGbiomin)
	canaerch4bios_sat(c,j) 		= max(canaerch4bios_sat(c,j),MFGbiomin)
	
!write(iulog,*) "xiaofeng here3 ", c,j,cacebios_sat(c,j),cco2bios_sat(c,j),caerch4bios_sat(c,j),canaerch4bios_sat(c,j)	
!write(iulog,*) "grow, dead: ",	AceMethanogenGrowth," ", AceMethanogenDying," ", H2MethanogenGrowth," ", H2MethanogenDying
!write(iulog,*) "grow, dead: ",	MethanotrophGrowth," ", MethanotrophDying," ", AOMMethanotrophGrowth," ", AOMMethanotrophDying

	cdocs_sat(c,j) 						= ACConcentration
	caces_unsat_prod(c,j)					= AceProd
	caces_unsat_prod_h2(c,j)				= H2AceProd	
	ch4_prod_ace_depth_sat(c,j) 			= CH4Prod 
	ch4_prod_co2_depth_sat(c,j) 			= H2CH4Prod
	ch4_oxid_o2_depth_sat(c,j) 			= CH4Oxid
	ch4_oxid_aom_depth_sat(c,j) 			= AOMCH4Oxid
	ch4_aere_depth_sat(c,j) 				= CH4PlantFlux
	ch4_dif_depth_sat(c,j) 				= 0._r8
	ch4_ebul_depth_sat(c,j) 				= CH4Ebull

	co2_prod_ace_depth_sat(c,j) 			= ACCO2Prod
!	co2_decomp_depth_sat(c,j) 			! calcuated in first part of this subroutine
	co2_cons_depth_sat(c,j) 				= 0._r8
	co2_ebul_depth_sat(c,j) 				= 0._r8
	co2_aere_depth_sat(c,j) 				= 0._r8
	co2_dif_depth_sat(c,j) 				= 0._r8

	o2_cons_depth_sat(c,j) 				= AerO2Cons + CH4O2Cons
	o2_aere_depth_sat(c,j) 				= O2PlantFlux
	o2_aere_oxid_depth_sat(c,j) 			= PlantO2Cons
	o2_decomp_depth_sat(c,j) 			= co2_decomp_depth_sat (c,j)		! mole / m3
	o2_dif_depth_sat(c,j) 				= 0._r8

	h2_prod_depth_sat(c,j) 				= ACH2Prod
	h2_cons_depth_sat(c,j) 				= H2Cons
	h2_aere_depth_sat(c,j) 				= H2PlantFlux
	h2_diff_depth_sat(c,j) 				= 0.0
	h2_ebul_depth_sat(c,j) 				= 0.0
       end do
! average to grid level fluxes or state variables    
	tem1 = Henry_kHpc_w(1) * exp (Henry_C_w(1) * (1/soiltemp(c,1) - 1/kh_tbase))
	tem2 = Henry_kHpc_w(2) * exp (Henry_C_w(2) * (1/soiltemp(c,1) - 1/kh_tbase))
	tem3 = Henry_kHpc_w(3) * exp (Henry_C_w(3) * (1/soiltemp(c,1) - 1/kh_tbase))
	tem4 = Henry_kHpc_w(4) * exp (Henry_C_w(4) * (1/soiltemp(c,1) - 1/kh_tbase))
	
	tem1 = 1. / tem1 * atmch4 * 1000. / 16.
	tem2 = 1. / tem2 * atmo2 * 1000. / 32.
	tem3 = 1. / tem3 * atmco2 * 1000. / 44.
	tem4 = 1. / tem4 * atmh2 * 1000. / 2.
	
	ch4_surf_aere_sat(c)				= 0._r8
	ch4_surf_ebul_sat(c) 			= 0._r8
	ch4_surf_dif_sat(c) 				= 0._r8
	ch4_surf_netflux_sat(c)			= 0._r8
	
	co2_surf_aere_sat(c)				= 0._r8
	co2_surf_ebul_sat(c) 			= 0._r8
	co2_surf_dif_sat(c) 				= 0._r8
	co2_surf_netflux_sat(c)			= 0._r8
	
	o2_surf_aere_sat(c)				= 0._r8
	o2_surf_dif_sat(c) 				= 0._r8
	o2_surf_netflux_sat(c)			= 0._r8
	
	h2_surf_aere_sat(c)				= 0._r8
	h2_surf_ebul_sat(c) 				= 0._r8
	h2_surf_dif_sat(c) 				= 0._r8
	h2_surf_netflux_sat(c)			= 0._r8
	
if(soiltemp(c,1) < SHR_CONST_TKFRZ) then
!	ch4_surf_aere_sat(c)				= 0._r8
	ch4_surf_ebul_sat(c) 			= 0._r8
	ch4_surf_dif_sat(c) 				= 0._r8
		
!	co2_surf_aere_sat(c)				= 0._r8
	co2_surf_ebul_sat(c) 			= 0._r8
	co2_surf_dif_sat(c) 				= 0._r8
		
!	o2_surf_aere_sat(c)				= 0._r8
	o2_surf_dif_sat(c) 				= 0._r8
		
!	h2_surf_aere_sat(c)				= 0._r8
	h2_surf_ebul_sat(c) 				= 0._r8
	h2_surf_dif_sat(c) 				= 0._r8
else
	ch4_surf_dif_sat(c) = (ccon_ch4s_sat(c,1) - tem1) / dt !* (2.0 ** ((soiltemp(c,j) - 286.65) / 10.)) 		!max(0._r8, (ccon_ch4s_sat(c,1) - tem1))  ! solubility of ch4 is 0.0000227g/L
	ccon_ch4s_sat(c,1) = ccon_ch4s_sat(c,1) - ch4_surf_dif_sat(c) * dt
	ch4_surf_dif_sat(c) = ch4_surf_dif_sat(c) * dz(c,1)

	o2_surf_dif_sat(c) = (ccon_o2s_sat(c,1) - tem2) / dt !* (2.0 ** ((soiltemp(c,j) - 286.65) / 10.)) 			!max(0._r8, (ccon_o2s_sat(c,1) - tem4))  ! solubility of ch4 is 0.00004/L
	ccon_o2s_sat(c,1) = ccon_o2s_sat(c,1) - o2_surf_dif_sat(c) * dt
	o2_surf_dif_sat(c) = o2_surf_dif_sat(c) * dz(c,1)

	co2_surf_dif_sat(c) = (ccon_co2s_sat(c,1) - tem3) / dt !* (2.0 ** ((soiltemp(c,j) - 286.65) / 10.)) 		!max(0._r8, (ccon_co2s_sat(c,1) - tem3))  ! solubility of ch4 is 0.002/L
	ccon_co2s_sat(c,1) = ccon_co2s_sat(c,1) - co2_surf_dif_sat(c) * dt
	co2_surf_dif_sat(c) = co2_surf_dif_sat(c) * dz(c,1)

	h2_surf_dif_sat(c) = (ccon_h2s_sat(c,1) - tem4) / dt !* (2.0 ** ((soiltemp(c,j) - 286.65) / 10.))			!max(0._r8, (ccon_h2s_sat(c,1) - tem2))  ! solubility of ch4 is 0.0000015/L
	ccon_h2s_sat(c,1) = ccon_h2s_sat(c,1) - h2_surf_dif_sat(c) * dt
	h2_surf_dif_sat(c) = h2_surf_dif_sat(c) * dz(c,1)

	do j = 1,nlevsoi
	ch4_surf_aere_sat(c)				= ch4_surf_aere_sat(c) + ch4_aere_depth_sat(c,j) * dz(c,j)
!	ch4_surf_ebul_sat(c) 			= ch4_surf_ebul_sat(c) + ch4_ebul_depth_sat(c,j)  * dz(c,j)
		
	co2_surf_aere_sat(c)				= co2_surf_aere_sat(c) + co2_aere_depth_sat(c,j) * dz(c,j)
!	co2_surf_ebul_sat(c) 			= co2_surf_ebul_sat(c) + co2_ebul_depth_sat(c,j) * dz(c,j)
		
	o2_surf_aere_sat(c)				= o2_surf_aere_sat(c) + o2_aere_depth_sat(c,j) * dz(c,j)
	
	h2_surf_aere_sat(c)				= h2_surf_aere_sat(c) + h2_aere_depth_sat(c,j) * dz(c,j)
!	h2_surf_ebul_sat(c) 				= h2_surf_ebul_sat(c) + h2_ebul_depth_sat(c,j) * dz(c,j)
	end do

! gas move up
	l = clandunit(c)     
	if (ltype(l) == istsoil .or. ltype(l) == istcrop) then 
	do j = 3,nlevsoi
	if(j > jwaterhead_unsat(c) .and. j < nlevsoi) then
		ch4_ebul_depth_unsat(c,j) = max(0._r8, ch4_ebul_depth_unsat(c,j))
!		o2_ebul_depth_unsat(c,j) = max(0._r8, o2_ebul_depth_unsat(c,j))
		co2_ebul_depth_unsat(c,j) = max(0._r8, co2_ebul_depth_unsat(c,j))
		h2_ebul_depth_unsat(c,j) = max(0._r8, h2_ebul_depth_unsat(c,j))
		
		ccon_ch4s_unsat(c,j-1) = (ccon_ch4s_unsat(c,j-1) * dz(c,j-1) + ch4_ebul_depth_unsat(c,j)) / dz(c,j-1)
		ccon_ch4s_unsat(c,j) = (ccon_ch4s_unsat(c,j) * dz(c,j) - ch4_ebul_depth_unsat(c,j)) / dz(c,j)
		
!		ccon_o2s_unsat(c,j-1) = (ccon_o2s_unsat(c,j-1) * dz(c,j-1) + o2_ebul_depth_unsat(c,j)) / dz(c,j-1)
!		ccon_o2s_unsat(c,j) = (ccon_o2s_unsat(c,j) * dz(c,j) - o2_ebul_depth_unsat(c,j)) / dz(c,j)

		ccon_co2s_unsat(c,j-1) = (ccon_co2s_unsat(c,j-1) * dz(c,j-1) + co2_ebul_depth_unsat(c,j)) / dz(c,j-1)
		ccon_co2s_unsat(c,j) = (ccon_co2s_unsat(c,j) * dz(c,j) - co2_ebul_depth_unsat(c,j)) / dz(c,j)

		ccon_h2s_unsat(c,j-1) = (ccon_h2s_unsat(c,j-1) * dz(c,j-1) + h2_ebul_depth_unsat(c,j)) / dz(c,j-1)
		ccon_h2s_unsat(c,j) = (ccon_h2s_unsat(c,j) * dz(c,j) - h2_ebul_depth_unsat(c,j)) / dz(c,j)
	end if
	! for the saturation portion
		ch4_ebul_depth_sat(c,j) = max(0._r8, ch4_ebul_depth_sat(c,j))
!		o2_ebul_depth_sat(c,j) = max(0._r8, o2_ebul_depth_sat(c,j))
		co2_ebul_depth_sat(c,j) = max(0._r8, co2_ebul_depth_sat(c,j))
		h2_ebul_depth_sat(c,j) = max(0._r8, h2_ebul_depth_sat(c,j))
		
		ccon_ch4s_sat(c,j-1) = (ccon_ch4s_sat(c,j-1) * dz(c,j-1) + ch4_ebul_depth_sat(c,j) * dz(c,j)) / dz(c,j-1)
		ccon_ch4s_sat(c,j) = (ccon_ch4s_sat(c,j) * dz(c,j) - ch4_ebul_depth_sat(c,j) * dz(c,j)) / dz(c,j)
		
!		ccon_o2s_sat(c,j-1) = (ccon_o2s_sat(c,j-1) * dz(c,j-1) + o2_ebul_depth_sat(c,j)) / dz(c,j-1)
!		ccon_o2s_sat(c,j) = (ccon_o2s_sat(c,j) * dz(c,j) - o2_ebul_depth_sat(c,j)) / dz(c,j)

		ccon_co2s_sat(c,j-1) = (ccon_co2s_sat(c,j-1) * dz(c,j-1) + co2_ebul_depth_sat(c,j)* dz(c,j)) / dz(c,j-1)
		ccon_co2s_sat(c,j) = (ccon_co2s_sat(c,j) * dz(c,j) - co2_ebul_depth_sat(c,j)* dz(c,j)) / dz(c,j)

		ccon_h2s_sat(c,j-1) = (ccon_h2s_sat(c,j-1) * dz(c,j-1) + h2_ebul_depth_sat(c,j)* dz(c,j)) / dz(c,j-1)
		ccon_h2s_sat(c,j) = (ccon_h2s_sat(c,j) * dz(c,j) - h2_ebul_depth_sat(c,j)* dz(c,j)) / dz(c,j)
!	write(iulog,*) "after ", ccon_ch4s_sat(c,j-1), ccon_ch4s_sat(c,j)
	end do
        end if
! end gas move up
	do j = 1,nlevsoi
	ch4_surf_ebul_sat(c) 			= ch4_surf_ebul_sat(c) + ch4_ebul_depth_sat(c,j) * dz(c,j) !+ ch4_ebul_depth_sat(c,2) * dz(c,2)
!	o2_surf_ebul_sat(c) 				= o2_surf_ebul_sat(c) + o2_ebul_depth_sat(c,1) * dz(c,1) + o2_ebul_depth_sat(c,2) * dz(c,2)
	co2_surf_ebul_sat(c) 			= co2_surf_ebul_sat(c) + co2_ebul_depth_sat(c,j) * dz(c,j) !+ co2_ebul_depth_sat(c,2) * dz(c,2)
	h2_surf_ebul_sat(c) 				= h2_surf_ebul_sat(c) + h2_ebul_depth_sat(c,j) * dz(c,j) !+ h2_ebul_depth_sat(c,2) * dz(c,2)
	end do
!	ch4_surf_ebul_sat(c) 			= ch4_surf_ebul_sat(c) + ch4_ebul_depth_sat(c,1) * dz(c,1) + ch4_ebul_depth_sat(c,2) * dz(c,2)
!!	o2_surf_ebul_sat(c) 				= o2_surf_ebul_sat(c) + o2_ebul_depth_sat(c,1) * dz(c,1) + o2_ebul_depth_sat(c,2) * dz(c,2)
!	co2_surf_ebul_sat(c) 			= co2_surf_ebul_sat(c) + co2_ebul_depth_sat(c,1) * dz(c,1) + co2_ebul_depth_sat(c,2) * dz(c,2)
!	h2_surf_ebul_sat(c) 				= h2_surf_ebul_sat(c) + h2_ebul_depth_sat(c,1) * dz(c,1) + h2_ebul_depth_sat(c,2) * dz(c,2)
	
end if

	ch4_surf_netflux_sat(c) 			= ch4_surf_netflux_sat(c) + ch4_surf_dif_sat(c) + ch4_surf_aere_sat(c) + ch4_surf_ebul_sat(c)
	co2_surf_netflux_sat(c) 			= co2_surf_netflux_sat(c) + co2_surf_dif_sat(c) + co2_surf_aere_sat(c) + co2_surf_ebul_sat(c)
	o2_surf_netflux_sat(c) 			= o2_surf_netflux_sat(c) + o2_surf_dif_sat(c) + o2_surf_aere_sat(c) !+ o2_surf_ebul_sat(c)
	h2_surf_netflux_sat(c) 			= h2_surf_netflux_sat(c) + h2_surf_dif_sat(c) + h2_surf_aere_sat(c) + h2_surf_ebul_sat(c)
! average to grid level fluxes or state variables    
end do
! end of the saturated fraction simulations

	do fc = 1,num_soilc
        c = filter_soilc(fc)
	l = clandunit(c)     
	if (ltype(l) == istsoil .or. ltype(l) == istcrop) then 
	do j = 1,nlevsoi
	cdocs_sat(c,j) 			= cdocs_sat(c,j) * 12. 		!/ dz(c,j)
	caces_sat(c,j) 			= caces_sat(c,j) * 12. 		!/ dz(c,j)
	cacebios_sat(c,j) 		= cacebios_sat(c,j) * 12. 		!/ dz(c,j)
	cco2bios_sat(c,j) 		= cco2bios_sat(c,j) * 12. 		!/ dz(c,j)
	caerch4bios_sat(c,j) 		= caerch4bios_sat(c,j) * 12. 	!/ dz(c,j)
	canaerch4bios_sat(c,j) 	= canaerch4bios_sat(c,j) * 12. 	!/ dz(c,j)
	ccon_o2s_sat(c,j) 		= ccon_o2s_sat(c,j) * 32. 		!/ dz(c,j)
	ccon_ch4s_sat(c,j) 		= ccon_ch4s_sat(c,j) * 12. 		!/ dz(c,j)
	ccon_h2s_sat(c,j) 		= ccon_h2s_sat(c,j) * 2. 		!/ dz(c,j)
	ccon_co2s_sat(c,j) 		= ccon_co2s_sat(c,j) * 12. 		!/ dz(c,j)
	end do
	end if
	end do
	
! vertical diffusion of ACE 
	do fc = 1,num_soilc
        c = filter_soilc(fc)
	l = clandunit(c)     
	if (ltype(l) == istsoil .or. ltype(l) == istcrop) then 
	do j = 2,nlevsoi
	if(j > jwaterhead_unsat(c) .and. j <= nlevsoi) then
		caces_unsat_temp(c,j) = (caces_unsat(c,j-1) - caces_unsat(c,j)) * dom_diffus * (soiltemp(c,j)/298)**1.87 * get_step_size() / (z(c,j) - z(c,j-1)) !CH4_dif

		if(abs(caces_unsat_temp(c,j)) < abs(caces_unsat(c,j-1) - (caces_unsat(c,j-1) * dz(c,j-1) + caces_unsat(c,j) * dz(c,j)) / (dz(c,j) + dz(c,j-1))))  then
		caces_unsat_temp(c,j) = caces_unsat_temp(c,j)
		else
		caces_unsat_temp(c,j) = caces_unsat(c,j-1) - (caces_unsat(c,j-1) * dz(c,j-1) + caces_unsat(c,j) * dz(c,j)) / (dz(c,j) + dz(c,j-1))
		end if
		!write(*,*)"before: unsat ",j, caces_unsat(c,j-1),caces_unsat(c,j)
		caces_unsat(c,j-1) = caces_unsat(c,j-1) - caces_unsat_temp(c,j) * dz(c,j-1)
		caces_unsat(c,j) = caces_unsat(c,j) + caces_unsat_temp(c,j) * dz(c,j)	
		!write(*,*)"after: unsat ",j, caces_unsat(c,j-1),caces_unsat(c,j)
	end if
		caces_sat_temp(c,j) = (caces_sat(c,j-1) - caces_sat(c,j)) * dom_diffus * (soiltemp(c,j)/298)**1.87 * get_step_size() / (z(c,j) - z(c,j-1)) !CH4_dif

		if(abs(caces_sat_temp(c,j)) < abs(caces_sat(c,j-1) - (caces_sat(c,j-1) * dz(c,j-1) + caces_sat(c,j) * dz(c,j)) / (dz(c,j) + dz(c,j-1))))  then
		caces_sat_temp(c,j) = caces_sat_temp(c,j)
		else
		caces_sat_temp(c,j) = caces_sat(c,j-1) - (caces_sat(c,j-1) * dz(c,j-1) + caces_sat(c,j) * dz(c,j)) / (dz(c,j) + dz(c,j-1))
		end if
		!write(*,*)"befre: sat ",j, caces_sat(c,j-1),caces_sat(c,j)
		caces_sat(c,j-1) = caces_sat(c,j-1) - caces_sat_temp(c,j) * dz(c,j-1)
		caces_sat(c,j) = caces_sat(c,j) + caces_sat_temp(c,j) * dz(c,j)
		!write(*,*)"after: sat ",j, caces_sat(c,j-1),caces_sat(c,j)
	end do
        end if
	end do
! vertical diffusion of ACE

! keeping ACE consistent with DOC    
	do fc = 1,num_soilc
        c = filter_soilc(fc)
	l = clandunit(c)     
	if (ltype(l) == istsoil .or. ltype(l) == istcrop) then 
	do j = 2,nlevsoi
	if(j > jwaterhead_unsat(c) .and. j <= nlevsoi) then
		cdocs_unsat_temp(c,j) = (cdocs_unsat(c,j-1) - cdocs_unsat(c,j)) * dom_diffus * (soiltemp(c,j)/298)**1.87 * get_step_size() / (z(c,j) - z(c,j-1)) !CH4_dif

		if(abs(cdocs_unsat_temp(c,j)) < abs(cdocs_unsat(c,j-1) - (cdocs_unsat(c,j-1) * dz(c,j-1) + cdocs_unsat(c,j) * dz(c,j)) / (dz(c,j) + dz(c,j-1))))  then
		cdocs_unsat_temp(c,j) = cdocs_unsat_temp(c,j)
		else
		cdocs_unsat_temp(c,j) = cdocs_unsat(c,j-1) - (cdocs_unsat(c,j-1) * dz(c,j-1) + cdocs_unsat(c,j) * dz(c,j)) / (dz(c,j) + dz(c,j-1))
		end if
		!write(*,*)"before: unsat ",j, cdocs_unsat(c,j-1),cdocs_unsat(c,j)
		cdocs_unsat(c,j-1) = cdocs_unsat(c,j-1) - cdocs_unsat_temp(c,j) * dz(c,j-1)
		cdocs_unsat(c,j) = cdocs_unsat(c,j) + cdocs_unsat_temp(c,j) * dz(c,j)	
		!write(*,*)"after: unsat ",j, cdocs_unsat(c,j-1),cdocs_unsat(c,j)
	end if
		cdocs_sat_temp(c,j) = (cdocs_sat(c,j-1) - cdocs_sat(c,j)) * dom_diffus * (soiltemp(c,j)/298)**1.87 * get_step_size() / (z(c,j) - z(c,j-1)) !CH4_dif

		if(abs(cdocs_sat_temp(c,j)) < abs(cdocs_sat(c,j-1) - (cdocs_sat(c,j-1) * dz(c,j-1) + cdocs_sat(c,j) * dz(c,j)) / (dz(c,j) + dz(c,j-1))))  then
		cdocs_sat_temp(c,j) = cdocs_sat_temp(c,j)
		else
		cdocs_sat_temp(c,j) = cdocs_sat(c,j-1) - (cdocs_sat(c,j-1) * dz(c,j-1) + cdocs_sat(c,j) * dz(c,j)) / (dz(c,j) + dz(c,j-1))
		end if
		!write(*,*)"befre: sat ",j, cdocs_sat(c,j-1),cdocs_sat(c,j)
		cdocs_sat(c,j-1) = cdocs_sat(c,j-1) - cdocs_sat_temp(c,j) * dz(c,j-1)
		cdocs_sat(c,j) = cdocs_sat(c,j) + cdocs_sat_temp(c,j) * dz(c,j)
		!write(*,*)"after: sat ",j, cdocs_sat(c,j-1),cdocs_sat(c,j)
	end do
        end if
	end do
! end of keeping ACE pattern consistent with DOC

!#if (defined HUM_HOL)
!   real(r8), public :: som_diffus = 1e-4_r8 / (secspday * 365._r8)  ! [m^2/sec] = 1 cm^2 / yr
!   real(r8), public :: dom_diffus = 10._r8				! times to som diffus 
!#endif

#if (defined HUM_HOL)  
	do c=1,2
        jwt(c) = nlevsoi
	! allow jwt to equal zero when zwt is in top layer
	do j = 1,nlevsoi
          if(zwt(c) <= zi(c,j)) then
             jwt(c) = j - 1
             exit
          end if
	enddo
	enddo
       
	! beginning of lateral transport of BGC variables in association with lateral water flow
        !do fc = 1, 2
	 !   c = fc
	 
	!~ do j = 1,nlevsoi
	!~ if(j >= jwt(1) .and. j <= (jwt(1)+1)) then
		!~ hum_frac = 0.75_r8
		!~ hol_frac = 0.25_r8
		!~ lat_flux_factor1 		= (hol_frac/hum_frac)
		!~ lat_flux_factor2 		= (hum_frac/hol_frac)
	!~ !if(jwt(c) < nlevsoi) then  ! does not consider the scenario of water table below soil column
!~ !	write(iulog,*) "LBGC debug 5", cdocs_unsat(1,j), " ", qflx_lat_aqu(1), " ", cdocs_unsat(2,j), " ", qflx_lat_aqu(2)
	!~ if(wa(1) > 1e-10 .and. abs(qflx_lat_aqu(1)) > 1e-10) then
		!~ lxdomunsat 		= qflx_lat_aqu(1) / sqrt(lat_flux_factor1) / vwc(1,j) * cdocs_unsat(1,j)! * get_step_size()
		!~ lxdomsat 			= qflx_lat_aqu(1) / sqrt(lat_flux_factor1) / vwc(1,j) * cdocs_sat(1,j)! * get_step_size()
		!~ lxaceunsat 			= qflx_lat_aqu(1) / sqrt(lat_flux_factor1) / vwc(1,j) * cdocs_unsat(1,j)! * get_step_size()
		!~ lxacesat 			= qflx_lat_aqu(1) / sqrt(lat_flux_factor1) / vwc(1,j) * caces_sat(1,j)! * get_step_size()
	!~ else
		!~ lxdomunsat 		= 0.
		!~ lxdomsat 			= 0.
		!~ lxaceunsat 			= 0.
		!~ lxacesat 			= 0.
	!~ end if
		!~ if(abs(lxdomunsat) < abs(cdocs_unsat(1,j) - (cdocs_unsat(1,j) * sqrt(lat_flux_factor2) + cdocs_unsat(2,j) * sqrt(lat_flux_factor1))) / (sqrt(lat_flux_factor2) + sqrt(lat_flux_factor1)))  then
		!~ lxdomunsat = lxdomunsat
		!~ else
		!~ lxdomunsat = cdocs_unsat(1,j) - ((cdocs_unsat(1,j) * sqrt(lat_flux_factor2) + cdocs_unsat(2,j) * sqrt(lat_flux_factor1))) / (sqrt(lat_flux_factor2) + sqrt(lat_flux_factor1))
		!~ end if
		!~ cdocs_unsat(1,j) 	= cdocs_unsat(1,j) - lxdomunsat * sqrt(lat_flux_factor1)
		!~ cdocs_unsat(2,j) 	= cdocs_unsat(2,j) + lxdomunsat * sqrt(lat_flux_factor2) 
		
		!~ if(abs(lxdomsat) < abs(cdocs_sat(1,j) - (cdocs_sat(1,j) * sqrt(lat_flux_factor2) + cdocs_sat(2,j) * sqrt(lat_flux_factor1))) / (sqrt(lat_flux_factor2) + sqrt(lat_flux_factor1)))  then
		!~ lxdomsat = lxdomsat
		!~ else
		!~ lxdomsat = cdocs_sat(1,j) - ((cdocs_sat(1,j) * sqrt(lat_flux_factor2) + cdocs_sat(2,j) * sqrt(lat_flux_factor1))) / (sqrt(lat_flux_factor2) + sqrt(lat_flux_factor1))
		!~ end if		
		!~ cdocs_sat(1,j) 		= cdocs_sat(1,j) - lxdomsat * sqrt(lat_flux_factor1)
		!~ cdocs_sat(2,j) 		= cdocs_sat(2,j) + lxdomsat * sqrt(lat_flux_factor2)

		!~ if(abs(lxaceunsat) < abs(caces_unsat(1,j) - (caces_unsat(1,j) * sqrt(lat_flux_factor2) + caces_unsat(2,j) * sqrt(lat_flux_factor1))) / (sqrt(lat_flux_factor2) + sqrt(lat_flux_factor1)))  then
		!~ lxaceunsat = lxaceunsat
		!~ else
		!~ lxaceunsat = caces_unsat(1,j) - ((caces_unsat(1,j) * sqrt(lat_flux_factor2) + caces_unsat(2,j) * sqrt(lat_flux_factor1))) / (sqrt(lat_flux_factor2) + sqrt(lat_flux_factor1))
		!~ end if
		!~ caces_unsat(1,j) 	= caces_unsat(1,j) - lxaceunsat * sqrt(lat_flux_factor1)
		!~ caces_unsat(2,j) 	= caces_unsat(2,j) + lxaceunsat * sqrt(lat_flux_factor2)

		!~ if(abs(lxacesat) < abs(caces_sat(1,j) - (caces_sat(1,j) * sqrt(lat_flux_factor2) + caces_sat(2,j) * sqrt(lat_flux_factor1))) / (sqrt(lat_flux_factor2) + sqrt(lat_flux_factor1)))  then
		!~ lxacesat = lxacesat
		!~ else
		!~ lxacesat = caces_sat(1,j) - ((caces_sat(1,j) * sqrt(lat_flux_factor2) + caces_sat(2,j) * sqrt(lat_flux_factor1))) / (sqrt(lat_flux_factor2) + sqrt(lat_flux_factor1))
		!~ end if
		!~ caces_sat(1,j) 		= caces_sat(1,j) - lxacesat * sqrt(lat_flux_factor1)
		!~ caces_sat(2,j) 		= caces_sat(2,j) + lxacesat * sqrt(lat_flux_factor2)
	!~ !write(iulog,*) "LBGC debug 6", cdocs_unsat(1,j), " ", qflx_lat_aqu(1), " ", cdocs_unsat(2,j), " ", qflx_lat_aqu(2)
		!~ !endif
	!~ end if
	!~ end do
	
	!end do
	! end of lateral transport of BGC variables in association with lateral water flow

	! beginning of lateral diffusion of BGC variables
	!~ do j = 1,nlevsoi
	!~ if(jwt(1)<=j) then
		!~ hum_frac 			= 0.75_r8
		!~ hol_frac 			= 0.25_r8
		!~ lat_flux_factor1 	= (hol_frac/hum_frac)
		!~ lat_flux_factor2 	= (hum_frac/hol_frac)
		
!~ !	if(wa(1) > 1e-10 .and. wa(2) > 1e-10) then
		!~ lxdomunsat 		= (cdocs_unsat(1,j) - cdocs_unsat(2,j)) * dom_diffus * get_step_size()
		!~ lxdomsat 			= (cdocs_sat(1,j) - cdocs_sat(2,j)) * dom_diffus * get_step_size()
		!~ lxaceunsat 		= (caces_unsat(1,j) - caces_unsat(2,j)) * dom_diffus *100* get_step_size()
		!~ lxacesat 			= (caces_sat(1,j) - caces_sat(2,j)) * dom_diffus *100* get_step_size()
!~ !write(iulog,*) "LBGC debug 7", lxdomunsat, cdocs_unsat(1,j), cdocs_unsat(2,j)
		!~ if(abs(lxdomunsat) < abs(cdocs_unsat(1,j) - (cdocs_unsat(1,j) * sqrt(lat_flux_factor2) + cdocs_unsat(2,j) * sqrt(lat_flux_factor1))) / (sqrt(lat_flux_factor2) + sqrt(lat_flux_factor1)))  then
		!~ lxdomunsat = lxdomunsat
		!~ else
		!~ lxdomunsat = cdocs_unsat(1,j) - ((cdocs_unsat(1,j) * sqrt(lat_flux_factor2) + cdocs_unsat(2,j) * sqrt(lat_flux_factor1))) / (sqrt(lat_flux_factor2) + sqrt(lat_flux_factor1))
		!~ end if
		!~ cdocs_unsat(1,j) 	= cdocs_unsat(1,j) - lxdomunsat * sqrt(lat_flux_factor1)
		!~ cdocs_unsat(2,j) 	= cdocs_unsat(2,j) + lxdomunsat * sqrt(lat_flux_factor2)

		!~ if(abs(lxdomsat) < abs(cdocs_sat(1,j) - (cdocs_sat(1,j) * sqrt(lat_flux_factor2) + cdocs_sat(2,j) * sqrt(lat_flux_factor1))) / (sqrt(lat_flux_factor2) + sqrt(lat_flux_factor1)))  then
		!~ lxdomsat = lxdomsat
		!~ else
		!~ lxdomsat = cdocs_sat(1,j) - ((cdocs_sat(1,j) * sqrt(lat_flux_factor2) + cdocs_sat(2,j) * sqrt(lat_flux_factor1))) / (sqrt(lat_flux_factor2) + sqrt(lat_flux_factor1))
		!~ end if
		!~ cdocs_sat(1,j) 		= cdocs_sat(1,j) - lxdomsat * sqrt(lat_flux_factor1)
		!~ cdocs_sat(2,j) 		= cdocs_sat(2,j) + lxdomsat * sqrt(lat_flux_factor2)

		!~ if(abs(lxaceunsat) < abs(caces_unsat(1,j) - (caces_unsat(1,j) * sqrt(lat_flux_factor2) + caces_unsat(2,j) * sqrt(lat_flux_factor1))) / (sqrt(lat_flux_factor2) + sqrt(lat_flux_factor1)))  then
		!~ lxaceunsat = lxaceunsat
		!~ else
		!~ lxaceunsat = caces_unsat(1,j) - ((caces_unsat(1,j) * sqrt(lat_flux_factor2) + caces_unsat(2,j) * sqrt(lat_flux_factor1))) / (sqrt(lat_flux_factor2) + sqrt(lat_flux_factor1))
		!~ end if
		!~ caces_unsat(1,j) 	= caces_unsat(1,j) - lxaceunsat * sqrt(lat_flux_factor1)
		!~ caces_unsat(2,j) 	= caces_unsat(2,j) + lxaceunsat * sqrt(lat_flux_factor2)

		!~ if(abs(lxacesat) < abs(caces_sat(1,j) - (caces_sat(1,j) * sqrt(lat_flux_factor2) + caces_sat(2,j) * sqrt(lat_flux_factor1))) / (sqrt(lat_flux_factor2) + sqrt(lat_flux_factor1)))  then
		!~ lxacesat = lxacesat
		!~ else
		!~ lxacesat = caces_sat(1,j) - ((caces_sat(1,j) * sqrt(lat_flux_factor2) + caces_sat(2,j) * sqrt(lat_flux_factor1))) / (sqrt(lat_flux_factor2) + sqrt(lat_flux_factor1))
		!~ end if
		!~ caces_sat(1,j) 		= caces_sat(1,j) - lxacesat * sqrt(lat_flux_factor1)
		!~ caces_sat(2,j) 		= caces_sat(2,j)  + lxacesat * sqrt(lat_flux_factor2)

!~ !write(iulog,*) "LBGC debug 8", lxdomunsat, cdocs_unsat(1,j), cdocs_unsat(2,j)
	!~ end if
	!~ end do
	! end of lateral diffusion of BGC variables
	
#endif

      do fc = 1,num_soilc
        c = filter_soilc(fc)
      do j = 1,nlevsoi
      if(j > jwaterhead_unsat(c)) then
      micfinundated = 0.99
      else
      micfinundated = finundated(c)
#if (defined HUM_HOL) 
      micfinundated = 0.99
#endif
      end if      
	cdocs(c,j) 				= cdocs_unsat(c,j) * (1.0 - micfinundated) + cdocs_sat(c,j) * micfinundated
	caces(c,j) 				= caces_unsat(c,j) * (1.0 - micfinundated)  + caces_sat(c,j) * micfinundated
	caces_prod(c,j) 		= caces_unsat_prod(c,j) * (1.0 - micfinundated)  + caces_sat_prod(c,j) * micfinundated
	caces_prod_h2(c,j) 		= caces_unsat_prod_h2(c,j) * (1.0 - micfinundated)  + caces_sat_prod_h2(c,j) * micfinundated
	cacebios(c,j) 			= cacebios_unsat(c,j) * (1.0 - micfinundated)  + cacebios_sat(c,j) * micfinundated
	cco2bios(c,j) 			= cco2bios_unsat(c,j) * (1.0 - micfinundated)  + cco2bios_sat(c,j) * micfinundated
	caerch4bios(c,j) 		= caerch4bios_unsat(c,j) * (1.0 - micfinundated)  + caerch4bios_sat(c,j) * micfinundated
	canaerch4bios(c,j) 		= canaerch4bios_unsat(c,j) * (1.0 - micfinundated)  + canaerch4bios_sat(c,j) * micfinundated 
	ccon_o2s(c,j) 			= ccon_o2s_unsat(c,j) * (1.0 - micfinundated)  + ccon_o2s_sat(c,j) * micfinundated
	ccon_ch4s(c,j) 			= ccon_ch4s_unsat(c,j) * (1.0 - micfinundated)  + ccon_ch4s_sat(c,j) * micfinundated
	ccon_h2s(c,j) 			= ccon_h2s_unsat(c,j) * (1.0 - micfinundated)  + ccon_h2s_sat(c,j) * micfinundated
	ccon_co2s(c,j)			= ccon_co2s_unsat(c,j) * (1.0 - micfinundated)  + ccon_co2s_sat(c,j) * micfinundated

	ch4_prod_ace_depth(c,j)	= ch4_prod_ace_depth_unsat(c,j) * (1.0 - micfinundated) + ch4_prod_ace_depth_sat(c,j) * micfinundated
	ch4_prod_co2_depth(c,j)	= ch4_prod_co2_depth_unsat(c,j) * (1.0 - micfinundated) + ch4_prod_co2_depth_sat(c,j) * micfinundated
	ch4_oxid_o2_depth(c,j)	= ch4_oxid_o2_depth_unsat(c,j) * (1.0 - micfinundated) + ch4_oxid_o2_depth_sat(c,j) * micfinundated
	ch4_oxid_aom_depth(c,j)	= ch4_oxid_aom_depth_unsat(c,j) * (1.0 - micfinundated) + ch4_oxid_aom_depth_sat(c,j) * micfinundated
	ch4_aere_depth(c,j)		= ch4_aere_depth_unsat(c,j) * (1.0 - micfinundated) + ch4_aere_depth_sat(c,j) * micfinundated
	ch4_dif_depth(c,j)		= ch4_dif_depth_unsat(c,j) * (1.0 - micfinundated) + ch4_dif_depth_sat(c,j) * micfinundated
	ch4_ebul_depth(c,j)		= ch4_ebul_depth_unsat(c,j) * (1.0 - micfinundated) + ch4_ebul_depth_sat(c,j) * micfinundated
	
	co2_prod_ace_depth(c,j)	= co2_prod_ace_depth_unsat(c,j) * (1.0 - micfinundated) + co2_prod_ace_depth_sat(c,j) * micfinundated
	co2_decomp_depth(c,j)	= co2_decomp_depth_unsat(c,j) * (1.0 - micfinundated) + co2_decomp_depth_sat(c,j) * micfinundated 
	co2_cons_depth(c,j)		= co2_cons_depth_unsat(c,j) * (1.0 - micfinundated) + co2_cons_depth_sat(c,j) * micfinundated
	co2_ebul_depth(c,j)		= co2_ebul_depth_unsat(c,j) * (1.0 - micfinundated) + co2_ebul_depth_sat(c,j) * micfinundated
	co2_aere_depth(c,j)		= co2_aere_depth_unsat(c,j) * (1.0 - micfinundated) + co2_aere_depth_sat(c,j) * micfinundated
	co2_dif_depth(c,j)		= co2_dif_depth_unsat(c,j) * (1.0 - micfinundated) + co2_dif_depth_sat(c,j) * micfinundated
	
	o2_cons_depth(c,j)		= o2_cons_depth_unsat(c,j) * (1.0 - micfinundated) + o2_cons_depth_sat(c,j) * micfinundated
	o2_aere_depth(c,j)		= o2_aere_depth_unsat(c,j) * (1.0 - micfinundated) + o2_aere_depth_sat(c,j) * micfinundated
	o2_aere_oxid_depth(c,j)	= o2_aere_oxid_depth_unsat(c,j) * (1.0 - micfinundated) + o2_aere_oxid_depth_sat(c,j) * micfinundated
	o2_decomp_depth(c,j)	= o2_decomp_depth_unsat(c,j) * (1.0 - micfinundated) + o2_decomp_depth_sat(c,j) * micfinundated
	o2_dif_depth(c,j)		= o2_dif_depth_unsat(c,j) * (1.0 - micfinundated) + o2_dif_depth_sat(c,j) * micfinundated
	
	h2_prod_depth(c,j)		= h2_prod_depth_unsat(c,j) * (1.0 - micfinundated) + h2_prod_depth_sat(c,j) * micfinundated
	h2_cons_depth(c,j)		= h2_cons_depth_unsat(c,j) * (1.0 - micfinundated) + h2_cons_depth_sat(c,j) * micfinundated
	h2_aere_depth(c,j)		= h2_aere_depth_unsat(c,j) * (1.0 - micfinundated) + h2_aere_depth_sat(c,j) * micfinundated
	h2_diff_depth(c,j)		= h2_diff_depth_unsat(c,j) * (1.0 - micfinundated) + h2_diff_depth_sat(c,j) * micfinundated
	h2_ebul_depth(c,j)		= h2_ebul_depth_unsat(c,j) * (1.0 - micfinundated) + h2_ebul_depth_sat(c,j) * micfinundated

!	cdons_min(c,j)			= (cdocs_pre(c,j) - cdocs(c,j)) / cn_dom check the code in ~30 lines below 
	
      end do
      
      ch4_surf_aere(c)			= ch4_surf_aere_unsat(c) * (1.0 - micfinundated) + ch4_surf_aere_sat(c) * micfinundated
      ch4_surf_ebul(c)			= ch4_surf_ebul_unsat(c) * (1.0 - micfinundated) + ch4_surf_ebul_sat(c) * micfinundated
      ch4_surf_dif(c)			= ch4_surf_dif_unsat(c) * (1.0 - micfinundated) + ch4_surf_dif_sat(c) * micfinundated
      ch4_surf_netflux(c)		= ch4_surf_netflux_unsat(c) * (1.0 - micfinundated) + ch4_surf_netflux_sat(c) * micfinundated
      co2_surf_aere(c)			= co2_surf_aere_unsat(c) * (1.0 - micfinundated) + co2_surf_aere_sat(c) * micfinundated
      co2_surf_ebul(c)			= co2_surf_ebul_unsat(c) * (1.0 - micfinundated) + co2_surf_ebul_sat(c) * micfinundated
      co2_surf_dif(c)			= co2_surf_dif_unsat(c) * (1.0 - micfinundated) + co2_surf_dif_sat(c) * micfinundated
      co2_surf_netflux(c)		= co2_surf_netflux_unsat(c) * (1.0 - micfinundated) + co2_surf_netflux_sat(c) * micfinundated
      o2_surf_aere(c)			= o2_surf_aere_unsat(c) * (1.0 - micfinundated) + o2_surf_aere_sat(c) * micfinundated
      o2_surf_dif(c)			= o2_surf_dif_unsat(c) * (1.0 - micfinundated) + o2_surf_dif_sat(c) * micfinundated
      o2_surf_netflux(c)		= o2_surf_netflux_unsat(c) * (1.0 - micfinundated) + o2_surf_netflux_sat(c) * micfinundated
      h2_surf_aere(c)			= h2_surf_aere_unsat(c) * (1.0 - micfinundated) + h2_surf_aere_sat(c) * micfinundated
      h2_surf_ebul(c)			= h2_surf_ebul_unsat(c) * (1.0 - micfinundated) + h2_surf_ebul_sat(c) * micfinundated
      h2_surf_dif(c)			= h2_surf_dif_unsat(c) * (1.0 - micfinundated) + h2_surf_dif_sat(c) * micfinundated
      h2_surf_netflux(c)		= h2_surf_netflux_unsat(c) * (1.0 - micfinundated) + h2_surf_netflux_sat(c) * micfinundated
      end do

! added by xiaofeng on July 21 
       do j = 1,nlevdecomp
            do fc = 1,num_soilc
               c = filter_soilc(fc)
!	      write(iulog,*) "microbial mod: ",decomp_cpools_vr(c,j,i_dom), cdocs(c,j), decomp_npools_vr(c,j,i_dom), cdons(c,j)
!      cdons_min(c,j) = cdons_min(c,j) + decomp_npools_vr(c,j,i_dom) - (cdocs(c,j)) / cn_dom
      decomp_cpools_vr(c,j,i_dom) = (cdocs(c,j)) ! gC/m3
      decomp_npools_vr(c,j,i_dom) = (cdocs(c,j)) / cn_dom
           end do
      end do
      
end subroutine microbech4


subroutine get_waterhead(lbc, ubc, num_c, filter_c, waterhead)

	use clmtype
implicit none
	integer, intent(in) :: lbc, ubc 				! column-index bounds
	integer, intent(in) :: num_c          			! number of column soil points in column filter
	integer, intent(in) :: filter_c(ubc-lbc+1)		! column filter for soil points
	integer, intent(out)  :: waterhead(lbc:ubc)		! index of the soil layer right above the water table (-)!

	real(r8), pointer :: h2osoi_vol(:,:) 				! volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
	real(r8), pointer :: watsat(:,:)     				!volumetric soil water at saturation (porosity)
	real(r8), pointer :: t_soisno(:,:)   				! soil temperature (Kelvin)  (-nlevsno+1:nlevsoi)

	integer :: c,j,perch							! indices
	integer :: fc       							! filter column index

! borrowed from CH4 module created by Bill Riley
! to caluate the water head for unsaturated fraction of soil column

	h2osoi_vol				=> cws%h2osoi_vol
	watsat				=> cps%watsat
	t_soisno				=> ces%t_soisno

	do fc = 1, num_c
	c = filter_c(fc)
	perch = nlevsoi
	do j = nlevsoi, 1, -1
		if (t_soisno(c,j) < tfrz .and. h2osoi_vol(c,j) > f_sat * watsat(c,j)) then
		perch = j-1
		end if
	end do
	waterhead(c) = perch

	do j = perch, 2, -1
	if(h2osoi_vol(c,j) > f_sat * watsat(c,j) .and. h2osoi_vol(c,j-1) < f_sat * watsat(c,j-1)) then
		waterhead(c) = j-1
		exit
        end if
	enddo
	if (waterhead(c) == perch .and. h2osoi_vol(c,1) > f_sat * watsat(c,1)) then 
	waterhead(c) = 0
	endif
	end do
end subroutine get_waterhead

! Function for calculation of seasonality
!subroutine seasonality(lbc, ubc, lbp, ubp, num_micbioc, filter_micbioc, num_micbiop, filter_micbiop)
!~ !
!~ ! !DESCRIPTION: Annual mean fields.
!~ !
!~ ! !USES:
	!~ use clmtype
	!~ use clm_time_manager, only: get_step_size, get_days_per_year, get_nstep
	!~ use clm_varcon      , only: secspday
!~ !
!~ ! !ARGUMENTS:
!~ implicit none
	!~ integer, intent(in) :: lbc, ubc                		! column bounds
	!~ integer, intent(in) :: lbp, ubp                		! pft bounds
	!~ integer, intent(in) :: num_micbioc               	! number of soil columns in filter
	!~ integer, intent(in) :: filter_micbioc(ubc-lbc+1) 	! filter for soil columns
	!~ integer, intent(in) :: num_micbiop               	! number of soil points in pft filter
	!~ integer, intent(in) :: filter_micbiop(ubp-lbp+1) 	! pft filter for soil points
!~ !
!~ ! !LOCAL VARIABLES:
!~ ! local pointers to implicit in scalars
	!~ real(r8), pointer :: somhr(:)          			! (gC/m2/s) soil organic matter heterotrophic respiration
	!~ real(r8), pointer :: finundated(:)     			! fractional inundated area in soil column
	!~ real(r8), pointer :: agnpp(:)          			! (gC/m2/s) aboveground NPP
	!~ real(r8), pointer :: bgnpp(:)          			! (gC/m2/s) belowground NPP
	!~ integer , pointer :: pcolumn(:)        			! index into column level quantities
!~ !
!~ ! local pointers to implicit in/out scalars
!~ !
	!~ real(r8), pointer :: annsum_counter(:)  		! seconds since last annual accumulator turnover
!~ ! This will point to a CH4 version, not the CN version.
	!~ real(r8), pointer :: tempavg_somhr(:)   		! temporary average SOM heterotrophic resp. (gC/m2/s)
	!~ real(r8), pointer :: annavg_somhr(:)    		! annual average SOM heterotrophic resp. (gC/m2/s)
	!~ real(r8), pointer :: tempavg_finrw(:)   		! respiration-weighted annual average of finundated
	!~ real(r8), pointer :: annavg_finrw(:)    		! respiration-weighted annual average of finundated
	!~ real(r8), pointer :: tempavg_agnpp(:)   		! temporary average above-ground NPP (gC/m2/s)
	!~ real(r8), pointer :: annavg_agnpp(:)    		! annual average above-ground NPP (gC/m2/s)
	!~ real(r8), pointer :: tempavg_bgnpp(:)   		! temporary average below-ground NPP (gC/m2/s)
	!~ real(r8), pointer :: annavg_bgnpp(:)    		! annual average below-ground NPP (gC/m2/s)
	
	!~ real(r8), pointer :: bgnpp_timestep(:)    		! annual average below-ground NPP (gC/m2/s)
	!~ real(r8), pointer :: bgnpp_avg(:)    		! annual average below-ground NPP (gC/m2/s)
	!~ real(r8), pointer :: cgridcell(:)
!~ !
!~ ! local pointers to implicit out scalars
!~ !
!~ ! !OTHER LOCAL VARIABLES:
	!~ integer :: c,p,g       ! indices
	!~ integer :: fc        ! soil column filter indices
	!~ integer :: fp        ! soil pft filter indices
	!~ real(r8):: dt        ! time step (seconds)
	!~ real(r8):: secsperyear
	!~ logical :: newrun

!~ !EOP
!~ !-----------------------------------------------------------------------
!~ ! assign local pointers to derived type arrays
	!~ annsum_counter    			=> cmic%annsum_counter
	!~ tempavg_somhr     			=> cmic%tempavg_somhr
	!~ annavg_somhr      			=> cmic%annavg_somhr
	!~ tempavg_finrw     			=> cmic%tempavg_finrw
	!~ annavg_finrw      			=> cmic%annavg_finrw
	!~ pcolumn           				=> p%column
	!~ agnpp             				=> p%pcf%agnpp
	!~ bgnpp             				=> p%pcf%bgnpp
	!~ tempavg_agnpp     			=> p%pcf%tempavg_agnpp
	!~ annavg_agnpp      			=> p%pcf%annavg_agnpp
	!~ tempavg_bgnpp     			=> p%pcf%tempavg_bgnpp
	!~ annavg_bgnpp      			=> p%pcf%annavg_bgnpp
	!~ somhr             				=> ccf%somhr
	!~ finundated        				=> cws%finundated
	
	!~ bgnpp_timestep      			=> cmic%bgnpp_timestep
	!~ bgnpp_avg      				=> cmic%bgnpp_avg

!~ ! set time steps
	!~ dt = real(get_step_size(), r8)
	!~ secsperyear = real( get_days_per_year() * secspday, r8)

	!~ newrun = .false.

!~ ! column loop
	!~ do fc = 1,num_micbioc
	!~ c = filter_micbioc(fc)
	!~ if (annsum_counter(c) == spval) then
!~ ! These variables are now in restart files for completeness, but might not be in inicFile and are not.
!~ ! set for arbinit.
        !~ newrun = .true.
        !~ annsum_counter(c)    = 0._r8
        !~ tempavg_somhr(c)     = 0._r8
        !~ tempavg_finrw(c)     = 0._r8
	!~ end if

	!~ annsum_counter(c) = annsum_counter(c) + dt
	!~ end do

!~ ! pft loop
	!~ do fp = 1,num_micbiop
	!~ p = filter_micbiop(fp)
	!~ if (newrun .or. tempavg_agnpp(p) == spval) then ! Extra check needed because for back-compatibility
        !~ tempavg_agnpp(p) = 0._r8
        !~ tempavg_bgnpp(p) = 0._r8
	!~ end if
	!~ end do

	!~ do fc = 1,num_micbioc
	!~ c = filter_micbioc(fc)
	!~ if (annsum_counter(c) >= secsperyear) then
!~ ! update annual average somhr
        !~ annavg_somhr(c)      =  tempavg_somhr(c)
        !~ tempavg_somhr(c)     = 0._r8

!~ ! update annual average finrw
        !~ if (annavg_somhr(c) > 0._r8) then
        !~ annavg_finrw(c)      =  tempavg_finrw(c) / annavg_somhr(c)
        !~ else
        !~ annavg_finrw(c)      = 0._r8
        !~ end if
        !~ tempavg_finrw(c)     = 0._r8
	!~ else
        !~ tempavg_somhr(c)     = tempavg_somhr(c) + dt/secsperyear * somhr(c)
        !~ tempavg_finrw(c)     = tempavg_finrw(c) + dt/secsperyear * finundated(c) * somhr(c)
	!~ end if
	!~ end do
	
	!~ bgnpp_timestep(:) = 0._r8
	!~ bgnpp_avg(:) = 0._r8
	!~ do fp = 1,num_micbiop
	!~ p = filter_micbiop(fp)
	!~ c = pcolumn(p)
	!~ if (annsum_counter(c) >= secsperyear) then
        !~ annavg_agnpp(p) = tempavg_agnpp(p)
        !~ tempavg_agnpp(p) = 0._r8
        !~ annavg_bgnpp(p) = tempavg_bgnpp(p)
        !~ tempavg_bgnpp(p) = 0._r8
	!~ else
        !~ tempavg_agnpp(p) = tempavg_agnpp(p) + dt/secsperyear * agnpp(p)
        !~ tempavg_bgnpp(p) = tempavg_bgnpp(p) + dt/secsperyear * bgnpp(p)
	!~ end if
	!~ bgnpp_timestep(c) = bgnpp_timestep(c) + tempavg_bgnpp(p)
	!~ bgnpp_avg(c) = bgnpp_avg(c) + max(min(annavg_bgnpp(p),1e3), 1e-10)
!~ !write(iulog,*) c, "bgnpp_timestep(c)",bgnpp_timestep(c),"tempavg_bgnpp(p)",tempavg_bgnpp(p),"bgnpp_avg(c)",bgnpp_avg(c)
	!~ end do
   
!~ ! column loop
	!~ do fc = 1,num_micbioc
	!~ c = filter_micbioc(fc)
	!~ if (annsum_counter(c) >= secsperyear) annsum_counter(c) = 0._r8
	!~ end do
!end subroutine seasonality
! end of calculation of seasonality



subroutine update_finundated(lbc, ubc,num_micbioc,filter_micbioc)
!
! !DESCRIPTION: Annual mean fields.
!
! !USES:
	use clmtype
	use clm_time_manager, only: get_step_size, get_days_per_year, get_nstep
	use clm_varcon      , only: secspday
	use clm_varcon, only : secspday, istwet, istsoil, istdlak, spval, istcrop
!
! !ARGUMENTS:
implicit none
	integer, intent(in) :: lbc, ubc                		! column bounds
	integer, intent(in) :: num_micbioc               		! number of soil columns in filter
	integer, intent(in) :: filter_micbioc(ubc-lbc+1) 	! filter for soil columns
!
! !LOCAL VARIABLES:
! local pointers to implicit in scalars
	real(r8), pointer :: finundated(:)     			! fractional inundated area in soil column
	real(r8), pointer :: finundated_pre(:)     		! fractional inundated area in soil column in previous time step
	
	real(r8), pointer :: cdocs_unsat(:,:)                   	! column-level dissolved organic carbon in unsatured fraction
	real(r8), pointer :: cdocs_sat(:,:)                      	! column-level dissolved organic carbon in saturated fraction
   
	real(r8), pointer :: caces_unsat(:,:)			! column-level acetate in unsaturated fraction      
	real(r8), pointer :: caces_sat(:,:)				! column-level acetate in saturated fraction        
	real(r8), pointer :: cacebios_unsat(:,:)			! column-level biomass of methanogen from acetate in unsaturated fraction
	real(r8), pointer :: cacebios_sat(:,:)			! column-level biomass of methanogen from acetate in saturated fraction
	real(r8), pointer :: cco2bios_unsat(:,:)			! column-level biomass of methanogen from CO2 and H2 in unsaturated fraction
	real(r8), pointer :: cco2bios_sat(:,:)			! column-level biomass of methanogen from CO2 and H2 in saturated fraction
	real(r8), pointer :: caerch4bios_unsat(:,:)		! column-level biomass of aerobic methanotrophy 
	real(r8), pointer :: caerch4bios_sat(:,:)			! column-level biomass of aerobic methanotrophy in saturated fraction
	real(r8), pointer :: canaerch4bios_unsat(:,:)		! column-level biomass of anaerobic methanotrophy in unsaturated fraction
	real(r8), pointer :: canaerch4bios_sat(:,:)		! column-level biomass of anaerobic methanotrophy in saturated fraction
   
	real(r8), pointer :: ccon_ch4s_unsat(:,:)		! column-level concentration of CH4 in unsaturated fraction
	real(r8), pointer :: ccon_ch4s_sat(:,:)			! column-level concentration of CH4 in saturated fraction
	real(r8), pointer :: ccon_o2s_unsat(:,:)			! column-level concentration of O2 in unsaturated fraction
	real(r8), pointer :: ccon_o2s_sat(:,:)			! column-level concentration of O2 in saturated fraction
	real(r8), pointer :: ccon_co2s_unsat(:,:)			! column-level concentration of CO2 in unsaturated fraction
	real(r8), pointer :: ccon_co2s_sat(:,:)			! column-level concentration of CO2 in saturated fraction
	real(r8), pointer :: ccon_h2s_unsat(:,:)			! column-level concentration of H2 in unsaturated fraction
	real(r8), pointer :: ccon_h2s_sat(:,:)			! column-level concentration of H2 in saturated fraction
	
	integer, pointer :: ltype(:)
	integer, pointer :: clandunit(:)

!
! local pointers to implicit in/out scalars
! !OTHER LOCAL VARIABLES:
	integer :: c, j, l      ! indices
	integer :: fc        ! soil column filter indices
	real(r8) :: dlt
!EOP
!-----------------------------------------------------------------------
! assign local pointers to derived type arrays
	finundated        				=> cws%finundated
	finundated_pre        			=> cmic%fsat_pre

	cdocs_unsat 				=> cmic%cdocs_unsat
	cdocs_sat					=> cmic%cdocs_sat
	caces_unsat 				=> cmic%caces_unsat
	caces_sat 					=> cmic%caces_sat
	cacebios_unsat 				=> cmic%cacebios_unsat
	cacebios_sat 				=> cmic%cacebios_sat
	cco2bios_unsat 				=> cmic%cco2bios_unsat
	cco2bios_sat 				=> cmic%cco2bios_sat
	caerch4bios_unsat 			=> cmic%caerch4bios_unsat
	caerch4bios_sat 				=> cmic%caerch4bios_sat
	canaerch4bios_unsat 			=> cmic%canaerch4bios_unsat
	canaerch4bios_sat 			=> cmic%canaerch4bios_sat
     
	ccon_ch4s_unsat    			=> cmic%ccon_ch4s_unsat
	ccon_ch4s_sat        			=> cmic%ccon_ch4s_sat
	ccon_co2s_unsat    			=> cmic%ccon_co2s_unsat
	ccon_co2s_sat        			=> cmic%ccon_co2s_sat
	ccon_o2s_unsat      			=> cmic%ccon_o2s_unsat
	ccon_o2s_sat          			=> cmic%ccon_o2s_sat
	ccon_h2s_unsat      			=> cmic%ccon_h2s_unsat
	ccon_h2s_sat          			=> cmic%ccon_h2s_sat
	
	ltype               				=> lun%itype
	clandunit       				=> col%landunit
! set time steps

	do fc = 1,num_micbioc
        c = filter_micbioc(fc)
	if(finundated_pre(c) > 1.0) then
	finundated_pre(c) = finundated(c)
	end if
!write(iulog,*) " c: ", c
	l = clandunit(c)     
	if (ltype(l) == istsoil .or. ltype(l) == istcrop) then 
	do j = 1,nlevsoi
!write(iulog,*)	"before cdocs_sat", cdocs_sat(c,j), "cdocs_unsat", cdocs_unsat(c,j), "finundate", finundated(c), "finundated_pre(c)", finundated_pre(c)
	if(finundated(c) < finundated_pre(c)) then
	dlt 					= (finundated_pre(c) - finundated(c)) / finundated_pre(c)
	cdocs_unsat(c,j) 		= cdocs_unsat(c,j) + cdocs_sat(c,j) * dlt
	cdocs_sat(c,j) 			= cdocs_sat(c,j) - cdocs_sat(c,j) * dlt
	caces_unsat(c,j) 		= caces_unsat(c,j) + caces_sat(c,j) * dlt
	caces_sat(c,j) 			= caces_sat(c,j) - caces_sat(c,j) * dlt
	cacebios_unsat(c,j) 		= cacebios_unsat(c,j) + cacebios_sat(c,j) * dlt
	cacebios_sat(c,j) 		= cacebios_sat(c,j) - cacebios_sat(c,j) * dlt
	cco2bios_unsat(c,j) 		= cco2bios_unsat(c,j) + cco2bios_sat(c,j) * dlt
	cco2bios_sat(c,j) 		= cco2bios_sat(c,j) - cco2bios_sat(c,j) * dlt
	caerch4bios_unsat(c,j) 	= caerch4bios_unsat(c,j) + caerch4bios_sat(c,j) * dlt
	caerch4bios_sat(c,j) 		= caerch4bios_sat(c,j) - caerch4bios_sat(c,j) * dlt
	canaerch4bios_unsat(c,j) 	= canaerch4bios_unsat(c,j) + canaerch4bios_sat(c,j) * dlt
	canaerch4bios_sat(c,j) 	= canaerch4bios_sat(c,j) - canaerch4bios_sat(c,j) * dlt
	ccon_ch4s_unsat(c,j) 		= ccon_ch4s_unsat(c,j) + ccon_ch4s_sat(c,j) * dlt
	ccon_ch4s_sat(c,j) 		= ccon_ch4s_sat(c,j) - ccon_ch4s_sat(c,j) * dlt
	ccon_o2s_unsat(c,j) 		= ccon_o2s_unsat(c,j) + ccon_o2s_sat(c,j) * dlt
	ccon_o2s_sat(c,j) 		= ccon_o2s_sat(c,j) - ccon_o2s_sat(c,j) * dlt
	ccon_co2s_unsat(c,j) 		= ccon_co2s_unsat(c,j) + ccon_co2s_sat(c,j) * dlt
	ccon_co2s_sat(c,j) 		= ccon_co2s_sat(c,j) - ccon_co2s_sat(c,j) * dlt
	ccon_h2s_unsat(c,j) 		= ccon_h2s_unsat(c,j) + ccon_h2s_sat(c,j) * dlt
	ccon_h2s_sat 			= ccon_h2s_sat - ccon_h2s_sat * dlt
	else
	dlt 					= (finundated(c) - finundated_pre(c)) / (1. - finundated_pre(c))
	cdocs_sat(c,j) 			= cdocs_sat(c,j) + cdocs_unsat(c,j) * dlt
	cdocs_unsat(c,j) 		= cdocs_unsat(c,j) - cdocs_unsat(c,j) * dlt
	caces_sat(c,j) 			= caces_sat(c,j) + caces_unsat(c,j) * dlt
	caces_unsat(c,j) 		= caces_unsat(c,j) - caces_unsat(c,j) * dlt
	cacebios_sat(c,j) 		= cacebios_sat(c,j) + cacebios_unsat(c,j) * dlt
	cacebios_unsat(c,j) 		= cacebios_unsat(c,j) - cacebios_unsat(c,j) * dlt
	cco2bios_sat(c,j) 		= cco2bios_sat(c,j) + cco2bios_unsat(c,j) * dlt
	cco2bios_unsat(c,j) 		= cco2bios_unsat(c,j) - cco2bios_unsat(c,j) * dlt
	caerch4bios_sat(c,j) 		= caerch4bios_sat(c,j) + caerch4bios_unsat(c,j) * dlt
	caerch4bios_unsat(c,j) 	= caerch4bios_unsat(c,j) - caerch4bios_unsat(c,j) * dlt
	canaerch4bios_sat(c,j) 	= canaerch4bios_sat(c,j) + canaerch4bios_unsat(c,j) * dlt
	canaerch4bios_unsat(c,j) 	= canaerch4bios_unsat(c,j) - canaerch4bios_unsat(c,j) * dlt
	ccon_ch4s_sat(c,j) 		= ccon_ch4s_sat(c,j) + ccon_ch4s_unsat(c,j) * dlt
	ccon_ch4s_unsat(c,j) 		= ccon_ch4s_unsat(c,j) - ccon_ch4s_unsat(c,j) * dlt
	ccon_o2s_sat(c,j) 		= ccon_o2s_sat(c,j) + ccon_o2s_unsat(c,j) * dlt
	ccon_o2s_unsat(c,j) 		= ccon_o2s_unsat(c,j) - ccon_o2s_unsat(c,j) * dlt
	ccon_co2s_sat(c,j) 		= ccon_co2s_sat(c,j) + ccon_co2s_unsat(c,j) * dlt
	ccon_co2s_unsat(c,j) 		= ccon_co2s_unsat(c,j) - ccon_co2s_unsat(c,j) * dlt
	ccon_h2s_sat(c,j) 		= ccon_h2s_sat(c,j) + ccon_h2s_unsat(c,j) * dlt
	ccon_h2s_unsat(c,j) 		= ccon_h2s_unsat(c,j) - ccon_h2s_unsat(c,j) * dlt
	end if
!write(iulog,*)	"after cdocs_sat", cdocs_sat(c,j), "cdocs_unsat", cdocs_unsat(c,j), "finundate", finundated(c), "finundated_pre(c)", finundated_pre(c)

	end do
        end if
	finundated_pre(c) = finundated(c)
	end do

end subroutine update_finundated


subroutine gas_diffusion(lbc, ubc,num_micbioc,filter_micbioc)
! The gas diffusion was simulated based on Fick's law; the concentration gradient is the key driver for gas diffusion along soil profile
! The gas diffusion occur in saturated layers, the unsaturated layers was not considered
! !USES:
	use clmtype
	use clm_time_manager, only: get_step_size, get_days_per_year, get_nstep
	use clm_varcon      , only: secspday
	use clm_varcon, only : secspday, istwet, istsoil, istdlak, spval, istcrop
	use microbevarcon
	use clm_varpar, only : nlevgrnd, nlevdecomp
!
! !ARGUMENTS:
implicit none
	integer, intent(in) :: lbc, ubc                		! column bounds
	integer, intent(in) :: num_micbioc               		! number of soil columns in filter
	integer, intent(in) :: filter_micbioc(ubc-lbc+1) 	! filter for soil columns
!
! !LOCAL VARIABLES:
! local pointers to implicit in scalars
	real(r8), pointer :: z(:,:)         				! layer depth (m) (-nlevsno+1:nlevsoi)
	real(r8), pointer :: dz(:,:)         				! layer thickness (m) (-nlevsno+1:nlevgrnd)
	real(r8), pointer :: finundated(:)     			! fractional inundated area in soil column
	real(r8), pointer :: finundated_pre(:)     		! fractional inundated area in soil column in previous time step
  
	real(r8), pointer :: ccon_ch4s_unsat(:,:)		! column-level concentration of CH4 in unsaturated fraction
	real(r8), pointer :: ccon_ch4s_sat(:,:)			! column-level concentration of CH4 in saturated fraction
	real(r8), pointer :: ccon_o2s_unsat(:,:)			! column-level concentration of O2 in unsaturated fraction
	real(r8), pointer :: ccon_o2s_sat(:,:)			! column-level concentration of O2 in saturated fraction
	real(r8), pointer :: ccon_co2s_unsat(:,:)			! column-level concentration of CO2 in unsaturated fraction
	real(r8), pointer :: ccon_co2s_sat(:,:)			! column-level concentration of CO2 in saturated fraction
	real(r8), pointer :: ccon_h2s_unsat(:,:)			! column-level concentration of H2 in unsaturated fraction
	real(r8), pointer :: ccon_h2s_sat(:,:)			! column-level concentration of H2 in saturated fraction
	real(r8), pointer :: t_soisno(:,:)
	
	real(r8) :: ccon_ch4s_unsat_temp(lbc:ubc,1:nlevgrnd)	! temporary array 
	real(r8) :: ccon_ch4s_sat_temp(lbc:ubc,1:nlevgrnd)		! temporary array 
	real(r8) :: ccon_o2s_unsat_temp(lbc:ubc,1:nlevgrnd)	! temporary array 
	real(r8) :: ccon_o2s_sat_temp(lbc:ubc,1:nlevgrnd) 		! temporary array 
	real(r8) :: ccon_co2s_unsat_temp(lbc:ubc,1:nlevgrnd)	! temporary array 
	real(r8) :: ccon_co2s_sat_temp(lbc:ubc,1:nlevgrnd) 		! temporary array 
	real(r8) :: ccon_h2s_unsat_temp(lbc:ubc,1:nlevgrnd)	! temporary array 
	real(r8) :: ccon_h2s_sat_temp(lbc:ubc,1:nlevgrnd) 		! temporary array 
	real(r8) :: temp
	integer :: jwaterhead_unsat(lbc:ubc)			! layer of the water head in unsaturated fraction

	integer, pointer :: ltype(:)
	integer, pointer :: clandunit(:)
#if (defined HUM_HOL)
	real(r8) :: hum_frac
	real(r8) :: hol_frac

	real(r8) :: lxch4unsat
	real(r8) :: lxch4sat
	real(r8) :: lxco2unsat
	real(r8) :: lxco2sat
	real(r8) :: lxh2unsat
	real(r8) :: lxh2sat
	real(r8) :: lxo2unsat
	real(r8) :: lxo2sat
#endif
!
! local pointers to implicit in/out scalars
! !OTHER LOCAL VARIABLES:
	integer :: c, j, l      ! indices
	integer :: fc        ! soil column filter indices
	real(r8) :: dlt
!EOP
!-----------------------------------------------------------------------
! assign local pointers to derived type arrays
	z                      				=> cps%z
	dz                      			=> cps%dz
	finundated        				=> cws%finundated
	finundated_pre        			=> cmic%fsat_pre

	ccon_ch4s_unsat    			=> cmic%ccon_ch4s_unsat
	ccon_ch4s_sat        			=> cmic%ccon_ch4s_sat
	ccon_co2s_unsat    			=> cmic%ccon_co2s_unsat
	ccon_co2s_sat        			=> cmic%ccon_co2s_sat
	ccon_o2s_unsat      			=> cmic%ccon_o2s_unsat
	ccon_o2s_sat          			=> cmic%ccon_o2s_sat
	ccon_h2s_unsat      			=> cmic%ccon_h2s_unsat
	ccon_h2s_sat          			=> cmic%ccon_h2s_sat
	
	ltype               				=> lun%itype
	clandunit       				=> col%landunit
	t_soisno					=> ces%t_soisno
	
	ccon_ch4s_unsat_temp 		= 0._r8
	ccon_ch4s_sat_temp 			= 0._r8
	ccon_o2s_unsat_temp 		= 0._r8
	ccon_o2s_sat_temp 			= 0._r8
	ccon_co2s_unsat_temp 		= 0._r8
	ccon_co2s_sat_temp 			= 0._r8
	ccon_h2s_unsat_temp 		= 0._r8
	ccon_h2s_sat_temp 			= 0._r8

	call get_waterhead(lbc, ubc, num_micbioc, filter_micbioc,jwaterhead_unsat)
	
	do fc = 1,num_micbioc
        c = filter_micbioc(fc)
	l = clandunit(c)     
	if (ltype(l) == istsoil .or. ltype(l) == istcrop) then 
	do j = 2,nlevsoi
	if(j > jwaterhead_unsat(c) .and. j < nlevsoi) then
		ccon_ch4s_unsat_temp(c,j) = (ccon_ch4s_unsat(c,j-1) - ccon_ch4s_unsat(c,j)) * Fick_D_w(1) * m_Fick_ad * 1.0e-3 * (t_soisno(c,j)/298)**1.87 * get_step_size() / (z(c,j) - z(c,j-1)) !CH4_dif
		ccon_o2s_unsat_temp(c,j) = (ccon_o2s_unsat(c,j-1) - ccon_o2s_unsat(c,j)) * Fick_D_w(2) * m_Fick_ad * 1.0e-3 * (t_soisno(c,j)/298)**1.87 * get_step_size() / (z(c,j) - z(c,j-1)) !O2_dif
		ccon_co2s_unsat_temp(c,j) = (ccon_co2s_unsat(c,j-1) - ccon_co2s_unsat(c,j)) * Fick_D_w(3) * m_Fick_ad * 1.0e-3 *  (t_soisno(c,j)/298)**1.87 * get_step_size() / (z(c,j) - z(c,j-1)) !CO2_dif
		ccon_h2s_unsat_temp(c,j) = (ccon_h2s_unsat(c,j-1) - ccon_h2s_unsat(c,j)) * Fick_D_w(4) * m_Fick_ad * 1.0e-3 * (t_soisno(c,j)/298)**1.87 * get_step_size() / (z(c,j) - z(c,j-1)) !H2_dif
		
		ccon_ch4s_unsat(c,j-1) = (ccon_ch4s_unsat(c,j-1) * dz(c,j-1) - ccon_ch4s_unsat_temp(c,j)) / dz(c,j-1)
		ccon_ch4s_unsat(c,j) = (ccon_ch4s_unsat(c,j) * dz(c,j) + ccon_ch4s_unsat_temp(c,j)) / dz(c,j)
		
		ccon_o2s_unsat(c,j-1) = (ccon_o2s_unsat(c,j-1) * dz(c,j-1) - ccon_o2s_unsat_temp(c,j)) / dz(c,j-1)
		ccon_o2s_unsat(c,j) = (ccon_o2s_unsat(c,j) * dz(c,j) + ccon_o2s_unsat_temp(c,j)) / dz(c,j)

		ccon_co2s_unsat(c,j-1) = (ccon_co2s_unsat(c,j-1) * dz(c,j-1) - ccon_co2s_unsat_temp(c,j)) / dz(c,j-1)
		ccon_co2s_unsat(c,j) = (ccon_co2s_unsat(c,j) * dz(c,j) + ccon_co2s_unsat_temp(c,j)) / dz(c,j)

		ccon_h2s_unsat(c,j-1) = (ccon_h2s_unsat(c,j-1) * dz(c,j-1) - ccon_h2s_unsat_temp(c,j)) / dz(c,j-1)
		ccon_h2s_unsat(c,j) = (ccon_h2s_unsat(c,j) * dz(c,j) + ccon_h2s_unsat_temp(c,j)) / dz(c,j)
	end if
	! for the saturation portion
	!write(iulog,*) "before ", j, ccon_ch4s_sat(c,j-1), ccon_ch4s_sat(c,j)
		ccon_ch4s_sat_temp(c,j) = (ccon_ch4s_sat(c,j-1) - ccon_ch4s_sat(c,j)) * Fick_D_w(1) * m_Fick_ad * 1.0e-3 * (t_soisno(c,j)/298)**1.87 * get_step_size() / (z(c,j) - z(c,j-1)) !CH4_dif
		ccon_o2s_sat_temp(c,j) = (ccon_o2s_sat(c,j-1) - ccon_o2s_sat(c,j)) * Fick_D_w(2) * m_Fick_ad * 1.0e-3 * (t_soisno(c,j)/298)**1.87 * get_step_size() / (z(c,j) - z(c,j-1)) !O2_dif
		ccon_co2s_sat_temp(c,j) = (ccon_co2s_sat(c,j-1) - ccon_co2s_sat(c,j)) * Fick_D_w(3) * m_Fick_ad * 1.0e-3 *  (t_soisno(c,j)/298)**1.87 * get_step_size() / (z(c,j) - z(c,j-1)) !CO2_dif
		ccon_h2s_sat_temp(c,j) = (ccon_h2s_sat(c,j-1) - ccon_h2s_sat(c,j)) * Fick_D_w(4) * m_Fick_ad * 1.0e-3 * (t_soisno(c,j)/298)**1.87 * get_step_size() / (z(c,j) - z(c,j-1)) !H2_dif
!	!write(iulog,*) "temp ", ccon_ch4s_sat_temp(c,j)
		ccon_ch4s_sat(c,j-1) = (ccon_ch4s_sat(c,j-1) * dz(c,j-1) - ccon_ch4s_sat_temp(c,j)) / dz(c,j-1)
		ccon_ch4s_sat(c,j) = (ccon_ch4s_sat(c,j) * dz(c,j) + ccon_ch4s_sat_temp(c,j)) / dz(c,j)
		
		ccon_o2s_sat(c,j-1) = (ccon_o2s_sat(c,j-1) * dz(c,j-1) - ccon_o2s_sat_temp(c,j)) / dz(c,j-1)
		ccon_o2s_sat(c,j) = (ccon_o2s_sat(c,j) * dz(c,j) + ccon_o2s_sat_temp(c,j)) / dz(c,j)

		ccon_co2s_sat(c,j-1) = (ccon_co2s_sat(c,j-1) * dz(c,j-1) - ccon_co2s_sat_temp(c,j)) / dz(c,j-1)
		ccon_co2s_sat(c,j) = (ccon_co2s_sat(c,j) * dz(c,j) + ccon_co2s_sat_temp(c,j)) / dz(c,j)

		ccon_h2s_sat(c,j-1) = (ccon_h2s_sat(c,j-1) * dz(c,j-1) - ccon_h2s_sat_temp(c,j)) / dz(c,j-1)
		ccon_h2s_sat(c,j) = (ccon_h2s_sat(c,j) * dz(c,j) + ccon_h2s_sat_temp(c,j)) / dz(c,j)
	!write(iulog,*) "after ", ccon_ch4s_sat(c,j-1), ccon_ch4s_sat(c,j)
	end do
        end if
	end do

#if (defined HUM_HOL)
        hum_frac = 0.75
	hol_frac = 0.25

	do j = 1,nlevsoi
	if(j < nlevsoi) then
!write(iulog, *) "LBGC debug 1", lxch4unsat, ccon_ch4s_unsat(1,j), ccon_ch4s_unsat(2,j),  Fick_D_w(1), t_soisno(c,j), get_step_size()
		lxch4unsat = (ccon_ch4s_unsat(1,j) - ccon_ch4s_unsat(2,j)) * Fick_D_w(1) * m_Fick_ad * 1.0e-4 * (t_soisno(c,j) / 298)**1.87 * get_step_size() !CH4_dif
!write(iulog, *) "LBGC debug 2", lxch4unsat, ccon_ch4s_unsat(1,j), ccon_ch4s_unsat(2,j),  Fick_D_w(1), t_soisno(c,j), get_step_size()

		lxo2unsat = (ccon_o2s_unsat(1,j) - ccon_o2s_unsat(2,j)) * Fick_D_w(2) * m_Fick_ad * 1.0e-4 * (t_soisno(c,j) / 298)**1.87 * get_step_size() !O2_dif
		lxco2unsat = (ccon_co2s_unsat(1,j) - ccon_co2s_unsat(2,j)) * Fick_D_w(3) * m_Fick_ad * 1.0e-4 *  (t_soisno(c,j) / 298)**1.87 * get_step_size() !CO2_dif
		lxh2unsat = (ccon_h2s_unsat(1,j) - ccon_h2s_unsat(2,j)) * Fick_D_w(4) * m_Fick_ad * 1.0e-4 * (t_soisno(c,j) / 298)**1.87 * get_step_size() !H2_dif
		
!write(iulog, *) "LBGC debug 3",ccon_ch4s_unsat(1,j),ccon_ch4s_unsat(2,j)
		ccon_ch4s_unsat(1,j) = (ccon_ch4s_unsat(1,j) * hum_frac - lxch4unsat) / hum_frac
		ccon_ch4s_unsat(2,j) = (ccon_ch4s_unsat(2,j) * hol_frac + lxch4unsat) / hol_frac
!write(iulog, *) "LBGC debug 4",ccon_ch4s_unsat(1,j),ccon_ch4s_unsat(2,j)

		ccon_o2s_unsat(1,j) = (ccon_o2s_unsat(1,j) * hum_frac - lxo2unsat) / hum_frac
		ccon_o2s_unsat(2,j) = (ccon_o2s_unsat(2,j) * hol_frac + lxo2unsat) / hol_frac

		ccon_co2s_unsat(1,j) = (ccon_co2s_unsat(1,j) * hum_frac - lxco2unsat) / hum_frac
		ccon_co2s_unsat(2,j) = (ccon_co2s_unsat(2,j) * hol_frac + lxco2unsat) / hol_frac

		ccon_h2s_unsat(1,j) = (ccon_h2s_unsat(1,j) * hum_frac - lxh2unsat) / hum_frac
		ccon_h2s_unsat(2,j) = (ccon_h2s_unsat(2,j) * hol_frac + lxh2unsat) / hol_frac

	end if
		lxch4sat = (ccon_ch4s_sat(1,j) - ccon_ch4s_sat(2,j)) * Fick_D_w(1) * m_Fick_ad * 1.0e-4 * (t_soisno(c,j) / 298)**1.87 * get_step_size() !CH4_dif
		lxo2sat = (ccon_o2s_sat(1,j) - ccon_o2s_sat(2,j)) * Fick_D_w(2) * m_Fick_ad * 1.0e-4 * (t_soisno(c,j) / 298)**1.87 * get_step_size() !O2_dif
		lxco2sat = (ccon_co2s_sat(1,j) - ccon_co2s_sat(2,j)) * Fick_D_w(3) * m_Fick_ad * 1.0e-4 *  (t_soisno(c,j) / 298)**1.87 * get_step_size() !CO2_dif
		lxh2sat = (ccon_h2s_sat(1,j) - ccon_h2s_sat(2,j)) * Fick_D_w(4) * m_Fick_ad * 1.0e-4 * (t_soisno(c,j) / 298)**1.87 * get_step_size() !H2_dif
		
		ccon_ch4s_sat(1,j) = (ccon_ch4s_sat(1,j) * hum_frac - lxch4sat) / hum_frac
		ccon_ch4s_sat(2,j) = (ccon_ch4s_sat(2,j) * hol_frac + lxch4sat) / hol_frac

		ccon_o2s_sat(1,j) = (ccon_o2s_sat(1,j) * hum_frac - lxo2sat) / hum_frac
		ccon_o2s_sat(2,j) = (ccon_o2s_sat(2,j) * hol_frac + lxo2sat) / hol_frac

		ccon_co2s_sat(1,j) = (ccon_co2s_sat(1,j) * hum_frac - lxco2sat) / hum_frac
		ccon_co2s_sat(2,j) = (ccon_co2s_sat(2,j) * hol_frac + lxco2sat) / hol_frac

		ccon_h2s_sat(1,j) = (ccon_h2s_sat(1,j) * hum_frac - lxh2sat) / hum_frac
		ccon_h2s_sat(2,j) = (ccon_h2s_sat(2,j) * hol_frac + lxh2sat) / hol_frac
	end do

#endif

end subroutine gas_diffusion
! end of subroutine for simulating gas diffusion along water column in saturated soil profile

#if (defined HUM_HOL)
subroutine lateral_bgc(lbc, ubc,num_micbioc,filter_micbioc)
	use clmtype
	use clm_time_manager, only: get_step_size, get_days_per_year, get_nstep
	use clm_varcon , only: secspday
	use clm_varcon, only : secspday, istwet, istsoil, istdlak, spval, istcrop
	use microbevarcon
	use clm_varpar, only: nlevsoi, nlevgrnd, nlevdecomp
!
! !ARGUMENTS:
implicit none
	integer, intent(in) :: lbc, ubc                		! column bounds
	integer, intent(in) :: num_micbioc               		! number of soil columns in filter
	integer, intent(in) :: filter_micbioc(ubc-lbc+1) 	! filter for soil columns
	
	real(r8), pointer :: cdocs_sat(:,:)		! doc concetrnation in saturated fraction 
	real(r8), pointer :: cdons_sat(:,:)		! don concetrnation in saturated fraction 
	real(r8), pointer :: caces_sat(:,:)		! ace concetrnation in saturated fraction 
	real(r8), pointer :: ccon_o2s_sat(:,:)		! o2 concetrnation in saturated fraction 
	real(r8), pointer :: ccon_h2s_sat(:,:)		! h2 concetrnation in saturated fraction 
	real(r8), pointer :: ccon_co2s_sat(:,:)		! co2 concetrnation in saturated fraction 
	real(r8), pointer :: ccon_ch4s_sat(:,:)		! ch4 concetrnation in saturated fraction 
	real(r8), pointer :: h2osoi_liq(:,:)
	real(r8), pointer :: dz(:,:)              			! layer thickness (m)
	real(r8), pointer :: qflx_lat_aqu_layer(:,:)		! lateral transport for each layer
	real(r8), pointer :: t_soisno(:,:)

	real(r8) :: qflx_lat_layer1(1:20)
	real(r8) :: qflx_lat_layer2(1:14)
	real(r8) :: lxdocsat, lxdonsat, lxacesat, lxo2sat, lxh2sat, lxco2sat, lxch4sat
	real(r8) :: lx_o2s_sat_temp, lx_h2s_sat_temp, lx_co2s_sat_temp, lx_ch4s_sat_temp
	integer :: j
	
	real(r8) :: hum_frac
	real(r8) :: hol_frac
	real(r8) :: lat_flux_factor1
	real(r8) :: lat_flux_factor2
	integer :: hu_soil_id(1:20)
	real(r8) :: hu_soil_tk(1:20)
	real(r8) :: hu_soil_top(1:20)
	real(r8) :: hu_soil_bot(1:20)
	real(r8) :: hu_h2o(1:20)
	
	real(r8) :: cdocs_sat_spruce1(1:20)
	real(r8) :: cdons_sat_spruce1(1:20)
	real(r8) :: caces_sat_spruce1(1:20)
	real(r8) :: ccon_o2s_sat_spruce1(1:20)
	real(r8) :: ccon_h2s_sat_spruce1(1:20)
	real(r8) :: ccon_co2s_sat_spruce1(1:20)
	real(r8) :: ccon_ch4s_sat_spruce1(1:20)

	integer :: ho_soil_id(1:14)
	real(r8) :: ho_soil_tk(1:14)
	real(r8) :: ho_soil_top(1:14)
	real(r8) :: ho_soil_bot(1:14)
	real(r8) :: ho_h2o(1:14)

	real(r8) :: cdocs_sat_spruce2(1:20)
	real(r8) :: cdons_sat_spruce2(1:20)
	real(r8) :: caces_sat_spruce2(1:20)
	real(r8) :: ccon_o2s_sat_spruce2(1:20)
	real(r8) :: ccon_h2s_sat_spruce2(1:20)
	real(r8) :: ccon_co2s_sat_spruce2(1:20)
	real(r8) :: ccon_ch4s_sat_spruce2(1:20)

	real(r8) :: som_diffus = 1e-4_r8 / (secspday * 365._r8)  ! [m^2/sec] = 1 cm^2 / yr
!	real(r8) :: dom_diffus = 1e8_r8							! times to som diffus 

	dz                				=> cps%dz
	qflx_lat_aqu_layer 			=> cwf%qflx_lat_aqu_layer
	cdocs_sat					=> cmic%cdocs_sat
	cdons_sat					=> cmic%cdons_sat
	caces_sat					=> cmic%caces_sat
	ccon_ch4s_sat        			=> cmic%ccon_ch4s_sat
	ccon_co2s_sat        			=> cmic%ccon_co2s_sat
	ccon_o2s_sat          			=> cmic%ccon_o2s_sat
	ccon_h2s_sat          			=> cmic%ccon_h2s_sat
	h2osoi_liq        				=> cws%h2osoi_liq
	t_soisno 					=> ces%t_soisno
	
! new structure to old structure
	hu_soil_id(1) = 1
	hu_soil_id(2) = 2
	hu_soil_id(3) = 3
	hu_soil_id(4) = 4
	hu_soil_id(5) = 5
	hu_soil_id(6) = 6
	hu_soil_id(7) = 6
	hu_soil_id(8) = 6
	hu_soil_id(9) = 6
	hu_soil_id(10) = 6
	hu_soil_id(11) = 6
	hu_soil_id(12) = 7
	hu_soil_id(13) = 7
	hu_soil_id(14) = 7
	hu_soil_id(15) = 8
	hu_soil_id(16) = 8
	hu_soil_id(17) = 9
	hu_soil_id(18) = 9
	hu_soil_id(19) = 10
	hu_soil_id(20) = 10
! keep ten digital to keep the profile partitioning more accurate
	hu_soil_tk(1) = 0.0175128179162552
	hu_soil_tk(2) = 0.0275789692596763
	hu_soil_tk(3) = 0.0454700332424131
	hu_soil_tk(4) = 0.0749674109862083
	hu_soil_tk(5) = 0.1236003651022810
	hu_soil_tk(6) = 0.0108704034931660
	hu_soil_tk(7) = 0.0175128179162550
	hu_soil_tk(8) = 0.0275789692596760
	hu_soil_tk(9) = 0.0454700332424130
	hu_soil_tk(10) = 0.0749674109862090
	hu_soil_tk(11) = 0.0273829161127130
	hu_soil_tk(12) = 0.0962174489895681
	hu_soil_tk(13) = 0.2037825510104320
	hu_soil_tk(14) = 0.0359806264484320
	hu_soil_tk(15) = 0.3000000000000000
	hu_soil_tk(16) = 0.2539384053686900
	hu_soil_tk(17) = 0.3000000000000000
	hu_soil_tk(18) = 0.6132900315890600
	hu_soil_tk(19) = 0.3000000000000000
	hu_soil_tk(20) = 1.2057607013992800

	hu_soil_top(1) = 0.0
	hu_soil_top(2) = 0.0175128179162552
	hu_soil_top(3) = 0.0450917871759315
	hu_soil_top(4) = 0.0905618204183447
	hu_soil_top(5) = 0.165529231404553
	hu_soil_top(6) = 0.289129596506834
	hu_soil_top(7) = 0.300000000000000
	hu_soil_top(8) = 0.317512817916255
	hu_soil_top(9) = 0.345091787175931
	hu_soil_top(10) = 0.390561820418344
	hu_soil_top(11) = 0.465529231404553
	hu_soil_top(12) = 0.492912147517266
	hu_soil_top(13) = 0.589129596506834
	hu_soil_top(14) = 0.792912147517266
	hu_soil_top(15) = 0.828892773965698
	hu_soil_top(16) = 1.12889277396569
	hu_soil_top(17) = 1.38283117933438
	hu_soil_top(18) = 1.68283117933438
	hu_soil_top(19) = 2.29612121092344
	hu_soil_top(20) = 2.59612121092344

	hu_soil_bot(1) = 0.0175128179162552
	hu_soil_bot(2) = 0.0450917871759315
	hu_soil_bot(3) = 0.0905618204183447
	hu_soil_bot(4) = 0.165529231404553
	hu_soil_bot(5) = 0.289129596506834
	hu_soil_bot(6) = 0.300000000000000
	hu_soil_bot(7) = 0.317512817916255
	hu_soil_bot(8) = 0.345091787175931
	hu_soil_bot(9) = 0.390561820418344
	hu_soil_bot(10) = 0.465529231404553
	hu_soil_bot(11) = 0.492912147517266
	hu_soil_bot(12) = 0.589129596506834
	hu_soil_bot(13) = 0.792912147517266
	hu_soil_bot(14) = 0.828892773965698
	hu_soil_bot(15) = 1.12889277396569
	hu_soil_bot(16) = 1.38283117933438
	hu_soil_bot(17) = 1.68283117933438
	hu_soil_bot(18) = 2.29612121092344
	hu_soil_bot(19) = 2.59612121092344
	hu_soil_bot(20) = 3.80188191232272

	ho_soil_id(1) = 1
	ho_soil_id(2) = 2
	ho_soil_id(3) = 3
	ho_soil_id(4) = 4
	ho_soil_id(5) = 5
	ho_soil_id(6) = 5
	ho_soil_id(7) = 6
	ho_soil_id(8) = 7
	ho_soil_id(9) = 7
	ho_soil_id(10) = 8
	ho_soil_id(11) = 8
	ho_soil_id(12) = 9
	ho_soil_id(13) = 9
	ho_soil_id(14) = 10

	ho_soil_tk(1) = 0.0175128179162552
	ho_soil_tk(2) = 0.0275789692596763
	ho_soil_tk(3) = 0.0454700332424131
	ho_soil_tk(4) = 0.0749674109862083
	ho_soil_tk(5) = 0.0273829161127130
	ho_soil_tk(6) = 0.0962174489895680
	ho_soil_tk(7) = 0.2037825510104320
	ho_soil_tk(8) = 0.0359806264484320
	ho_soil_tk(9) = 0.3000000000000000
	ho_soil_tk(10) = 0.2539384053686820
	ho_soil_tk(11) = 0.3000000000000000
	ho_soil_tk(12) = 0.6132900315890600
	ho_soil_tk(13) = 0.3000000000000000
	ho_soil_tk(14) = 1.5057607013992800

	ho_soil_top(1) = 0
	ho_soil_top(2) = 0.0175128179162552
	ho_soil_top(3) = 0.0450917871759315
	ho_soil_top(4) = 0.0905618204183447
	ho_soil_top(5) = 0.165529231404553
	ho_soil_top(6) = 0.192912147517266
	ho_soil_top(7) = 0.289129596506834
	ho_soil_top(8) = 0.492912147517266
	ho_soil_top(9) = 0.528892773965698
	ho_soil_top(10) = 0.828892773965698
	ho_soil_top(11) = 1.08283117933438
	ho_soil_top(12) = 1.38283117933438
	ho_soil_top(13) = 1.99612121092344
	ho_soil_top(14) = 2.29612121092344

	ho_soil_bot(1) = 0.0175128179162552
	ho_soil_bot(2) = 0.0450917871759315
	ho_soil_bot(3) = 0.0905618204183447
	ho_soil_bot(4) = 0.165529231404553
	ho_soil_bot(5) = 0.192912147517266
	ho_soil_bot(6) = 0.289129596506834
	ho_soil_bot(7) = 0.492912147517266
	ho_soil_bot(8) = 0.528892773965698
	ho_soil_bot(9) = 0.828892773965698
	ho_soil_bot(10) = 1.08283117933438
	ho_soil_bot(11) = 1.38283117933438
	ho_soil_bot(12) = 1.99612121092344
	ho_soil_bot(13) = 2.29612121092344
	ho_soil_bot(14) = 3.80188191232272

	qflx_lat_layer1(:) = 0._r8
	qflx_lat_layer2(:) = 0._r8
	hu_h2o(:) = 0._r8
	ho_h2o(:) = 0._r8
	
	cdocs_sat_spruce1(:) = 0._r8
	cdons_sat_spruce1(:) = 0._r8
	caces_sat_spruce1(:) = 0._r8
	ccon_o2s_sat_spruce1(:) = 0._r8
	ccon_h2s_sat_spruce1(:) = 0._r8
	ccon_co2s_sat_spruce1(:) = 0._r8
	ccon_ch4s_sat_spruce1(:) = 0._r8
	
	cdocs_sat_spruce2(:) = 0._r8
	cdons_sat_spruce2(:) = 0._r8
	caces_sat_spruce2(:) = 0._r8
	ccon_o2s_sat_spruce2(:) = 0._r8
	ccon_h2s_sat_spruce2(:) = 0._r8
	ccon_co2s_sat_spruce2(:) = 0._r8
	ccon_ch4s_sat_spruce2(:) = 0._r8
	
	lxdocsat = 0._r8
	lxdonsat = 0._r8
	lxacesat = 0._r8
	lxo2sat = 0._r8
	lxh2sat = 0._r8
	lxco2sat = 0._r8
	lxch4sat = 0._r8
	lx_o2s_sat_temp = 0._r8
	lx_h2s_sat_temp = 0._r8
	lx_co2s_sat_temp = 0._r8
	lx_ch4s_sat_temp = 0._r8

	hum_frac = 0.75_r8
	hol_frac = 0.25_r8
	lat_flux_factor1 		= (hol_frac/hum_frac)
	lat_flux_factor2 		= (hum_frac/hol_frac)
	
!write(iulog,*) "herebefore hummock 10", cdocs_sat(1,1), cdocs_sat(1,2), cdocs_sat(1,3), cdocs_sat(1,4), cdocs_sat(1,5), cdocs_sat(1,6), cdocs_sat(1,7), cdocs_sat(1,8), cdocs_sat(1,9), cdocs_sat(1,10)
!write(iulog,*) "herebefore hollow 10", cdocs_sat(2,1), cdocs_sat(2,2), cdocs_sat(2,3), cdocs_sat(2,4), cdocs_sat(2,5), cdocs_sat(2,6), cdocs_sat(2,7), cdocs_sat(2,8), cdocs_sat(2,9), cdocs_sat(2,10)

do j = 1, 6
		qflx_lat_layer1(j) = qflx_lat_aqu_layer(1,hu_soil_id(j)) * hu_soil_tk(j) / dz(1,hu_soil_id(j))
		cdocs_sat_spruce1(j) = cdocs_sat(1,hu_soil_id(j)) !* hu_soil_tk(j) / dz(1,hu_soil_id(j))
		cdons_sat_spruce1(j) = cdons_sat(1,hu_soil_id(j)) !* hu_soil_tk(j) / dz(1,hu_soil_id(j))
		caces_sat_spruce1(j) = caces_sat(1,hu_soil_id(j)) !* hu_soil_tk(j) / dz(1,hu_soil_id(j))
		ccon_o2s_sat_spruce1(j) = ccon_o2s_sat(1,hu_soil_id(j)) !* hu_soil_tk(j) / dz(1,hu_soil_id(j))
		ccon_h2s_sat_spruce1(j) = ccon_h2s_sat(1,hu_soil_id(j)) !* hu_soil_tk(j) / dz(1,hu_soil_id(j))
		ccon_co2s_sat_spruce1(j) = ccon_co2s_sat(1,hu_soil_id(j)) !* hu_soil_tk(j) / dz(1,hu_soil_id(j))
		ccon_ch4s_sat_spruce1(j) = ccon_ch4s_sat(1,hu_soil_id(j)) !* hu_soil_tk(j) / dz(1,hu_soil_id(j))
end do

do j = 1, 14 ! allocate the 10 layers bgc variables to 14 layers in hollow and 20 layers in hummock
!write(iulog,*) "here", j, qflx_lat_aqu_layer(1,hu_soil_id(j + 6))
		qflx_lat_layer1(j+6) = qflx_lat_aqu_layer(1,hu_soil_id(j + 6)) * hu_soil_tk(j+6) / dz(1,hu_soil_id(j + 6))
		qflx_lat_layer2(j) = qflx_lat_aqu_layer(2,ho_soil_id(j)) * ho_soil_tk(j) / dz(2,ho_soil_id(j))
!write(iulog,*) "here2", j, qflx_lat_layer1(j+6), qflx_lat_layer2(j)		
		cdocs_sat_spruce1(j+6) = cdocs_sat(1,hu_soil_id(j + 6)) !* hu_soil_tk(j+6) / dz(1,hu_soil_id(j + 6))
		cdons_sat_spruce1(j+6) = cdons_sat(1,hu_soil_id(j + 6)) !* hu_soil_tk(j+6) / dz(1,hu_soil_id(j + 6))
		caces_sat_spruce1(j+6) = caces_sat(1,hu_soil_id(j + 6)) !* hu_soil_tk(j+6) / dz(1,hu_soil_id(j + 6))
		cdocs_sat_spruce2(j) = cdocs_sat(2,ho_soil_id(j)) !* ho_soil_tk(j) / dz(2,ho_soil_id(j))
		cdons_sat_spruce2(j) = cdons_sat(2,ho_soil_id(j)) !* ho_soil_tk(j) / dz(2,ho_soil_id(j))
		caces_sat_spruce2(j) = caces_sat(2,ho_soil_id(j)) !* ho_soil_tk(j) / dz(2,ho_soil_id(j))
		
		ccon_o2s_sat_spruce1(j+6) = ccon_o2s_sat(1,hu_soil_id(j + 6)) !* hu_soil_tk(j+6) / dz(1,hu_soil_id(j + 6))
		ccon_h2s_sat_spruce1(j+6) = ccon_h2s_sat(1,hu_soil_id(j + 6)) !* hu_soil_tk(j+6) / dz(1,hu_soil_id(j + 6))
		ccon_co2s_sat_spruce1(j+6) = ccon_co2s_sat(1,hu_soil_id(j + 6)) !* hu_soil_tk(j+6) / dz(1,hu_soil_id(j + 6))
		ccon_ch4s_sat_spruce1(j+6) = ccon_ch4s_sat(1,hu_soil_id(j + 6)) !* hu_soil_tk(j+6) / dz(1,hu_soil_id(j + 6))
		ccon_o2s_sat_spruce2(j) = ccon_o2s_sat(2,ho_soil_id(j)) !* ho_soil_tk(j) / dz(2,ho_soil_id(j))
		ccon_h2s_sat_spruce2(j) = ccon_h2s_sat(2,ho_soil_id(j)) !* ho_soil_tk(j) / dz(2,ho_soil_id(j))
		ccon_co2s_sat_spruce2(j) = ccon_co2s_sat(2,ho_soil_id(j)) !* ho_soil_tk(j) / dz(2,ho_soil_id(j))
		ccon_ch4s_sat_spruce2(j) = ccon_ch4s_sat(2,ho_soil_id(j)) !* ho_soil_tk(j) / dz(2,ho_soil_id(j))		
		
		hu_h2o(j+6) = h2osoi_liq(1,hu_soil_id(j+6)) * hu_soil_tk(j+6) / dz(1,hu_soil_id(j+6))
		ho_h2o(j) = h2osoi_liq(2,ho_soil_id(j)) * ho_soil_tk(j) / dz(2,ho_soil_id(j))
end do

!~ do j = 1, 6
		!~ cdocs_sat_spruce1(j) = cdocs_sat(1,hu_soil_id(j)) * hu_soil_tk(j) / dz(1,hu_soil_id(j))
		!~ cdons_sat_spruce1(j) = cdons_sat(1,hu_soil_id(j)) * hu_soil_tk(j) / dz(1,hu_soil_id(j))
		!~ caces_sat_spruce1(j) = caces_sat(1,hu_soil_id(j)) * hu_soil_tk(j) / dz(1,hu_soil_id(j))
		!~ ccon_o2s_sat_spruce1(j) = ccon_o2s_sat(1,hu_soil_id(j)) * hu_soil_tk(j) / dz(1,hu_soil_id(j))
		!~ ccon_h2s_sat_spruce1(j) = ccon_h2s_sat(1,hu_soil_id(j)) * hu_soil_tk(j) / dz(1,hu_soil_id(j))
		!~ ccon_co2s_sat_spruce1(j) = ccon_co2s_sat(1,hu_soil_id(j)) * hu_soil_tk(j) / dz(1,hu_soil_id(j))
		!~ ccon_ch4s_sat_spruce1(j) = ccon_ch4s_sat(1,hu_soil_id(j)) * hu_soil_tk(j) / dz(1,hu_soil_id(j))
!~ end do

!~ do j = 1, 14 ! allocate the 10 layers bgc variables to 14 layers in hollow and 20 layers in hummock
!~ !write(iulog,*) "here", j, qflx_lat_aqu_layer(1,hu_soil_id(j + 6))
		!~ qflx_lat_layer1(j+6) = qflx_lat_aqu_layer(1,hu_soil_id(j + 6)) * hu_soil_tk(j+6) / dz(1,hu_soil_id(j + 6))
		!~ qflx_lat_layer2(j) = qflx_lat_aqu_layer(2,ho_soil_id(j)) * ho_soil_tk(j) / dz(2,ho_soil_id(j))
!~ !write(iulog,*) "here2", j, qflx_lat_layer1(j+6), qflx_lat_layer2(j)		
		!~ cdocs_sat_spruce1(j+6) = cdocs_sat(1,hu_soil_id(j + 6)) * hu_soil_tk(j+6) / dz(1,hu_soil_id(j + 6))
		!~ cdons_sat_spruce1(j+6) = cdons_sat(1,hu_soil_id(j + 6)) * hu_soil_tk(j+6) / dz(1,hu_soil_id(j + 6))
		!~ caces_sat_spruce1(j+6) = caces_sat(1,hu_soil_id(j + 6)) * hu_soil_tk(j+6) / dz(1,hu_soil_id(j + 6))
		!~ cdocs_sat_spruce2(j) = cdocs_sat(2,ho_soil_id(j)) * ho_soil_tk(j) / dz(2,ho_soil_id(j))
		!~ cdons_sat_spruce2(j) = cdons_sat(2,ho_soil_id(j)) * ho_soil_tk(j) / dz(2,ho_soil_id(j))
		!~ caces_sat_spruce2(j) = caces_sat(2,ho_soil_id(j)) * ho_soil_tk(j) / dz(2,ho_soil_id(j))
		
		!~ ccon_o2s_sat_spruce1(j+6) = ccon_o2s_sat(1,hu_soil_id(j + 6)) * hu_soil_tk(j+6) / dz(1,hu_soil_id(j + 6))
		!~ ccon_h2s_sat_spruce1(j+6) = ccon_h2s_sat(1,hu_soil_id(j + 6)) * hu_soil_tk(j+6) / dz(1,hu_soil_id(j + 6))
		!~ ccon_co2s_sat_spruce1(j+6) = ccon_co2s_sat(1,hu_soil_id(j + 6)) * hu_soil_tk(j+6) / dz(1,hu_soil_id(j + 6))
		!~ ccon_ch4s_sat_spruce1(j+6) = ccon_ch4s_sat(1,hu_soil_id(j + 6)) * hu_soil_tk(j+6) / dz(1,hu_soil_id(j + 6))
		!~ ccon_o2s_sat_spruce2(j) = ccon_o2s_sat(2,ho_soil_id(j)) * ho_soil_tk(j) / dz(2,ho_soil_id(j))
		!~ ccon_h2s_sat_spruce2(j) = ccon_h2s_sat(2,ho_soil_id(j)) * ho_soil_tk(j) / dz(2,ho_soil_id(j))
		!~ ccon_co2s_sat_spruce2(j) = ccon_co2s_sat(2,ho_soil_id(j)) * ho_soil_tk(j) / dz(2,ho_soil_id(j))
		!~ ccon_ch4s_sat_spruce2(j) = ccon_ch4s_sat(2,ho_soil_id(j)) * ho_soil_tk(j) / dz(2,ho_soil_id(j))		
		
		!~ hu_h2o(j+6) = h2osoi_liq(1,hu_soil_id(j+6)) * hu_soil_tk(j+6) / dz(1,hu_soil_id(j+6))
		!~ ho_h2o(j) = h2osoi_liq(2,ho_soil_id(j)) * ho_soil_tk(j) / dz(2,ho_soil_id(j))
!~ end do

!write(iulog,*) "herebefore hummock 20", cdocs_sat_spruce1(1), cdocs_sat_spruce1(2), cdocs_sat_spruce1(3), cdocs_sat_spruce1(4), cdocs_sat_spruce1(5), cdocs_sat_spruce1(6), &
!cdocs_sat_spruce1(7), cdocs_sat_spruce1(8), cdocs_sat_spruce1(9), cdocs_sat_spruce1(10), cdocs_sat_spruce1(11), cdocs_sat_spruce1(12), cdocs_sat_spruce1(13), cdocs_sat_spruce1(14), &
!cdocs_sat_spruce1(15), cdocs_sat_spruce1(16), cdocs_sat_spruce1(17), cdocs_sat_spruce1(18), cdocs_sat_spruce1(19), cdocs_sat_spruce1(20)
!write(iulog,*) "herebefore hollow 14", cdocs_sat_spruce2(1), cdocs_sat_spruce2(2), cdocs_sat_spruce2(3), cdocs_sat_spruce2(4), cdocs_sat_spruce2(5), cdocs_sat_spruce2(6), &
!cdocs_sat_spruce2(7), cdocs_sat_spruce2(8), cdocs_sat_spruce2(9), cdocs_sat_spruce2(10), cdocs_sat_spruce2(11), cdocs_sat_spruce2(12), cdocs_sat_spruce2(13), cdocs_sat_spruce2(14), &

!write(iulog,*) "herebefore", cdocs_sat_spruce2(1)+cdocs_sat_spruce2(2)+cdocs_sat_spruce2(3)+cdocs_sat_spruce2(4)+cdocs_sat_spruce2(5)+cdocs_sat_spruce2(6)+&
!cdocs_sat_spruce2(7)+cdocs_sat_spruce2(8)+cdocs_sat_spruce2(9)+cdocs_sat_spruce2(10)+cdocs_sat_spruce2(11)+cdocs_sat_spruce2(12)+cdocs_sat_spruce2(13)+&
!cdocs_sat_spruce2(14)

! lateral transport
!~ if(jwt(c) >= nlevsoi) then
!~ else
	do j = 1, 14
		if(abs(qflx_lat_layer1(j+6))> 0) then
		lxdocsat = cdocs_sat_spruce1(j+6) * qflx_lat_layer1(j+6)/hu_h2o(j+6)/1E3	! mm/s for one time step in qflx_lat_layer to consistent with m
	if(lxdocsat > 0 .and. abs(lxdocsat) > (cdocs_sat_spruce1(j+6) / sqrt(lat_flux_factor1))) then
	lxdocsat = cdocs_sat_spruce1(j+6) / sqrt(lat_flux_factor1) * 0.9
	end if
	if(lxdocsat < 0 .and. abs(lxdocsat) > (cdocs_sat_spruce2(j) / sqrt(lat_flux_factor2))) then
	lxdocsat = -cdocs_sat_spruce2(j) / sqrt(lat_flux_factor2) * 0.9
	end if
	
		cdocs_sat_spruce1(j+6) = cdocs_sat_spruce1(j+6) - lxdocsat * sqrt(lat_flux_factor1)
		cdocs_sat_spruce2(j) = cdocs_sat_spruce2(j) + lxdocsat * sqrt(lat_flux_factor2)
		
		lxdonsat = cdons_sat_spruce1(j+6) * qflx_lat_layer1(j+6)/hu_h2o(j+6)/1E3
	if(lxdonsat > 0 .and. abs(lxdonsat) > (cdons_sat_spruce1(j+6) / sqrt(lat_flux_factor1))) then
	lxdonsat = cdons_sat_spruce1(j+6) / sqrt(lat_flux_factor1) * 0.9
	end if
	if(lxdonsat < 0 .and. abs(lxdonsat) > (cdons_sat_spruce2(j)/ sqrt(lat_flux_factor2))) then
	lxdonsat = -cdons_sat_spruce2(j) / sqrt(lat_flux_factor2) * 0.9
	end if		
		cdons_sat_spruce1(j+6) = cdons_sat_spruce1(j+6) - lxdonsat * sqrt(lat_flux_factor1)
		cdons_sat_spruce2(j) = cdons_sat_spruce2(j) + lxdonsat * sqrt(lat_flux_factor2)
		
		lxacesat = caces_sat_spruce1(j+6) * qflx_lat_layer1(j+6)/hu_h2o(j+6)/1E3
	if(lxacesat > 0 .and. abs(lxacesat) > (caces_sat_spruce1(j+6) / sqrt(lat_flux_factor1))) then
	lxacesat = caces_sat_spruce1(j+6) / sqrt(lat_flux_factor1) * 0.9
	end if
	if(lxacesat < 0 .and. abs(lxacesat) > (caces_sat_spruce2(j) /sqrt(lat_flux_factor2))) then
	lxacesat = -caces_sat_spruce2(j) / sqrt(lat_flux_factor2) * 0.9
	end if		
		caces_sat_spruce1(j+6) = caces_sat_spruce1(j+6) - lxacesat * sqrt(lat_flux_factor1)
		caces_sat_spruce2(j) = caces_sat_spruce2(j) + lxacesat * sqrt(lat_flux_factor2)
			
		lxo2sat = ccon_o2s_sat_spruce1(j+6) * qflx_lat_layer1(j+6)/hu_h2o(j+6)/1E3
	if(lx_o2s_sat_temp > 0 .and. abs(lx_o2s_sat_temp) > (ccon_o2s_sat_spruce1(j+6) / sqrt(lat_flux_factor1))) then
	lx_o2s_sat_temp = ccon_o2s_sat_spruce1(j+6) / sqrt(lat_flux_factor1) * 0.9
	end if
	if(lx_o2s_sat_temp < 0 .and. abs(lx_o2s_sat_temp) > (ccon_o2s_sat_spruce2(j) / sqrt(lat_flux_factor2))) then
	lx_o2s_sat_temp = -ccon_o2s_sat_spruce2(j) / sqrt(lat_flux_factor2) * 0.9
	end if	
		ccon_o2s_sat_spruce1(j+6) = ccon_o2s_sat_spruce1(j+6) - lxo2sat * sqrt(lat_flux_factor1)
		ccon_o2s_sat_spruce2(j) = ccon_o2s_sat_spruce2(j) + lxo2sat * sqrt(lat_flux_factor2)
		
		lxh2sat = ccon_h2s_sat_spruce1(j+6) * qflx_lat_layer1(j+6)/hu_h2o(j+6)/1E3
	if(lx_h2s_sat_temp > 0 .and. abs(lx_h2s_sat_temp) > (ccon_h2s_sat_spruce1(j+6) / sqrt(lat_flux_factor1))) then
	lx_h2s_sat_temp = ccon_h2s_sat_spruce1(j+6) / sqrt(lat_flux_factor1) * 0.9
	end if
	if(lx_h2s_sat_temp < 0 .and. abs(lx_h2s_sat_temp) > (ccon_h2s_sat_spruce2(j) / sqrt(lat_flux_factor2))) then
	lx_h2s_sat_temp = -ccon_h2s_sat_spruce2(j) / sqrt(lat_flux_factor2) * 0.9
	end if	
		ccon_h2s_sat_spruce1(j+6) = ccon_h2s_sat_spruce1(j+6) - lxh2sat * sqrt(lat_flux_factor1)
		ccon_h2s_sat_spruce2(j) = ccon_h2s_sat_spruce2(j) + lxh2sat * sqrt(lat_flux_factor2)
		
		lxco2sat = ccon_co2s_sat_spruce1(j+6) * qflx_lat_layer1(j+6)/hu_h2o(j+6)/1E3
	if(lx_co2s_sat_temp > 0 .and. abs(lx_co2s_sat_temp) > (ccon_co2s_sat_spruce1(j+6) / sqrt(lat_flux_factor1))) then
	lx_co2s_sat_temp = ccon_co2s_sat_spruce1(j+6) / sqrt(lat_flux_factor1) * 0.9
	end if
	if(lx_co2s_sat_temp < 0 .and. abs(lx_co2s_sat_temp) > (ccon_co2s_sat_spruce2(j) / sqrt(lat_flux_factor2))) then
	lx_co2s_sat_temp = -ccon_co2s_sat_spruce2(j) / sqrt(lat_flux_factor2) * 0.9
	end if	
		ccon_co2s_sat_spruce1(j+6) = ccon_co2s_sat_spruce1(j+6) - lxco2sat * sqrt(lat_flux_factor1)
		ccon_co2s_sat_spruce2(j) = ccon_co2s_sat_spruce2(j) + lxco2sat * sqrt(lat_flux_factor2)

		lxch4sat = ccon_ch4s_sat_spruce1(j+6) * qflx_lat_layer1(j+6)/hu_h2o(j+6)/1E3
	if(lx_ch4s_sat_temp > 0 .and. abs(lx_ch4s_sat_temp) > (ccon_ch4s_sat_spruce1(j+6)/ sqrt(lat_flux_factor1))) then
	lx_ch4s_sat_temp = ccon_ch4s_sat_spruce1(j+6)/ sqrt(lat_flux_factor1) * 0.9
	end if
	if(lx_ch4s_sat_temp < 0 .and. abs(lx_ch4s_sat_temp) > (ccon_ch4s_sat_spruce2(j)/ sqrt(lat_flux_factor2))) then
	lx_ch4s_sat_temp = -ccon_ch4s_sat_spruce2(j)/ sqrt(lat_flux_factor2) * 0.9
	end if
		ccon_ch4s_sat_spruce1(j+6) = ccon_ch4s_sat_spruce1(j+6) - lxch4sat * sqrt(lat_flux_factor1)
		ccon_ch4s_sat_spruce2(j) = ccon_ch4s_sat_spruce2(j) + lxch4sat * sqrt(lat_flux_factor2)
		end if
	end do	
!~ end if

!~ write(iulog,*) "after laterer transport hummock", cdocs_sat_spruce1(1), cdocs_sat_spruce1(2), cdocs_sat_spruce1(3), cdocs_sat_spruce1(4), cdocs_sat_spruce1(5), &
!cdocs_sat_spruce1(6), cdocs_sat_spruce1(7), cdocs_sat_spruce1(8), cdocs_sat_spruce1(9), cdocs_sat_spruce1(10), cdocs_sat_spruce1(11), cdocs_sat_spruce1(12), &
!cdocs_sat_spruce1(13), cdocs_sat_spruce1(14), cdocs_sat_spruce1(15), cdocs_sat_spruce1(16), cdocs_sat_spruce1(17), cdocs_sat_spruce1(18), cdocs_sat_spruce1(19), cdocs_sat_spruce1(20)
!~ write(iulog,*) "after lateral transport hollow", cdocs_sat_spruce2(1), cdocs_sat_spruce2(2), cdocs_sat_spruce2(3), cdocs_sat_spruce2(4), cdocs_sat_spruce2(5), &
!cdocs_sat_spruce2(6), cdocs_sat_spruce2(7), cdocs_sat_spruce2(8), cdocs_sat_spruce2(9), cdocs_sat_spruce2(10), cdocs_sat_spruce2(11), cdocs_sat_spruce2(12), &
!cdocs_sat_spruce2(13), cdocs_sat_spruce2(14)

! lateral diffusion
do j = 1, 14
		lxdocsat 			= (cdocs_sat_spruce1(j+6) - cdocs_sat_spruce2(j)) * dom_diffus * get_step_size()
		lxdonsat 			= (cdons_sat_spruce1(j+6) - cdons_sat_spruce2(j)) * dom_diffus * get_step_size()
		lxacesat 			= (caces_sat_spruce1(j+6) - caces_sat_spruce2(j)) * dom_diffus * get_step_size()
		!~ lxdocsat 			= (cdocs_sat_spruce1(j+6)/hu_h2o(j+6) - cdocs_sat_spruce2(j)/ho_h2o(j)) * dom_diffus * get_step_size()
		!~ lxdonsat 			= (cdocs_sat_spruce1(j+6)/hu_h2o(j+6) - cdocs_sat_spruce2(j)/ho_h2o(j)) * dom_diffus * get_step_size()
		!~ lxacesat 			= (caces_sat_spruce1(j+6)/hu_h2o(j+6) - caces_sat_spruce2(j)/ho_h2o(j)) * dom_diffus * get_step_size()
!write(iulog,*) "inside0: ", j, lxdocsat, cdocs_sat_spruce1(j+6), cdocs_sat_spruce2(j), hu_h2o(j+6), ho_h2o(j)
	if(lxdocsat > 0 .and. abs(lxdocsat) > (cdocs_sat_spruce1(j+6)/sqrt(lat_flux_factor1))) then
	lxdocsat = cdocs_sat_spruce1(j+6)/sqrt(lat_flux_factor1)*0.9
	end if
	if(lxdocsat < 0 .and. abs(lxdocsat) > (cdocs_sat_spruce2(j)/sqrt(lat_flux_factor2))) then
	lxdocsat = -cdocs_sat_spruce2(j)/sqrt(lat_flux_factor2)*0.9
	end if
		cdocs_sat_spruce1(j+6) 	= (cdocs_sat_spruce1(j+6) - lxdocsat * sqrt(lat_flux_factor1))! * hu_h2o(j+6) 
		cdocs_sat_spruce2(j) 		= (cdocs_sat_spruce2(j) + lxdocsat * sqrt(lat_flux_factor2))! * hu_h2o(j+6) 
		
!write(iulog,*) "inside1: ", j, lxdocsat, cdocs_sat_spruce1(j+6), cdocs_sat_spruce2(j), hu_h2o(j+6), ho_h2o(j)
	if(lxdonsat > 0 .and. abs(lxdonsat) > (cdons_sat_spruce1(j+6)/sqrt(lat_flux_factor1))) then
	lxdonsat = cdons_sat_spruce1(j+6)/sqrt(lat_flux_factor1)*0.9
	end if
	if(lxdonsat < 0 .and. abs(lxdonsat) > (cdons_sat_spruce2(j)/sqrt(lat_flux_factor2))) then
	lxdonsat = -cdons_sat_spruce2(j)/sqrt(lat_flux_factor2)*0.9
	end if
		cdons_sat_spruce1(j+6) 	= (cdons_sat_spruce1(j+6) - lxdonsat * sqrt(lat_flux_factor1)) !* hu_h2o(j+6) 
		cdons_sat_spruce2(j) 		= (cdons_sat_spruce2(j) + lxdonsat * sqrt(lat_flux_factor2)) !* hu_h2o(j+6) 
		
	if(lxacesat > 0 .and. abs(lxacesat) > (caces_sat_spruce1(j+6)/sqrt(lat_flux_factor1))) then
	lxacesat = caces_sat_spruce1(j+6)/sqrt(lat_flux_factor1)*0.9
	end if
	if(lxacesat < 0 .and. abs(lxacesat) > (caces_sat_spruce2(j)/sqrt(lat_flux_factor2))) then
	lxacesat = -caces_sat_spruce2(j)/sqrt(lat_flux_factor2)*0.9
	end if
		caces_sat_spruce1(j+6) 	= (caces_sat_spruce1(j+6) - lxacesat * sqrt(lat_flux_factor1)) !* hu_h2o(j+6)  
		caces_sat_spruce2(j) 		= (caces_sat_spruce2(j) + lxacesat * sqrt(lat_flux_factor2)) !* hu_h2o(j+6) 
		
		lx_ch4s_sat_temp		= (ccon_ch4s_sat_spruce1(j+6)-ccon_ch4s_sat_spruce2(j)) * Fick_D_w(1) * 1e-4 * (t_soisno(1,j)/298)**1.87 * get_step_size() 	! CH4_dif
		lx_o2s_sat_temp 			= (ccon_o2s_sat_spruce1(j+6) - ccon_o2s_sat_spruce2(j)) * Fick_D_w(2) * 1e-4 * (t_soisno(1,j)/298)**1.87 * get_step_size() 	! O2_dif
		lx_co2s_sat_temp 		= (ccon_co2s_sat_spruce1(j+6)-ccon_co2s_sat_spruce2(j)) * Fick_D_w(3) * 1e-4 *  (t_soisno(1,j)/298)**1.87 * get_step_size() 	! CO2_dif
		lx_h2s_sat_temp 			= (ccon_h2s_sat_spruce1(j+6) - ccon_h2s_sat_spruce2(j)) * Fick_D_w(4) * 1e-4 * (t_soisno(1,j)/298)**1.87 * get_step_size() 	! H2_dif

		!~ lx_ch4s_sat_temp		= (ccon_ch4s_sat_spruce1(j+6)/hu_h2o(j+6) - ccon_ch4s_sat_spruce2(j)/ho_h2o(j)) * Fick_D_w(1) * 1e-4 * (t_soisno(1,j)/298)**1.87 * get_step_size() 	! CH4_dif
		!~ lx_o2s_sat_temp 			= (ccon_o2s_sat_spruce1(j+6)/hu_h2o(j+6) - ccon_o2s_sat_spruce2(j)/ho_h2o(j)) * Fick_D_w(2) * 1e-4 * (t_soisno(1,j)/298)**1.87 * get_step_size() 	! O2_dif
		!~ lx_co2s_sat_temp 		= (ccon_co2s_sat_spruce1(j+6)/hu_h2o(j+6) - ccon_co2s_sat_spruce2(j)/ho_h2o(j)) * Fick_D_w(3) * 1e-4 *  (t_soisno(1,j)/298)**1.87 * get_step_size() 	! CO2_dif
		!~ lx_h2s_sat_temp 			= (ccon_h2s_sat_spruce1(j+6)/hu_h2o(j+6) - ccon_h2s_sat_spruce2(j)/ho_h2o(j)) * Fick_D_w(4) * 1e-4 * (t_soisno(1,j)/298)**1.87 * get_step_size() 	! H2_dif

	if(lx_ch4s_sat_temp > 0 .and. abs(lx_ch4s_sat_temp) > (ccon_ch4s_sat_spruce1(j+6)/sqrt(lat_flux_factor1))) then
	lx_ch4s_sat_temp = ccon_ch4s_sat_spruce1(j+6)/sqrt(lat_flux_factor1)*0.9
	end if
	if(lx_ch4s_sat_temp < 0 .and. abs(lx_ch4s_sat_temp) > (ccon_ch4s_sat_spruce2(j)/sqrt(lat_flux_factor2))) then
	lx_ch4s_sat_temp = -ccon_ch4s_sat_spruce2(j)/sqrt(lat_flux_factor2)*0.9
	end if
		ccon_ch4s_sat_spruce1(j+6) = (ccon_ch4s_sat_spruce1(j+6) - lx_ch4s_sat_temp * sqrt(lat_flux_factor1)) !* hu_h2o(j+6)  
		ccon_ch4s_sat_spruce2(j) 	= (ccon_ch4s_sat_spruce2(j) + lx_ch4s_sat_temp * sqrt(lat_flux_factor2)) !* hu_h2o(j+6) 
		
	if(lx_o2s_sat_temp > 0 .and. abs(lx_o2s_sat_temp) > (ccon_o2s_sat_spruce1(j+6)/sqrt(lat_flux_factor1))) then
	lx_o2s_sat_temp = ccon_o2s_sat_spruce1(j+6)/sqrt(lat_flux_factor1)*0.9
	end if
	if(lx_o2s_sat_temp < 0 .and. abs(lx_o2s_sat_temp) > (ccon_o2s_sat_spruce2(j)/sqrt(lat_flux_factor2))) then
	lx_o2s_sat_temp = -ccon_o2s_sat_spruce2(j)/sqrt(lat_flux_factor2)*0.9
	end if
		ccon_o2s_sat_spruce1(j+6) = (ccon_o2s_sat_spruce1(j+6) - lx_o2s_sat_temp * sqrt(lat_flux_factor1)) !* hu_h2o(j+6) 
		ccon_o2s_sat_spruce2(j) 	= (ccon_o2s_sat_spruce2(j) + lx_o2s_sat_temp * sqrt(lat_flux_factor2)) !* hu_h2o(j+6) 
		
	if(lx_co2s_sat_temp > 0 .and. abs(lx_co2s_sat_temp) > (ccon_co2s_sat_spruce1(j+6)/sqrt(lat_flux_factor1))) then
	lx_co2s_sat_temp = ccon_co2s_sat_spruce1(j+6)/sqrt(lat_flux_factor1)*0.9
	end if
	if(lx_co2s_sat_temp < 0 .and. abs(lx_co2s_sat_temp) > (ccon_co2s_sat_spruce2(j)/sqrt(lat_flux_factor2))) then
	lx_co2s_sat_temp = -ccon_co2s_sat_spruce2(j)/sqrt(lat_flux_factor2)*0.9
	end if
		ccon_co2s_sat_spruce1(j+6) = (ccon_co2s_sat_spruce1(j+6) - lx_co2s_sat_temp * sqrt(lat_flux_factor1)) !* hu_h2o(j+6)  
		ccon_co2s_sat_spruce2(j) 	= (ccon_co2s_sat_spruce2(j) + lx_co2s_sat_temp * sqrt(lat_flux_factor2)) !* hu_h2o(j+6) 
		
	if(lx_h2s_sat_temp > 0 .and. abs(lx_h2s_sat_temp) > (ccon_h2s_sat_spruce1(j+6)/sqrt(lat_flux_factor1))) then
	lx_h2s_sat_temp = ccon_h2s_sat_spruce1(j+6)/sqrt(lat_flux_factor1)*0.9
	end if
	if(lx_h2s_sat_temp < 0 .and. abs(lx_h2s_sat_temp) > (ccon_h2s_sat_spruce2(j)/sqrt(lat_flux_factor2))) then
	lx_h2s_sat_temp = -ccon_h2s_sat_spruce2(j)/sqrt(lat_flux_factor2) * 0.9
	end if
		ccon_h2s_sat_spruce1(j+6) = (ccon_h2s_sat_spruce1(j+6) - lx_h2s_sat_temp * sqrt(lat_flux_factor1)) !* hu_h2o(j+6) 
		ccon_h2s_sat_spruce2(j) 	= (ccon_h2s_sat_spruce2(j) + lx_h2s_sat_temp * sqrt(lat_flux_factor2)) !* hu_h2o(j+6) 
end do

!write(iulog,*) "hereafter hummock 20", cdocs_sat_spruce1(1), cdocs_sat_spruce1(2), cdocs_sat_spruce1(3), cdocs_sat_spruce1(4), cdocs_sat_spruce1(5), cdocs_sat_spruce1(6), &
!cdocs_sat_spruce1(7), cdocs_sat_spruce1(8), cdocs_sat_spruce1(9), cdocs_sat_spruce1(10), cdocs_sat_spruce1(11), cdocs_sat_spruce1(12), cdocs_sat_spruce1(13), &
!cdocs_sat_spruce1(14), cdocs_sat_spruce1(15), cdocs_sat_spruce1(16), cdocs_sat_spruce1(17), cdocs_sat_spruce1(18), cdocs_sat_spruce1(19), cdocs_sat_spruce1(20)
!write(iulog,*) "hereafter holow 14", cdocs_sat_spruce2(1), cdocs_sat_spruce2(2), cdocs_sat_spruce2(3), cdocs_sat_spruce2(4), cdocs_sat_spruce2(5), cdocs_sat_spruce2(6), &
!cdocs_sat_spruce2(7), cdocs_sat_spruce2(8), cdocs_sat_spruce2(9), cdocs_sat_spruce2(10), cdocs_sat_spruce2(11), cdocs_sat_spruce2(12), cdocs_sat_spruce2(13), &
!cdocs_sat_spruce2(14)

	!~ cdocs_sat(1,1) = cdocs_sat_spruce1(1)
	!~ cdocs_sat(1,2) = cdocs_sat_spruce1(2)
	!~ cdocs_sat(1,3) = cdocs_sat_spruce1(3)
	!~ cdocs_sat(1,4) = cdocs_sat_spruce1(4)
	!~ cdocs_sat(1,5) = cdocs_sat_spruce1(5)
	!~ cdocs_sat(1,6) = cdocs_sat_spruce1(6)+cdocs_sat_spruce1(7)+cdocs_sat_spruce1(8)+cdocs_sat_spruce1(9)+cdocs_sat_spruce1(10)+cdocs_sat_spruce1(11)	!/dz(1,6)
	!~ cdocs_sat(1,7) = cdocs_sat_spruce1(12)+cdocs_sat_spruce1(13)+cdocs_sat_spruce1(14)	!/dz(1,7)
	!~ cdocs_sat(1,8) = cdocs_sat_spruce1(15)+cdocs_sat_spruce1(16)	!/dz(1,8)
	!~ cdocs_sat(1,9) = cdocs_sat_spruce1(17)+cdocs_sat_spruce1(18)	!/dz(1,9)
	!~ cdocs_sat(1,10)=cdocs_sat_spruce1(19)+cdocs_sat_spruce1(20)	!/dz(1,10)

	!~ cdocs_sat(2,1) = cdocs_sat_spruce2(1)
	!~ cdocs_sat(2,2) = cdocs_sat_spruce2(2)
	!~ cdocs_sat(2,3) = cdocs_sat_spruce2(3)
	!~ cdocs_sat(2,4) = cdocs_sat_spruce2(4)
	!~ cdocs_sat(2,5) = cdocs_sat_spruce2(5)+cdocs_sat_spruce2(6)	!/dz(2,5)
	!~ cdocs_sat(2,6) = cdocs_sat_spruce2(7)
	!~ cdocs_sat(2,7) = cdocs_sat_spruce2(8)+cdocs_sat_spruce2(9)	!/dz(2,7)
	!~ cdocs_sat(2,8) = cdocs_sat_spruce2(10)+cdocs_sat_spruce2(11)	!/dz(2,8)
	!~ cdocs_sat(2,9) = cdocs_sat_spruce2(12)+cdocs_sat_spruce2(13)	!/dz(2,9)
	!~ cdocs_sat(2,10) = cdocs_sat_spruce2(14)

	!~ cdons_sat(1,1) = cdons_sat_spruce1(1)
	!~ cdons_sat(1,2) = cdons_sat_spruce1(1)
	!~ cdons_sat(1,3) = cdons_sat_spruce1(1)
	!~ cdons_sat(1,4) = cdons_sat_spruce1(1)
	!~ cdons_sat(1,5) = cdons_sat_spruce1(1)
	!~ cdons_sat(1,6) = cdons_sat_spruce1(6)+cdons_sat_spruce1(7)+cdons_sat_spruce1(8)+cdons_sat_spruce1(9)+cdons_sat_spruce1(10)+cdons_sat_spruce1(11)	!/dz(1,6)
	!~ cdons_sat(1,7) = cdons_sat_spruce1(12)+cdons_sat_spruce1(13)+cdons_sat_spruce1(14)	!/dz(1,7)
	!~ cdons_sat(1,8) = cdons_sat_spruce1(15)+cdons_sat_spruce1(16)	!/dz(1,8)
	!~ cdons_sat(1,9) = cdons_sat_spruce1(17)+cdons_sat_spruce1(18)	!/dz(1,9)
	!~ cdons_sat(1,10)=cdons_sat_spruce1(19)+cdons_sat_spruce1(20)	!/dz(1,10)

	!~ cdons_sat(2,1) = cdons_sat_spruce2(1)
	!~ cdons_sat(2,2) = cdons_sat_spruce2(2)
	!~ cdons_sat(2,3) = cdons_sat_spruce2(3)
	!~ cdons_sat(2,4) = cdons_sat_spruce2(4)
	!~ cdons_sat(2,5) = cdons_sat_spruce2(5)+cdons_sat_spruce2(6)	!/dz(2,5)
	!~ cdons_sat(2,6) = cdons_sat_spruce2(7)
	!~ cdons_sat(2,7) = cdons_sat_spruce2(8)+cdons_sat_spruce2(9)	!/dz(2,7)
	!~ cdons_sat(2,8) = cdons_sat_spruce2(10)+cdons_sat_spruce2(11)	!/dz(2,8)
	!~ cdons_sat(2,9) = cdons_sat_spruce2(12)+cdons_sat_spruce2(13)	!/dz(2,9)
	!~ cdons_sat(2,10) =cdons_sat_spruce2(14)

	!~ caces_sat(1,1) = caces_sat_spruce1(1)
	!~ caces_sat(1,2) = caces_sat_spruce1(2)
	!~ caces_sat(1,3) = caces_sat_spruce1(3)
	!~ caces_sat(1,4) = caces_sat_spruce1(4)
	!~ caces_sat(1,5) = caces_sat_spruce1(5)
	!~ caces_sat(1,6) = caces_sat_spruce1(6)+caces_sat_spruce1(7)+caces_sat_spruce1(8)+caces_sat_spruce1(9)+caces_sat_spruce1(10)+caces_sat_spruce1(11)	!/dz(1,6)
	!~ caces_sat(1,7) = caces_sat_spruce1(12)+caces_sat_spruce1(13)+caces_sat_spruce1(14)	!/dz(1,7)
	!~ caces_sat(1,8) = caces_sat_spruce1(15)+caces_sat_spruce1(16)	!/dz(1,8)
	!~ caces_sat(1,9) = caces_sat_spruce1(17)+caces_sat_spruce1(18)	!/dz(1,9)
	!~ caces_sat(1,10)=caces_sat_spruce1(19)+caces_sat_spruce1(20)

	!~ caces_sat(2,1) = caces_sat_spruce2(1)
	!~ caces_sat(2,2) = caces_sat_spruce2(2)
	!~ caces_sat(2,3) = caces_sat_spruce2(3)
	!~ caces_sat(2,4) = caces_sat_spruce2(4)
	!~ caces_sat(2,5) = caces_sat_spruce2(5)+caces_sat_spruce2(6)	!/dz(2,5)
	!~ caces_sat(2,6) = caces_sat_spruce2(7)
	!~ caces_sat(2,7) = caces_sat_spruce2(8)+caces_sat_spruce2(9)	!/dz(2,7)
	!~ caces_sat(2,8) = caces_sat_spruce2(10)+caces_sat_spruce2(11)	!/dz(2,8)
	!~ caces_sat(2,9) = caces_sat_spruce2(12)+caces_sat_spruce2(13)	!/dz(2,9)
	!~ caces_sat(2,10)=caces_sat_spruce2(14)
	
	!~ ccon_o2s_sat(1,1) = ccon_o2s_sat_spruce1(1)
	!~ ccon_o2s_sat(1,2) = ccon_o2s_sat_spruce1(2)
	!~ ccon_o2s_sat(1,3) = ccon_o2s_sat_spruce1(3)
	!~ ccon_o2s_sat(1,4) = ccon_o2s_sat_spruce1(4)
	!~ ccon_o2s_sat(1,5) = ccon_o2s_sat_spruce1(5)
	!~ ccon_o2s_sat(1,6) = ccon_o2s_sat_spruce1(6)+ccon_o2s_sat_spruce1(7)+ccon_o2s_sat_spruce1(8)+ccon_o2s_sat_spruce1(9)+ccon_o2s_sat_spruce1(10)+ccon_o2s_sat_spruce1(11)
	!~ ccon_o2s_sat(1,7) = ccon_o2s_sat_spruce1(12)+ccon_o2s_sat_spruce1(13)+ccon_o2s_sat_spruce1(14)
	!~ ccon_o2s_sat(1,8) = ccon_o2s_sat_spruce1(15)+ccon_o2s_sat_spruce1(16)
	!~ ccon_o2s_sat(1,9) = ccon_o2s_sat_spruce1(17)+ccon_o2s_sat_spruce1(18)
	!~ ccon_o2s_sat(1,10)=ccon_o2s_sat_spruce1(19)+ccon_o2s_sat_spruce1(20)

	!~ ccon_o2s_sat(2,1) = ccon_o2s_sat_spruce2(1)
	!~ ccon_o2s_sat(2,2) = ccon_o2s_sat_spruce2(2)
	!~ ccon_o2s_sat(2,3) = ccon_o2s_sat_spruce2(3)
	!~ ccon_o2s_sat(2,4) = ccon_o2s_sat_spruce2(4)
	!~ ccon_o2s_sat(2,5) = ccon_o2s_sat_spruce2(5)+ccon_o2s_sat_spruce2(6)
	!~ ccon_o2s_sat(2,6) = ccon_o2s_sat_spruce2(7)
	!~ ccon_o2s_sat(2,7) = ccon_o2s_sat_spruce2(8)+ccon_o2s_sat_spruce2(9)
	!~ ccon_o2s_sat(2,8) = ccon_o2s_sat_spruce2(10)+ccon_o2s_sat_spruce2(11)
	!~ ccon_o2s_sat(2,9) = ccon_o2s_sat_spruce2(12)+ccon_o2s_sat_spruce2(13)
	!~ ccon_o2s_sat(2,10) =ccon_o2s_sat_spruce2(14)
	
	!~ ccon_h2s_sat(1,1) = ccon_h2s_sat_spruce1(1)
	!~ ccon_h2s_sat(1,2) = ccon_h2s_sat_spruce1(2)
	!~ ccon_h2s_sat(1,3) = ccon_h2s_sat_spruce1(3)
	!~ ccon_h2s_sat(1,4) = ccon_h2s_sat_spruce1(4)
	!~ ccon_h2s_sat(1,5) = ccon_h2s_sat_spruce1(5)
	!~ ccon_h2s_sat(1,6) = ccon_h2s_sat_spruce1(6)+ccon_h2s_sat_spruce1(7)+ccon_h2s_sat_spruce1(8)+ccon_h2s_sat_spruce1(9)+ccon_h2s_sat_spruce1(10)+ccon_h2s_sat_spruce1(11)
	!~ ccon_h2s_sat(1,7) = ccon_h2s_sat_spruce1(12)+ccon_h2s_sat_spruce1(13)+ccon_h2s_sat_spruce1(14)
	!~ ccon_h2s_sat(1,8) = ccon_h2s_sat_spruce1(15)+ccon_h2s_sat_spruce1(16)	
	!~ ccon_h2s_sat(1,9) = ccon_h2s_sat_spruce1(17)+ccon_h2s_sat_spruce1(18)
	!~ ccon_h2s_sat(1,10) = ccon_h2s_sat_spruce1(19)+ccon_h2s_sat_spruce1(20)

	!~ ccon_h2s_sat(2,1) = ccon_h2s_sat_spruce2(1)
	!~ ccon_h2s_sat(2,2) = ccon_h2s_sat_spruce2(2)
	!~ ccon_h2s_sat(2,3) = ccon_h2s_sat_spruce2(3)
	!~ ccon_h2s_sat(2,4) = ccon_h2s_sat_spruce2(4)
	!~ ccon_h2s_sat(2,5) = ccon_h2s_sat_spruce2(5)+ccon_h2s_sat_spruce2(6)
	!~ ccon_h2s_sat(2,6) = ccon_h2s_sat_spruce2(7)
	!~ ccon_h2s_sat(2,7) = ccon_h2s_sat_spruce2(8)+ccon_h2s_sat_spruce2(9)
	!~ ccon_h2s_sat(2,8) = ccon_h2s_sat_spruce2(10)+ccon_h2s_sat_spruce2(11)
	!~ ccon_h2s_sat(2,9) = ccon_h2s_sat_spruce2(12)+ccon_h2s_sat_spruce2(13)
	!~ ccon_h2s_sat(2,10) = ccon_h2s_sat_spruce2(14)
	
	!~ ccon_co2s_sat(1,1) = ccon_co2s_sat_spruce1(1)
	!~ ccon_co2s_sat(1,2) = ccon_co2s_sat_spruce1(2)
	!~ ccon_co2s_sat(1,3) = ccon_co2s_sat_spruce1(3)
	!~ ccon_co2s_sat(1,4) = ccon_co2s_sat_spruce1(4)
	!~ ccon_co2s_sat(1,5) = ccon_co2s_sat_spruce1(5)
	!~ ccon_co2s_sat(1,6) = ccon_co2s_sat_spruce1(6)+ccon_co2s_sat_spruce1(7)+ccon_co2s_sat_spruce1(8)+ccon_co2s_sat_spruce1(9)+ccon_co2s_sat_spruce1(10)+ccon_co2s_sat_spruce1(11)
	!~ ccon_co2s_sat(1,7) = ccon_co2s_sat_spruce1(12)+ccon_co2s_sat_spruce1(13)+ccon_co2s_sat_spruce1(14)
	!~ ccon_co2s_sat(1,8) = ccon_co2s_sat_spruce1(15)+ccon_co2s_sat_spruce1(16)
	!~ ccon_co2s_sat(1,9) = ccon_co2s_sat_spruce1(17)+ccon_co2s_sat_spruce1(18)
	!~ ccon_co2s_sat(1,10) = ccon_co2s_sat_spruce1(19)+ccon_co2s_sat_spruce1(20)
	
	!~ ccon_co2s_sat(2,1) = ccon_co2s_sat_spruce2(1)
	!~ ccon_co2s_sat(2,2) = ccon_co2s_sat_spruce2(2)
	!~ ccon_co2s_sat(2,3) = ccon_co2s_sat_spruce2(3)
	!~ ccon_co2s_sat(2,4) = ccon_co2s_sat_spruce2(4)
	!~ ccon_co2s_sat(2,5) = ccon_co2s_sat_spruce2(5)+ccon_co2s_sat_spruce2(6)
	!~ ccon_co2s_sat(2,6) = ccon_co2s_sat_spruce2(7)
	!~ ccon_co2s_sat(2,7) = ccon_co2s_sat_spruce2(8)+ccon_co2s_sat_spruce2(9)
	!~ ccon_co2s_sat(2,8) = ccon_co2s_sat_spruce2(10)+ccon_co2s_sat_spruce2(11)
	!~ ccon_co2s_sat(2,9) = ccon_co2s_sat_spruce2(12)+ccon_co2s_sat_spruce2(13)
	!~ ccon_co2s_sat(2,10) = ccon_co2s_sat_spruce2(14)
	
	!~ ccon_ch4s_sat(1,1) = ccon_ch4s_sat_spruce1(1)
	!~ ccon_ch4s_sat(1,2) = ccon_ch4s_sat_spruce1(2)
	!~ ccon_ch4s_sat(1,3) = ccon_ch4s_sat_spruce1(3)
	!~ ccon_ch4s_sat(1,4) = ccon_ch4s_sat_spruce1(4)
	!~ ccon_ch4s_sat(1,5) = ccon_ch4s_sat_spruce1(5)
	!~ ccon_ch4s_sat(1,6) = ccon_ch4s_sat_spruce1(6)+ccon_ch4s_sat_spruce1(7)+ccon_ch4s_sat_spruce1(8)+ccon_ch4s_sat_spruce1(9)+ccon_ch4s_sat_spruce1(10)+ccon_ch4s_sat_spruce1(11)
	!~ ccon_ch4s_sat(1,7) = ccon_ch4s_sat_spruce1(12)+ccon_ch4s_sat_spruce1(13)+ccon_ch4s_sat_spruce1(14)
	!~ ccon_ch4s_sat(1,8) = ccon_ch4s_sat_spruce1(15)+ccon_ch4s_sat_spruce1(16)
	!~ ccon_ch4s_sat(1,9) = ccon_ch4s_sat_spruce1(17)+ccon_ch4s_sat_spruce1(18)
	!~ ccon_ch4s_sat(1,10) = ccon_ch4s_sat_spruce1(19)+ccon_ch4s_sat_spruce1(20)

	!~ ccon_ch4s_sat(2,1) = ccon_ch4s_sat_spruce2(1)
	!~ ccon_ch4s_sat(2,2) = ccon_ch4s_sat_spruce2(2)
	!~ ccon_ch4s_sat(2,3) = ccon_ch4s_sat_spruce2(3)
	!~ ccon_ch4s_sat(2,4) = ccon_ch4s_sat_spruce2(4)
	!~ ccon_ch4s_sat(2,5) = ccon_ch4s_sat_spruce2(5)+ccon_ch4s_sat_spruce2(6)
	!~ ccon_ch4s_sat(2,6) = ccon_ch4s_sat_spruce2(7)
	!~ ccon_ch4s_sat(2,7) = ccon_ch4s_sat_spruce2(8)+ccon_ch4s_sat_spruce2(9)
	!~ ccon_ch4s_sat(2,8) = ccon_ch4s_sat_spruce2(10)+ccon_ch4s_sat_spruce2(11)
	!~ ccon_ch4s_sat(2,9) = ccon_ch4s_sat_spruce2(12)+ccon_ch4s_sat_spruce2(13)
	!~ ccon_ch4s_sat(2,10) = ccon_ch4s_sat_spruce2(14)	
	
!~ ! aggregate the 20/14 layers bgc variables to 10 layers
	cdocs_sat(1,1) = cdocs_sat_spruce1(1)
	cdocs_sat(1,2) = cdocs_sat_spruce1(2)
	cdocs_sat(1,3) = cdocs_sat_spruce1(3)
	cdocs_sat(1,4) = cdocs_sat_spruce1(4)
	cdocs_sat(1,5) = cdocs_sat_spruce1(5)
	cdocs_sat(1,6) = (cdocs_sat_spruce1(6)*hu_soil_tk(6)+cdocs_sat_spruce1(7)*hu_soil_tk(7)+cdocs_sat_spruce1(8)*hu_soil_tk(8)+cdocs_sat_spruce1(9)*hu_soil_tk(9)+&
	cdocs_sat_spruce1(10)*hu_soil_tk(10)+cdocs_sat_spruce1(11)*hu_soil_tk(11)) / dz(1,6)
	
	cdocs_sat(1,7) = (cdocs_sat_spruce1(12)*hu_soil_tk(12)+cdocs_sat_spruce1(13)*hu_soil_tk(13)+cdocs_sat_spruce1(14)*hu_soil_tk(14)) / dz(1,7)
	cdocs_sat(1,8) = (cdocs_sat_spruce1(15)*hu_soil_tk(15)+cdocs_sat_spruce1(16)*hu_soil_tk(16)) / dz(1,8)
	cdocs_sat(1,9) = (cdocs_sat_spruce1(17)*hu_soil_tk(17)+cdocs_sat_spruce1(18)*hu_soil_tk(18)) / dz(1,9)
	cdocs_sat(1,10) = (cdocs_sat_spruce1(19)*hu_soil_tk(19)+cdocs_sat_spruce1(20)*hu_soil_tk(20)) / dz(1,10)
!~ write(iulog,*) hu_soil_tk(6)+hu_soil_tk(7)+hu_soil_tk(8)+hu_soil_tk(9)+hu_soil_tk(10)+hu_soil_tk(11)
!~ write(iulog,*) dz(1,6)
!~ write(iulog,*) hu_soil_tk(12)+hu_soil_tk(13)+hu_soil_tk(14)
!~ write(iulog,*) dz(1,7)
!~ write(iulog,*) hu_soil_tk(15)+hu_soil_tk(16)
!~ write(iulog,*) dz(1,8)
!~ write(iulog,*) hu_soil_tk(17)+hu_soil_tk(18)
!~ write(iulog,*) dz(1,9)
!~ write(iulog,*) hu_soil_tk(19)+hu_soil_tk(20)
!~ write(iulog,*) dz(1,10)

	cdocs_sat(2,1) = cdocs_sat_spruce2(1)
	cdocs_sat(2,2) = cdocs_sat_spruce2(2)
	cdocs_sat(2,3) = cdocs_sat_spruce2(3)
	cdocs_sat(2,4) = cdocs_sat_spruce2(4)
	cdocs_sat(2,5) = (cdocs_sat_spruce2(5)*ho_soil_tk(5)+cdocs_sat_spruce2(6)*ho_soil_tk(6)) / dz(2,5)
	cdocs_sat(2,6) = cdocs_sat_spruce2(7)
	cdocs_sat(2,7) = (cdocs_sat_spruce2(8)*ho_soil_tk(8)+cdocs_sat_spruce2(9)*ho_soil_tk(9)) / dz(2,7)
	cdocs_sat(2,8) = (cdocs_sat_spruce2(10)*ho_soil_tk(10)+cdocs_sat_spruce2(11)*ho_soil_tk(11)) / dz(2,8)
	cdocs_sat(2,9) = (cdocs_sat_spruce2(12)*ho_soil_tk(12)+cdocs_sat_spruce2(13)*ho_soil_tk(13)) / dz(2,9)
	cdocs_sat(2,10) = cdocs_sat_spruce2(14)*ho_soil_tk(14) / dz(2,10)

	cdons_sat(1,1) = cdons_sat_spruce1(1)
	cdons_sat(1,2) = cdons_sat_spruce1(1)
	cdons_sat(1,3) = cdons_sat_spruce1(1)
	cdons_sat(1,4) = cdons_sat_spruce1(1)
	cdons_sat(1,5) = cdons_sat_spruce1(1)
	cdons_sat(1,6) = (cdons_sat_spruce1(6)*hu_soil_tk(6)+cdons_sat_spruce1(7)*hu_soil_tk(7)+cdons_sat_spruce1(8)*hu_soil_tk(8)+cdons_sat_spruce1(9)*hu_soil_tk(9)+&
	cdons_sat_spruce1(10)*hu_soil_tk(10)+cdons_sat_spruce1(11)*hu_soil_tk(11)) / dz(1,6)
	
	cdons_sat(1,7) = (cdons_sat_spruce1(12)*hu_soil_tk(12)+cdons_sat_spruce1(13)*hu_soil_tk(13)+cdons_sat_spruce1(14)*hu_soil_tk(14)) / dz(1,7)
	cdons_sat(1,8) = (cdons_sat_spruce1(15)*hu_soil_tk(15)+cdons_sat_spruce1(16)*hu_soil_tk(16)) / dz(1,8)
	cdons_sat(1,9) = (cdons_sat_spruce1(17)*hu_soil_tk(17)+cdons_sat_spruce1(18)*hu_soil_tk(18)) / dz(1,9)
	cdons_sat(1,10) = (cdons_sat_spruce1(19)*hu_soil_tk(19)+cdons_sat_spruce1(20)*hu_soil_tk(20)) / dz(1,10)

	cdons_sat(2,1) = cdons_sat_spruce2(1)
	cdons_sat(2,2) = cdons_sat_spruce2(2)
	cdons_sat(2,3) = cdons_sat_spruce2(3)
	cdons_sat(2,4) = cdons_sat_spruce2(4)
	cdons_sat(2,5) = (cdons_sat_spruce2(5)*ho_soil_tk(5)+cdons_sat_spruce2(6)*ho_soil_tk(6)) / dz(2,5)
	cdons_sat(2,6) = cdons_sat_spruce2(7)
	cdons_sat(2,7) = (cdons_sat_spruce2(8)*ho_soil_tk(8)+cdons_sat_spruce2(9)*ho_soil_tk(9)) / dz(2,7)
	cdons_sat(2,8) = (cdons_sat_spruce2(10)*ho_soil_tk(10)+cdons_sat_spruce2(11)*ho_soil_tk(11)) / dz(2,8)
	cdons_sat(2,9) = (cdons_sat_spruce2(12)*ho_soil_tk(12)+cdons_sat_spruce2(13)*ho_soil_tk(13)) / dz(2,9)
	cdons_sat(2,10) = cdons_sat_spruce2(14)*ho_soil_tk(14) / dz(2,10)

	caces_sat(1,1) = caces_sat_spruce1(1)
	caces_sat(1,2) = caces_sat_spruce1(2)
	caces_sat(1,3) = caces_sat_spruce1(3)
	caces_sat(1,4) = caces_sat_spruce1(4)
	caces_sat(1,5) = caces_sat_spruce1(5)
	caces_sat(1,6) = (caces_sat_spruce1(6)*hu_soil_tk(6)+caces_sat_spruce1(7)*hu_soil_tk(7)+caces_sat_spruce1(8)*hu_soil_tk(8)+caces_sat_spruce1(9)*hu_soil_tk(9)+&
	caces_sat_spruce1(10)*hu_soil_tk(10)+caces_sat_spruce1(11)*hu_soil_tk(11)) / dz(1,6)
	
	caces_sat(1,7) = (caces_sat_spruce1(12)*hu_soil_tk(12)+caces_sat_spruce1(13)*hu_soil_tk(13)+caces_sat_spruce1(14)*hu_soil_tk(14))	/ dz(1,7)
	caces_sat(1,8) = (caces_sat_spruce1(15)*hu_soil_tk(15)+caces_sat_spruce1(16)*hu_soil_tk(16)) / dz(1,8)
	caces_sat(1,9) = (caces_sat_spruce1(17)*hu_soil_tk(17)+caces_sat_spruce1(18)*hu_soil_tk(18)) / dz(1,9)
	caces_sat(1,10) = (caces_sat_spruce1(19)*hu_soil_tk(19)+caces_sat_spruce1(20)*hu_soil_tk(20)) / dz(1,10)

	caces_sat(2,1) = caces_sat_spruce2(1)
	caces_sat(2,2) = caces_sat_spruce2(2)
	caces_sat(2,3) = caces_sat_spruce2(3)
	caces_sat(2,4) = caces_sat_spruce2(4)
	caces_sat(2,5) = (caces_sat_spruce2(5)*ho_soil_tk(5)+caces_sat_spruce2(6)*ho_soil_tk(6)) / dz(2,5)
	caces_sat(2,6) = caces_sat_spruce2(7)
	caces_sat(2,7) = (caces_sat_spruce2(8)*ho_soil_tk(8)+caces_sat_spruce2(9)*ho_soil_tk(9)) / dz(2,7)
	caces_sat(2,8) = (caces_sat_spruce2(10)*ho_soil_tk(10)+caces_sat_spruce2(11)*ho_soil_tk(11)) / dz(2,8)
	caces_sat(2,9) = (caces_sat_spruce2(12)*ho_soil_tk(12)+caces_sat_spruce2(13)*ho_soil_tk(13)) / dz(2,9)
	caces_sat(2,10) = caces_sat_spruce2(14)*ho_soil_tk(14) / dz(2,10)
	
	ccon_o2s_sat(1,1) = ccon_o2s_sat_spruce1(1)
	ccon_o2s_sat(1,2) = ccon_o2s_sat_spruce1(2)
	ccon_o2s_sat(1,3) = ccon_o2s_sat_spruce1(3)
	ccon_o2s_sat(1,4) = ccon_o2s_sat_spruce1(4)
	ccon_o2s_sat(1,5) = ccon_o2s_sat_spruce1(5)
	ccon_o2s_sat(1,6) = (ccon_o2s_sat_spruce1(6)*hu_soil_tk(6)+ccon_o2s_sat_spruce1(7)*hu_soil_tk(7)+ccon_o2s_sat_spruce1(8)*hu_soil_tk(8)+ccon_o2s_sat_spruce1(9)*hu_soil_tk(9)+&
	ccon_o2s_sat_spruce1(10)*hu_soil_tk(10)+ccon_o2s_sat_spruce1(11)*hu_soil_tk(11)) / dz(1,6)
	
	ccon_o2s_sat(1,7) = (ccon_o2s_sat_spruce1(12)*hu_soil_tk(12)+ccon_o2s_sat_spruce1(13)*hu_soil_tk(13)+ccon_o2s_sat_spruce1(14)*hu_soil_tk(14)) / dz(1,7)
	ccon_o2s_sat(1,8) = (ccon_o2s_sat_spruce1(15)*hu_soil_tk(15)+ccon_o2s_sat_spruce1(16)*hu_soil_tk(16)) / dz(1,8)
	ccon_o2s_sat(1,9) = (ccon_o2s_sat_spruce1(17)*hu_soil_tk(17)+ccon_o2s_sat_spruce1(18)*hu_soil_tk(18)) / dz(1,9)
	ccon_o2s_sat(1,10) = (ccon_o2s_sat_spruce1(19)*hu_soil_tk(19)+ccon_o2s_sat_spruce1(20)*hu_soil_tk(20)) / dz(1,10)

	ccon_o2s_sat(2,1) = ccon_o2s_sat_spruce2(1)
	ccon_o2s_sat(2,2) = ccon_o2s_sat_spruce2(2)
	ccon_o2s_sat(2,3) = ccon_o2s_sat_spruce2(3)
	ccon_o2s_sat(2,4) = ccon_o2s_sat_spruce2(4)
	ccon_o2s_sat(2,5) = (ccon_o2s_sat_spruce2(5)*ho_soil_tk(5)+ccon_o2s_sat_spruce2(6)*ho_soil_tk(6)) / dz(2,5)
	ccon_o2s_sat(2,6) = ccon_o2s_sat_spruce2(7)
	ccon_o2s_sat(2,7) = (ccon_o2s_sat_spruce2(8)*ho_soil_tk(8)+ccon_o2s_sat_spruce2(9)*ho_soil_tk(9)) / dz(2,7)
	ccon_o2s_sat(2,8) = (ccon_o2s_sat_spruce2(10)*ho_soil_tk(10)+ccon_o2s_sat_spruce2(11)*ho_soil_tk(11)) / dz(2,8)
	ccon_o2s_sat(2,9) = (ccon_o2s_sat_spruce2(12)*ho_soil_tk(12)+ccon_o2s_sat_spruce2(13)*ho_soil_tk(13)) / dz(2,9)
	ccon_o2s_sat(2,10) = ccon_o2s_sat_spruce2(14)*ho_soil_tk(14) / dz(2,10)
	
	ccon_h2s_sat(1,1) = ccon_h2s_sat_spruce1(1)
	ccon_h2s_sat(1,2) = ccon_h2s_sat_spruce1(2)
	ccon_h2s_sat(1,3) = ccon_h2s_sat_spruce1(3)
	ccon_h2s_sat(1,4) = ccon_h2s_sat_spruce1(4)
	ccon_h2s_sat(1,5) = ccon_h2s_sat_spruce1(5)
	ccon_h2s_sat(1,6) = (ccon_h2s_sat_spruce1(6)*hu_soil_tk(6)+ccon_h2s_sat_spruce1(7)*hu_soil_tk(7)+ccon_h2s_sat_spruce1(8)*hu_soil_tk(8)+ccon_h2s_sat_spruce1(9)*hu_soil_tk(9)+&
	ccon_h2s_sat_spruce1(10)*hu_soil_tk(10)+ccon_h2s_sat_spruce1(11)*hu_soil_tk(11)) / dz(1,6)
	
	ccon_h2s_sat(1,7) = (ccon_h2s_sat_spruce1(12)*hu_soil_tk(12)+ccon_h2s_sat_spruce1(13)*hu_soil_tk(13)+ccon_h2s_sat_spruce1(14)*hu_soil_tk(14)) / dz(1,7)
	ccon_h2s_sat(1,8) = (ccon_h2s_sat_spruce1(15)*hu_soil_tk(15)+ccon_h2s_sat_spruce1(16)*hu_soil_tk(16)) / dz(1,8)
	ccon_h2s_sat(1,9) = (ccon_h2s_sat_spruce1(17)*hu_soil_tk(17)+ccon_h2s_sat_spruce1(18)*hu_soil_tk(18)) / dz(1,9)
	ccon_h2s_sat(1,10) = (ccon_h2s_sat_spruce1(19)*hu_soil_tk(19)+ccon_h2s_sat_spruce1(20)*hu_soil_tk(20)) / dz(1,10)

	ccon_h2s_sat(2,1) = ccon_h2s_sat_spruce2(1)
	ccon_h2s_sat(2,2) = ccon_h2s_sat_spruce2(2)
	ccon_h2s_sat(2,3) = ccon_h2s_sat_spruce2(3)
	ccon_h2s_sat(2,4) = ccon_h2s_sat_spruce2(4)
	ccon_h2s_sat(2,5) = (ccon_h2s_sat_spruce2(5)*ho_soil_tk(5)+ccon_h2s_sat_spruce2(6)*ho_soil_tk(6)) / dz(2,5)
	ccon_h2s_sat(2,6) = ccon_h2s_sat_spruce2(7)
	ccon_h2s_sat(2,7) = (ccon_h2s_sat_spruce2(8)*ho_soil_tk(8)+ccon_h2s_sat_spruce2(9)*ho_soil_tk(9)) / dz(2,7)
	ccon_h2s_sat(2,8) = (ccon_h2s_sat_spruce2(10)*ho_soil_tk(10)+ccon_h2s_sat_spruce2(11)*ho_soil_tk(11)) / dz(2,8)
	ccon_h2s_sat(2,9) = (ccon_h2s_sat_spruce2(12)*ho_soil_tk(12)+ccon_h2s_sat_spruce2(13)*ho_soil_tk(13)) / dz(2,9)
	ccon_h2s_sat(2,10) = ccon_h2s_sat_spruce2(14)*ho_soil_tk(14) / dz(2,10)
	
	ccon_co2s_sat(1,1) = ccon_co2s_sat_spruce1(1)
	ccon_co2s_sat(1,2) = ccon_co2s_sat_spruce1(2)
	ccon_co2s_sat(1,3) = ccon_co2s_sat_spruce1(3)
	ccon_co2s_sat(1,4) = ccon_co2s_sat_spruce1(4)
	ccon_co2s_sat(1,5) = ccon_co2s_sat_spruce1(5)
	ccon_co2s_sat(1,6) = (ccon_co2s_sat_spruce1(6)*hu_soil_tk(6)+ccon_co2s_sat_spruce1(7)*hu_soil_tk(7)+ccon_co2s_sat_spruce1(8)*hu_soil_tk(8)+ccon_co2s_sat_spruce1(9)*hu_soil_tk(9)+&
	ccon_co2s_sat_spruce1(10)*hu_soil_tk(10)+ccon_co2s_sat_spruce1(11)*hu_soil_tk(11)) / dz(1,6)
	
	ccon_co2s_sat(1,7) = (ccon_co2s_sat_spruce1(12)*hu_soil_tk(12)+ccon_co2s_sat_spruce1(13)*hu_soil_tk(13)+ccon_co2s_sat_spruce1(14)*hu_soil_tk(14)) / dz(1,7)
	ccon_co2s_sat(1,8) = (ccon_co2s_sat_spruce1(15)*hu_soil_tk(15)+ccon_co2s_sat_spruce1(16)*hu_soil_tk(16)) / dz(1,8)
	ccon_co2s_sat(1,9) = (ccon_co2s_sat_spruce1(17)*hu_soil_tk(17)+ccon_co2s_sat_spruce1(18)*hu_soil_tk(18)) / dz(1,9)
	ccon_co2s_sat(1,10) = (ccon_co2s_sat_spruce1(19)*hu_soil_tk(19)+ccon_co2s_sat_spruce1(20)*hu_soil_tk(20)) / dz(1,10)

	ccon_co2s_sat(2,1) = ccon_co2s_sat_spruce2(1)
	ccon_co2s_sat(2,2) = ccon_co2s_sat_spruce2(2)
	ccon_co2s_sat(2,3) = ccon_co2s_sat_spruce2(3)
	ccon_co2s_sat(2,4) = ccon_co2s_sat_spruce2(4)
	ccon_co2s_sat(2,5) = (ccon_co2s_sat_spruce2(5)*ho_soil_tk(5)+ccon_co2s_sat_spruce2(6)*ho_soil_tk(6)) / dz(2,5)
	ccon_co2s_sat(2,6) = ccon_co2s_sat_spruce2(7)
	ccon_co2s_sat(2,7) = (ccon_co2s_sat_spruce2(8)*ho_soil_tk(8)+ccon_co2s_sat_spruce2(9)*ho_soil_tk(9)) / dz(2,7)
	ccon_co2s_sat(2,8) = (ccon_co2s_sat_spruce2(10)*ho_soil_tk(10)+ccon_co2s_sat_spruce2(11)*ho_soil_tk(11)) / dz(2,8)
	ccon_co2s_sat(2,9) = (ccon_co2s_sat_spruce2(12)*ho_soil_tk(12)+ccon_co2s_sat_spruce2(13)*ho_soil_tk(13)) / dz(2,9)
	ccon_co2s_sat(2,10) = ccon_co2s_sat_spruce2(14)*ho_soil_tk(14) / dz(2,10)
	
	ccon_ch4s_sat(1,1) = ccon_ch4s_sat_spruce1(1)
	ccon_ch4s_sat(1,2) = ccon_ch4s_sat_spruce1(2)
	ccon_ch4s_sat(1,3) = ccon_ch4s_sat_spruce1(3)
	ccon_ch4s_sat(1,4) = ccon_ch4s_sat_spruce1(4)
	ccon_ch4s_sat(1,5) = ccon_ch4s_sat_spruce1(5)
	ccon_ch4s_sat(1,6) = (ccon_ch4s_sat_spruce1(6)*hu_soil_tk(6)+ccon_ch4s_sat_spruce1(7)*hu_soil_tk(7)+ccon_ch4s_sat_spruce1(8)*hu_soil_tk(8)+ccon_ch4s_sat_spruce1(9)*hu_soil_tk(9)+&
	ccon_ch4s_sat_spruce1(10)*hu_soil_tk(10)+ccon_ch4s_sat_spruce1(11)*hu_soil_tk(11)) / dz(1,6)
	
	ccon_ch4s_sat(1,7) = (ccon_ch4s_sat_spruce1(12)*hu_soil_tk(12)+ccon_ch4s_sat_spruce1(13)*hu_soil_tk(13)+ccon_ch4s_sat_spruce1(14)*hu_soil_tk(14)) / dz(1,7)
	ccon_ch4s_sat(1,8) = (ccon_ch4s_sat_spruce1(15)*hu_soil_tk(15)+ccon_ch4s_sat_spruce1(16)*hu_soil_tk(16)) / dz(1,8)
	ccon_ch4s_sat(1,9) = (ccon_ch4s_sat_spruce1(17)*hu_soil_tk(17)+ccon_ch4s_sat_spruce1(18)*hu_soil_tk(18)) / dz(1,9)
	ccon_ch4s_sat(1,10) = (ccon_ch4s_sat_spruce1(19)*hu_soil_tk(19)+ccon_ch4s_sat_spruce1(20)*hu_soil_tk(20)) / dz(1,10)

	ccon_ch4s_sat(2,1) = ccon_ch4s_sat_spruce2(1)
	ccon_ch4s_sat(2,2) = ccon_ch4s_sat_spruce2(2)
	ccon_ch4s_sat(2,3) = ccon_ch4s_sat_spruce2(3)
	ccon_ch4s_sat(2,4) = ccon_ch4s_sat_spruce2(4)
	ccon_ch4s_sat(2,5) = (ccon_ch4s_sat_spruce2(5)*ho_soil_tk(5)+ccon_ch4s_sat_spruce2(6)*ho_soil_tk(6)) / dz(2,5)
	ccon_ch4s_sat(2,6) = ccon_ch4s_sat_spruce2(7)
	ccon_ch4s_sat(2,7) = (ccon_ch4s_sat_spruce2(8)*ho_soil_tk(8)+ccon_ch4s_sat_spruce2(9)*ho_soil_tk(9)) / dz(2,7)
	ccon_ch4s_sat(2,8) = (ccon_ch4s_sat_spruce2(10)*ho_soil_tk(10)+ccon_ch4s_sat_spruce2(11)*ho_soil_tk(11)) / dz(2,8)
	ccon_ch4s_sat(2,9) = (ccon_ch4s_sat_spruce2(12)*ho_soil_tk(12)+ccon_ch4s_sat_spruce2(13)*ho_soil_tk(13)) / dz(2,9)
	ccon_ch4s_sat(2,10) = ccon_ch4s_sat_spruce2(14)*ho_soil_tk(14) / dz(2,10)	
	
!write(iulog,*) "hereafter hummock 10", cdocs_sat(1,1), cdocs_sat(1,2), cdocs_sat(1,3), cdocs_sat(1,4), cdocs_sat(1,5), cdocs_sat(1,6), cdocs_sat(1,7), cdocs_sat(1,8), cdocs_sat(1,9), cdocs_sat(1,10)
!write(iulog,*) "hereafter hollow 10", cdocs_sat(2,1), cdocs_sat(2,2), cdocs_sat(2,3), cdocs_sat(2,4), cdocs_sat(2,5), cdocs_sat(2,6), cdocs_sat(2,7), cdocs_sat(2,8), cdocs_sat(2,9), cdocs_sat(2,10)
end subroutine lateral_bgc
#endif

!#endif
! for the model testing

#endif
end module microbeMod

