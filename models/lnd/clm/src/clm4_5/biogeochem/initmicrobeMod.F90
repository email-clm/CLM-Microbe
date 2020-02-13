module initmicrobeMod
#ifdef MICROBE

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: initmicrobeMod
!
! !DESCRIPTION:
! Contains time constant (and flux / diagnostic vars) and time-varying (state vars) initialization code for CH4 scheme.
! 
!
! !PUBLIC TYPES:
  implicit none
  save
  private
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: initmicrobe ! driver
!
! !PRIVATE MEMBER FUNCTIONS:
  private :: initTimeConst_microbe        ! Set constant parameters.
  private :: makearbinit_microbe          ! Set time-variable parameters for spin up.
!
! !REVISION HISTORY:
! Created by Xiaofeng Xu, 2013.
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: initch4
!
! !INTERFACE:
  subroutine initmicrobe( arbinit )
!
! !DESCRIPTION:
! Calls initTimeConst_microbe.
! Calls makearbinit_microbe with logical arbinit. If arbinit == .true. OR if initial conditions file
! does not contain methane and oxygen concentrations, then initializes time varying values. This
! allows back-compatibility with initial condition files that have not been spun up with the new
! lake code. In future versions, this could be phased out.
! 
! !ARGUMENTS:
	implicit none
	logical, intent(in) ::  arbinit ! Whether mkarbinit has been called.
!
! !CALLED FROM:
! subroutine initialize2 in module initializeMod
!
! !REVISION HISTORY:
! Created by Xiaofeng Xu, 2013
!
! !LOCAL VARIABLES:
!
!
!EOP
!

	call initTimeConst_microbe()
! Attn EK
! For now
	call makearbinit_microbe(arbinit)
! For future versions always using initial condition files spun up with the new microbe code:

end subroutine initmicrobe

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: makearbinit_ch4
!
! !INTERFACE:
  subroutine makearbinit_microbe( arbinit )
!
! !DESCRIPTION:
! If arbinit == .true., or if methane & oxygen concentrations (or lake soil org matter) 
! have not been initialized, then sets time
! varying values.
! Initializes the following time varying variables:
! 
! USES:
	use shr_kind_mod , only : r8 => shr_kind_r8
	use clmtype
	use clm_varpar   , only : nlevsoi, nlevgrnd, ndecomp_cascade_transitions
	use clm_varcon   , only : istsoil, istdlak, spval, istcrop
	use spmdMod      , only : masterproc
	use decompMod , only : get_proc_bounds
	use clm_varctl , only : iulog
!
! !ARGUMENTS:
	implicit none
	logical, intent(in) :: arbinit ! Whether mkarbinit has been called.
!
! !CALLED FROM:
! subroutine initialize in module initializeMod
!
! !REVISION HISTORY:
! Created by Zack Subin, 2009.
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
	integer , pointer :: clandunit(:)      			! landunit index associated with each column
	integer , pointer :: ltype(:)          				! landunit type
!real(r8), pointer :: cellorg(:,:)      							! column 3D organic matter (kg/m^3, 58% by mass carbon) (nlevsoi)
	real(r8), pointer :: fsat_pre(:)       				! finundated from previous timestep
	real(r8), pointer :: waterhead_unsat(:)       			
!
! local pointers to implicit out arguments
!
	real(r8), pointer :: cmicbiocs(:,:)				! column-level biomass of all microbe carbon
	real(r8), pointer :: cdocs_pre(:,:)					! column-level concentration of DOC
	real(r8), pointer :: cdocs(:,:)					! column-level concentration of DOC
	real(r8), pointer :: cdocs_unsat(:,:)			! column-level concentration of DOC in unsaturated fraction
	real(r8), pointer :: cdocs_sat(:,:)				! column-level concentration of DOC in saturated fraction
	real(r8), pointer :: cmicbions(:,:)				! column-level biomass of all microbes nitrogen
	real(r8), pointer :: cdons(:,:)					! column-level concentration of DON
	real(r8), pointer :: cdons_unsat(:,:)			! column-level concentration of DON in unsaturated fraction
	real(r8), pointer :: cdons_sat(:,:)		
	real(r8), pointer :: cdons_min(:,:)		
	real(r8), pointer :: caces(:,:)					! column-level concentration of acetate
	real(r8), pointer :: cacebios(:,:)				! column-level biomass of methanogen based on acetate
	real(r8), pointer :: cco2bios(:,:)				! column-level biomass of methanogen based on CO2/H2
	real(r8), pointer :: caerch4bios(:,:)			! column-level biomass of aerobix methanotrophy
	real(r8), pointer :: canaerch4bios(:,:)			! column-level biomass of anaerobic methanotrophy
	real(r8), pointer :: cacebios_unsat(:,:)			! column-level biomass of methanogen from acetate in unsaturated fraction
	real(r8), pointer :: cacebios_sat(:,:)			! column-level biomass of methanogen from acetate in saturated fraction
	real(r8), pointer :: cco2bios_unsat(:,:)			! column-level biomass of methanogen from CO2 and H2 in unsaturated fraction
	real(r8), pointer :: cco2bios_sat(:,:)			! column-level biomass of methanogen from CO2 and H2 in saturated fraction
	real(r8), pointer :: caerch4bios_unsat(:,:)		! column-level biomass of aerobic methanotrophy in unsaturated fraction
	real(r8), pointer :: caerch4bios_sat(:,:)			! column-level biomass of aerobic methanotrophy in saturated fraction
	real(r8), pointer :: canaerch4bios_unsat(:,:)		! column-level biomass of anaerobic methanotrophy in unsaturated fraction
	real(r8), pointer :: canaerch4bios_sat(:,:)		! column-level biomass of anaerobic methanotrophy in saturated fraction
	real(r8), pointer :: caces_unsat(:,:)			! column-level acetate in unsaturated fraction      
	real(r8), pointer :: caces_sat(:,:)				! column-level acetate in saturated fraction        

	real(r8), pointer :: ccon_ch4s(:,:)				! column-level concentration of CH4
	real(r8), pointer :: ccon_o2s(:,:)				! column-level concentration of O2
	real(r8), pointer :: ccon_co2s(:,:)				! column-level concentration of CO2
	real(r8), pointer :: ccon_h2s(:,:)				! column-level concentrtation of H2  
	real(r8), pointer :: ccon_ch4s_unsat(:,:)		! column-level concentration of CH4 in unsaturated fraction
	real(r8), pointer :: ccon_ch4s_sat(:,:)			! column-level concentration of CH4 in saturated fraction
	real(r8), pointer :: ccon_o2s_unsat(:,:)			! column-level concentration of O2 in unsaturated fraction
	real(r8), pointer :: ccon_o2s_sat(:,:)			! column-level concentration of O2 in saturated fraction
	real(r8), pointer :: ccon_co2s_unsat(:,:)			! column-level concentration of CO2 in unsaturated fraction
	real(r8), pointer :: ccon_co2s_sat(:,:)			! column-level concentration of CO2 in saturated fraction
	real(r8), pointer :: ccon_h2s_unsat(:,:)			! column-level concentration of H2 in unsaturated fraction
	real(r8), pointer :: ccon_h2s_sat(:,:)			! column-level concentration of H2 in saturated fraction
    
	real(r8), pointer :: micbio_hr_vr(:,:,:)       			! vertically resoluved microbial respiration 
	real(r8), pointer :: micbio_hr(:,:)       			! vertically resoluved microbial respiration 
	real(r8), pointer :: micbiohr_col(:)				! (gC/m2) total column microbial respiration
	real(r8), pointer :: micbioc_col(:)				! (gC/m2) total column carbon in microbial biomass 
	real(r8), pointer :: doc_col(:)					! (gC/m2) total column dissolved organic carbon
	real(r8), pointer :: dochr_col(:)				! (gC/m2) total column decomposition of dissolved organic carbon 
	real(r8), pointer :: dochr_vr(:,:)       			! vertically resoluved dissolved organic carbon
	real(r8), pointer :: ace_col(:)					! (gC/m2) total column carbon in acetate acid 
	real(r8), pointer :: ace_prod_col(:)			! (gC/m2/s) productionrate of total column carbon in acetate acid 
	real(r8), pointer :: acebios_col(:)				! (gC/m2) total column microbial biomass carbon in methanogen based on acetate 
	real(r8), pointer :: co2bio_col(:)				! (gC/m2) total column microbial biomass carbon in methanogen based on CO2/H2 
	real(r8), pointer :: aerch4bio_col(:)			! (gC/m2) total column microbial biomass carbon in aerobic methanotrophy 
	real(r8), pointer :: anaerch4bio_col(:)			! (gC/m2) total column microbial biomass carbon in anaerobic methanotrophy
	real(r8), pointer :: ccon_ch4s_col(:)				! column-level concentration of CH4
	real(r8), pointer :: ccon_co2s_col(:)				! column-level concentration of CO2
	
	real(r8), Pointer :: froot_r(:,:)
	real(r8), Pointer :: root2doc(:,:)				! part of fine root exudate carbon to soil to develop dissolved organic carbon
											! we set it as part of fine root respiration
        real(r8), Pointer :: cn_microbe(:,:)
	real(r8), pointer :: caces_unsat_prod(:,:)			! column-level acetate in unsaturated fraction      
	real(r8), pointer :: caces_sat_prod(:,:)				! column-level acetate in saturated fraction        
	real(r8), pointer :: caces_prod(:,:)				! column-level acetate in saturated fraction        

	real(r8), pointer :: caces_unsat_prod_h2(:,:)			! column-level acetogenesis in unsaturated fraction      
	real(r8), pointer :: caces_sat_prod_h2(:,:)				! column-level acetogenesis in saturated fraction        
	real(r8), pointer :: caces_prod_h2(:,:)				! column-level acetogenesis in saturated fraction        
   
   	real(r8), pointer :: micbion_col(:)				! (gN/m2) total column nitrogen in microbial biomass 
	real(r8), pointer :: don_col(:)					! (gN/m2) total column dissolved organic nitrogen
	
    !real(r8), pointer :: qflx_surf_lag(:)        				! time-lagged surface runoff (mm H2O /s)
    !real(r8), pointer :: finundated_lag(:)       				! time-lagged fractional inundated area
    !real(r8), pointer :: o2stress_unsat(:,:)     				! Ratio of oxygen available to that demanded by roots, aerobes, & methanotrophs (nlevsoi)
    !real(r8), pointer :: o2stress_sat(:,:)       				! Ratio of oxygen available to that demanded by roots, aerobes, & methanotrophs (nlevsoi)
	real(r8), pointer :: finundated(:)           		! inundated gridcell fractional area (excluding dedicated wetland columns)
    !real(r8), pointer :: layer_sat_lag(:,:)      				! Lagged saturation status of soil layer in the unsaturated zone (1 = sat)

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

!
!EOP
!
! !OTHER LOCAL VARIABLES:
	integer :: j,l,c,p,k   					! indices
	integer :: begp, endp   				! per-proc beginning and ending pft indices
	integer :: begc, endc   				! per-proc beginning and ending column indices
	integer :: begl, endl   				! per-proc beginning and ending landunit indices
	integer :: begg, endg   				! per-proc gridcell ending gridcell indices
!-----------------------------------------------------------------------
    
	if ( masterproc ) write (iulog,*) 'Setting initial data to non-spun up values for CH4 process in microbial Mod,', &
                                  'if no inicFile or no valid values for concentrations,', &
                                  'for microbe Model.'

    ! Assign local pointers to derived subtypes components (landunit-level)

	ltype      					=> lun%itype

    ! Assign local pointers to derived subtypes components (column-level)

	clandunit          				=> col%landunit
    
	cmicbiocs          				=> cmic%cmicbiocs
	cdocs          				=> cmic%cdocs
	cdocs_pre          			=> cmic%cdocs_pre
	cmicbions          				=> cmic%cmicbions
	cdons          				=> cmic%cdons	
	caces          				=> cmic%caces
	cacebios          				=> cmic%cacebios
	cco2bios          				=> cmic%cco2bios
	caerch4bios          			=> cmic%caerch4bios
	canaerch4bios          			=> cmic%canaerch4bios
    
	cdocs_unsat          			=> cmic%cdocs_unsat
	cdocs_sat          			=> cmic%cdocs_sat
	cdons_unsat          			=> cmic%cdons_unsat
	cdons_sat          			=> cmic%cdons_sat
	cdons_min          			=> cmic%cdons_min
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
	
        cn_microbe      				=> cmic%cn_microbe
	root2doc           				=> cmic%root2doc
	froot_r           				=> cmic%froot_r
	micbio_hr_vr           			=> cmic%micbio_hr_vr
	micbio_hr           				=> cmic%micbio_hr
	micbiohr_col           			=> ccs%micbiohr_col
	micbioc_col           			=> ccs%micbioc_col
	doc_col           				=> ccs%doc_col
	dochr_vr           				=> ccs%dochr_vr
	dochr_col           			=> ccs%dochr_col
	ace_col           				=> ccs%ace_col
	ace_prod_col           			=> ccs%ace_prod_col
	acebios_col           			=> ccs%acebios_col
	co2bio_col           			=> ccs%co2bio_col
	aerch4bio_col           			=> ccs%aerch4bio_col
	anaerch4bio_col           		=> ccs%anaerch4bio_col
	ccon_ch4s_col				=> ccs%ccon_ch4s_col
	ccon_co2s_col				=> ccs%ccon_co2s_col

	micbion_col           			=> cns%micbion_col
	don_col           				=> cns%don_col
	
	caces_unsat_prod			=> cmic%caces_unsat_prod
	caces_sat_prod				=> cmic%caces_sat_prod
	caces_prod					=> cmic%caces_prod

	caces_unsat_prod_h2			=> cmic%caces_unsat_prod_h2
	caces_sat_prod_h2			=> cmic%caces_sat_prod_h2
	caces_prod_h2				=> cmic%caces_prod_h2
	
	ccon_ch4s          			=> cmic%ccon_ch4s
	ccon_co2s          			=> cmic%ccon_co2s
	ccon_o2s         				=> cmic%ccon_o2s
	ccon_h2s          				=> cmic%ccon_h2s
	ccon_ch4s_unsat    			=> cmic%ccon_ch4s_unsat
	ccon_ch4s_sat        			=> cmic%ccon_ch4s_sat
	ccon_co2s_unsat    			=> cmic%ccon_co2s_unsat
	ccon_co2s_sat        			=> cmic%ccon_co2s_sat
	ccon_o2s_unsat      			=> cmic%ccon_o2s_unsat
	ccon_o2s_sat          			=> cmic%ccon_o2s_sat
	ccon_h2s_unsat      			=> cmic%ccon_h2s_unsat
	ccon_h2s_sat          			=> cmic%ccon_h2s_sat
   
    !qflx_surf_lag      => cmic%qflx_surf_lag
    !finundated_lag     => cmic%finundated_lag
    !o2stress_sat       => cmic%o2stress_sat
    !o2stress_unsat     => cmic%o2stress_unsat
	finundated         				=> cws%finundated
	fsat_pre           				=> cmic%fsat_pre
	waterhead_unsat   			=> cmic%waterhead_unsat
    !layer_sat_lag      => cmic%layer_sat_lag

    
    	ch4_prod_ace_depth_unsat	=> cmic%ch4_prod_ace_depth_unsat
	ch4_prod_co2_depth_unsat	=> cmic%ch4_prod_co2_depth_unsat
	ch4_oxid_o2_depth_unsat 		=> cmic%ch4_oxid_o2_depth_unsat
	ch4_oxid_aom_depth_unsat 	=> cmic%ch4_oxid_aom_depth_unsat
	ch4_aere_depth_unsat 		=> cmic%ch4_aere_depth_unsat
	ch4_dif_depth_unsat 		=> cmic%ch4_dif_depth_unsat
	ch4_ebul_depth_unsat 		=> cmic%ch4_ebul_depth_unsat
	co2_prod_ace_depth_unsat	=> cmic%co2_prod_ace_depth_unsat
	co2_decomp_depth_unsat 		=> cmic%co2_decomp_depth_unsat
	co2_cons_depth_unsat 		=> cmic%co2_cons_depth_unsat
	co2_ebul_depth_unsat 		=> cmic%co2_ebul_depth_unsat
	co2_aere_depth_unsat 		=> cmic%co2_aere_depth_unsat
	co2_dif_depth_unsat 			=> cmic%co2_dif_depth_unsat
	o2_cons_depth_unsat 		=> cmic%o2_cons_depth_unsat
	o2_aere_depth_unsat 		=> cmic%o2_aere_depth_unsat
	o2_aere_oxid_depth_unsat 	=> cmic%o2_aere_oxid_depth_unsat
	o2_decomp_depth_unsat 		=> cmic%o2_decomp_depth_unsat
	o2_cons_depth_unsat 		=> cmic%o2_cons_depth_unsat
	o2_dif_depth_unsat 			=> cmic%o2_dif_depth_unsat
	h2_prod_depth_unsat 		=> cmic%h2_prod_depth_unsat
	h2_cons_depth_unsat 		=> cmic%h2_cons_depth_unsat
	h2_aere_depth_unsat 		=> cmic%h2_aere_depth_unsat
	h2_diff_depth_unsat 			=> cmic%h2_diff_depth_unsat
	h2_ebul_depth_unsat 		=> cmic%h2_ebul_depth_unsat
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
	o2_surf_netflux_unsat 		=> cmic%o2_surf_netflux_unsat
	h2_surf_aere_unsat			=> cmic%h2_surf_aere_unsat
	h2_surf_ebul_unsat 			=> cmic%h2_surf_ebul_unsat
	h2_surf_dif_unsat 			=> cmic%h2_surf_dif_unsat
	h2_surf_netflux_unsat 		=> cmic%h2_surf_netflux_unsat

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
	o2_decomp_depth_sat 		=> cmic%o2_decomp_depth_sat
	o2_cons_depth_sat 			=> cmic%o2_cons_depth_sat
	o2_dif_depth_sat 			=> cmic%o2_dif_depth_sat
	h2_prod_depth_sat 			=> cmic%h2_prod_depth_sat
	h2_cons_depth_sat 			=> cmic%h2_cons_depth_sat
	h2_aere_depth_sat 			=> cmic%h2_aere_depth_sat
	h2_diff_depth_sat 			=> cmic%h2_diff_depth_sat
	h2_ebul_depth_sat 			=> cmic%h2_ebul_depth_sat
	ch4_surf_aere_sat			=> cmic%ch4_surf_aere_sat
	ch4_surf_ebul_sat 			=> cmic%ch4_surf_ebul_sat
	ch4_surf_dif_sat 			=> cmic%ch4_surf_dif_sat
	ch4_surf_netflux_sat 			=> cmic%ch4_surf_netflux_sat
	co2_surf_aere_sat			=> cmic%co2_surf_aere_sat
	co2_surf_ebul_sat 			=> cmic%co2_surf_ebul_sat
	co2_surf_dif_sat 			=> cmic%co2_surf_dif_sat
	co2_surf_netflux_sat 			=> cmic%co2_surf_netflux_sat
	o2_surf_aere_sat			=> cmic%o2_surf_aere_sat
	o2_surf_dif_sat 				=> cmic%o2_surf_dif_sat
	o2_surf_netflux_sat 			=> cmic%o2_surf_netflux_sat
	h2_surf_aere_sat			=> cmic%h2_surf_aere_sat
	h2_surf_ebul_sat 			=> cmic%h2_surf_ebul_sat
	h2_surf_dif_sat 			=> cmic%h2_surf_dif_sat
	h2_surf_netflux_sat 			=> cmic%h2_surf_netflux_sat

	ch4_prod_ace_depth 			=> cmic%ch4_prod_ace_depth
	ch4_prod_co2_depth 			=> cmic%ch4_prod_co2_depth
	ch4_oxid_o2_depth 			=> cmic%ch4_oxid_o2_depth
	ch4_oxid_aom_depth 		=> cmic%ch4_oxid_aom_depth
	ch4_aere_depth 			=> cmic%ch4_aere_depth
	ch4_dif_depth 				=> cmic%ch4_dif_depth
	ch4_ebul_depth 			=> cmic%ch4_ebul_depth
	co2_prod_ace_depth 			=> cmic%co2_prod_ace_depth
	co2_decomp_depth 			=> cmic%co2_decomp_depth
	co2_cons_depth 			=> cmic%co2_cons_depth
	co2_ebul_depth 			=> cmic%co2_ebul_depth
	co2_aere_depth 			=> cmic%co2_aere_depth
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
	ch4_surf_netflux 			=> cmic%ch4_surf_netflux
	co2_surf_aere				=> cmic%co2_surf_aere
	co2_surf_ebul 				=> cmic%co2_surf_ebul
	co2_surf_dif 				=> cmic%co2_surf_dif
	co2_surf_netflux 			=> cmic%co2_surf_netflux
	o2_surf_aere				=> cmic%o2_surf_aere
	o2_surf_dif 				=> cmic%o2_surf_dif
	o2_surf_netflux 			=> cmic%o2_surf_netflux
	h2_surf_aere				=> cmic%h2_surf_aere
	h2_surf_ebul 				=> cmic%h2_surf_ebul
	h2_surf_dif 				=> cmic%h2_surf_dif
	h2_surf_netflux 			=> cmic%h2_surf_netflux

    ! Assign local pointers to derived subtypes components (pft-level)
    ! Determine subgrid bounds on this processor

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

	do c = begc,endc
	l = clandunit(c)
	if(micbiohr_col(c) == spval .or. arbinit) 					micbiohr_col(c) = 1e-15_r8
	if(micbioc_col(c) == spval .or. arbinit) 					micbioc_col(c) = 1e-15_r8
	if(doc_col(c) == spval .or. arbinit) 					doc_col(c) = 1e-15_r8
	if(dochr_col(c) == spval .or. arbinit) 					dochr_col(c) = 1e-15_r8
	if(ace_col(c) == spval .or. arbinit) 						ace_col(c) = 1e-15_r8
	if(ace_prod_col(c) == spval .or. arbinit) 					ace_prod_col(c) = 1e-15_r8
	if(acebios_col(c) == spval .or. arbinit) 					acebios_col(c) = 1e-15_r8
	if(co2bio_col(c) == spval .or. arbinit) 					co2bio_col(c) = 1e-15_r8
	if(aerch4bio_col(c) == spval .or. arbinit) 					aerch4bio_col(c) = 1e-15_r8
	if(anaerch4bio_col(c) == spval .or. arbinit) 				anaerch4bio_col(c) = 1e-15_r8

	if(ch4_surf_aere_unsat(c) == spval .or. arbinit) 			ch4_surf_aere_unsat(c) = 0._r8
	if(ch4_surf_ebul_unsat(c) == spval .or. arbinit) 			ch4_surf_ebul_unsat(c) = 0._r8
	if(ch4_surf_dif_unsat(c) == spval .or. arbinit) 				ch4_surf_dif_unsat(c) = 0._r8
	if(ch4_surf_netflux_unsat(c) == spval .or. arbinit) 			ch4_surf_netflux_unsat(c) = 0._r8
	if(co2_surf_aere_unsat(c) == spval .or. arbinit) 			co2_surf_aere_unsat(c) = 0._r8
	if(co2_surf_ebul_unsat(c) == spval .or. arbinit) 			co2_surf_ebul_unsat(c) = 0._r8
	if(co2_surf_dif_unsat(c) == spval .or. arbinit) 				co2_surf_dif_unsat(c) = 0._r8
	if(co2_surf_netflux_unsat(c) == spval .or. arbinit) 			co2_surf_netflux_unsat(c) = 0._r8
	if(o2_surf_aere_unsat(c) == spval .or. arbinit) 				o2_surf_aere_unsat(c) = 0._r8
	if(o2_surf_dif_unsat(c) == spval .or. arbinit) 				o2_surf_dif_unsat(c) = 0._r8
	if(o2_surf_netflux_unsat(c) == spval .or. arbinit) 			o2_surf_netflux_unsat(c) = 0._r8
	if(h2_surf_aere_unsat(c) == spval .or. arbinit) 				h2_surf_aere_unsat(c) = 0._r8
	if(h2_surf_ebul_unsat(c) == spval .or. arbinit) 				h2_surf_ebul_unsat(c) = 0._r8
	if(h2_surf_dif_unsat(c) == spval .or. arbinit) 				h2_surf_dif_unsat(c) = 0._r8
	if(h2_surf_netflux_unsat(c) == spval .or. arbinit) 			h2_surf_netflux_unsat(c) = 0._r8

	if(ch4_surf_aere_sat(c) == spval .or. arbinit) 				ch4_surf_aere_sat(c) = 0._r8
	if(ch4_surf_ebul_sat(c) == spval .or. arbinit) 				ch4_surf_ebul_sat(c) = 0._r8
	if(ch4_surf_dif_sat(c) == spval .or. arbinit) 				ch4_surf_dif_sat(c) = 0._r8
	if(ch4_surf_netflux_sat(c) == spval .or. arbinit) 			ch4_surf_netflux_sat(c) = 0._r8
	if(co2_surf_aere_sat(c) == spval .or. arbinit) 				co2_surf_aere_sat(c) = 0._r8
	if(co2_surf_ebul_sat(c) == spval .or. arbinit) 				co2_surf_ebul_sat(c) = 0._r8
	if(co2_surf_dif_sat(c) == spval .or. arbinit) 				co2_surf_dif_sat(c) = 0._r8
	if(co2_surf_netflux_sat(c) == spval .or. arbinit) 			co2_surf_netflux_sat(c) = 0._r8
	if(o2_surf_aere_sat(c) == spval .or. arbinit) 				o2_surf_aere_sat(c) = 0._r8
	if(o2_surf_dif_sat(c) == spval .or. arbinit) 				o2_surf_dif_sat(c) = 0._r8
	if(o2_surf_netflux_sat(c) == spval .or. arbinit) 				o2_surf_netflux_sat(c) = 0._r8
	if(h2_surf_aere_sat(c) == spval .or. arbinit) 				h2_surf_aere_sat(c) = 0._r8
	if(h2_surf_ebul_sat(c) == spval .or. arbinit) 				h2_surf_ebul_sat(c) = 0._r8
	if(h2_surf_dif_sat(c) == spval .or. arbinit) 				h2_surf_dif_sat(c) = 0._r8
	if(h2_surf_netflux_sat(c) == spval .or. arbinit) 			h2_surf_netflux_sat(c) = 0._r8

	if(ch4_surf_aere(c) == spval .or. arbinit) 					ch4_surf_aere(c) = 0._r8
	if(ch4_surf_ebul(c) == spval .or. arbinit) 					ch4_surf_ebul(c) = 0._r8
	if(ch4_surf_dif(c) == spval .or. arbinit) 					ch4_surf_dif(c) = 0._r8
	if(ch4_surf_netflux(c) == spval .or. arbinit) 				ch4_surf_netflux(c) = 0._r8
	if(co2_surf_aere(c) == spval .or. arbinit) 					co2_surf_aere(c) = 0._r8
	if(co2_surf_ebul(c) == spval .or. arbinit) 					co2_surf_ebul(c) = 0._r8
	if(co2_surf_dif(c) == spval .or. arbinit) 					co2_surf_dif(c) = 0._r8
	if(co2_surf_netflux(c) == spval .or. arbinit) 				co2_surf_netflux(c) = 0._r8
	if(o2_surf_aere(c) == spval .or. arbinit) 					o2_surf_aere(c) = 0._r8
	if(o2_surf_dif(c) == spval .or. arbinit) 					o2_surf_dif(c) = 0._r8
	if(o2_surf_netflux(c) == spval .or. arbinit) 				o2_surf_netflux(c) = 0._r8
	if(h2_surf_aere(c) == spval .or. arbinit) 					h2_surf_aere(c) = 0._r8
	if(h2_surf_ebul(c) == spval .or. arbinit) 					h2_surf_ebul(c) = 0._r8
	if(h2_surf_dif(c) == spval .or. arbinit) 					h2_surf_dif(c) = 0._r8
	if(h2_surf_netflux(c) == spval .or. arbinit) 				h2_surf_netflux(c) = 0._r8
	
	if (ltype(l) == istsoil .or. ltype(l) == istcrop) then      
        do j=1,nlevgrnd
	if(dochr_vr(c,j) == spval .or. arbinit) 					dochr_vr(c,j) = 0._r8
	if(root2doc(c,j) == spval .or. arbinit) 					root2doc(c,j) = 0._r8
	if(cn_microbe(c,j) == spval .or. arbinit) 					cn_microbe(c,j) = 10_r8
	if(froot_r(c,j) == spval .or. arbinit) 						froot_r(c,j) = 0._r8
	if(cmicbiocs(c,j) == spval .or. arbinit) 					cmicbiocs(c,j) = 0._r8
	if(cdocs_pre(c,j) == spval .or. arbinit) 					cdocs_pre(c,j) = 0._r8
	if(cdocs(c,j) == spval .or. arbinit) 						cdocs(c,j) = 0._r8
	if(cdocs_unsat(c,j) == spval .or. arbinit) 					cdocs_unsat(c,j) = 0._r8
	if(cdocs_sat(c,j) == spval .or. arbinit) 					cdocs_sat(c,j) = 0._r8
	if(cmicbions(c,j) == spval .or. arbinit) 					cmicbions(c,j) = 0._r8
	if(cdons(c,j) == spval .or. arbinit) 						cdons(c,j) = 0._r8
	if(cdons_unsat(c,j) == spval .or. arbinit) 					cdons_unsat(c,j) = 0._r8
	if(cdons_sat(c,j) == spval .or. arbinit) 					cdons_sat(c,j) = 0._r8
	if(cdons_min(c,j) == spval .or. arbinit) 					cdons_min(c,j) = 0._r8
	if(caces(c,j) == spval .or. arbinit) 						caces(c,j) = 0._r8
	if(caces_unsat(c,j) == spval .or. arbinit) 					caces_unsat(c,j) = 0._r8
	if(caces_sat(c,j) == spval .or. arbinit) 					caces_sat(c,j) = 0._r8
	if(caces_prod(c,j) == spval .or. arbinit) 					caces_prod(c,j) = 0._r8
	if(caces_unsat_prod(c,j) == spval .or. arbinit) 			caces_unsat_prod(c,j) = 0._r8
	if(caces_sat_prod(c,j) == spval .or. arbinit) 				caces_sat_prod(c,j) = 0._r8
	if(caces_prod_h2(c,j) == spval .or. arbinit) 				caces_prod_h2(c,j) = 0._r8
	if(caces_unsat_prod_h2(c,j) == spval .or. arbinit) 			caces_unsat_prod_h2(c,j) = 0._r8
	if(caces_sat_prod_h2(c,j) == spval .or. arbinit) 			caces_sat_prod_h2(c,j) = 0._r8
	if(cacebios(c,j) == spval .or. arbinit) 					cacebios(c,j) = 1e-15_r8
	if(cacebios_unsat(c,j) == spval .or. arbinit) 				cacebios_unsat(c,j) = 1e-15_r8
	if(cacebios_sat(c,j) == spval .or. arbinit) 				cacebios_sat(c,j) = 1e-15_r8
	if(cco2bios(c,j) == spval .or. arbinit) 					cco2bios(c,j) = 1e-15_r8
	if(cco2bios_unsat(c,j) == spval .or. arbinit) 				cco2bios_unsat(c,j) = 1e-15_r8
	if(cco2bios_sat(c,j) == spval .or. arbinit) 				cco2bios_sat(c,j) = 1e-15_r8
	if(caerch4bios(c,j) == spval .or. arbinit) 					caerch4bios(c,j) = 1e-15_r8
	if(caerch4bios_unsat(c,j) == spval .or. arbinit) 				caerch4bios_unsat(c,j) = 1e-15_r8
	if(caerch4bios_sat(c,j) == spval .or. arbinit) 				caerch4bios(c,j) = 1e-15_r8  
	if(canaerch4bios(c,j) == spval .or. arbinit) 				canaerch4bios(c,j) = 1e-15_r8
	if(canaerch4bios_unsat(c,j) == spval .or. arbinit) 			canaerch4bios_unsat(c,j) = 1e-15_r8
	if(canaerch4bios_sat(c,j) == spval .or. arbinit) 				canaerch4bios_sat(c,j) = 1e-15_r8

	if(ccon_ch4s(c,j) == spval .or. arbinit) 					ccon_ch4s(c,j) = 0._r8
	if(ccon_co2s(c,j) == spval .or. arbinit) 					ccon_co2s(c,j) = 0._r8
	if(ccon_h2s(c,j) == spval .or. arbinit) 					ccon_h2s(c,j) = 0._r8
	if(ccon_o2s(c,j) == spval .or. arbinit) 					ccon_o2s(c,j) =0._r8
	if(ccon_ch4s_unsat(c,j) == spval .or. arbinit) 				ccon_ch4s_unsat(c,j) = 0._r8
	if(ccon_co2s_unsat(c,j) == spval .or. arbinit) 				ccon_co2s_unsat(c,j) = 0._r8
	if(ccon_h2s_unsat(c,j) == spval .or. arbinit) 				ccon_h2s_unsat(c,j) = 0._r8
	if(ccon_o2s_unsat(c,j) == spval .or. arbinit) 				ccon_o2s_unsat(c,j) = 0._r8
	if(ccon_ch4s_sat(c,j) == spval .or. arbinit) 				ccon_ch4s_sat(c,j) = 0._r8
	if(ccon_co2s_sat(c,j) == spval .or. arbinit) 				ccon_co2s_sat(c,j) = 0._r8
	if(ccon_h2s_sat(c,j) == spval .or. arbinit) 				ccon_h2s_sat(c,j) = 0._r8
	if(ccon_o2s_sat(c,j) == spval .or. arbinit) 				ccon_o2s_sat(c,j) = 0._r8
	     
	if(ch4_prod_ace_depth_unsat(c,j) == spval .or. arbinit)		ch4_prod_ace_depth_unsat(c,j) = 0._r8
	if(ch4_prod_co2_depth_unsat(c,j) == spval .or. arbinit)		ch4_prod_co2_depth_unsat(c,j) = 0._r8
	if(ch4_oxid_o2_depth_unsat(c,j) == spval .or. arbinit)		ch4_oxid_o2_depth_unsat(c,j) = 0._r8
	if(ch4_oxid_aom_depth_unsat(c,j) == spval .or. arbinit)		ch4_oxid_aom_depth_unsat(c,j) = 0._r8
	if(ch4_aere_depth_unsat(c,j) == spval .or. arbinit)			ch4_aere_depth_unsat(c,j) = 0._r8
	if(ch4_dif_depth_unsat(c,j) == spval .or. arbinit)			ch4_dif_depth_unsat(c,j) = 0._r8
	if(ch4_ebul_depth_unsat(c,j) == spval .or. arbinit)			ch4_ebul_depth_unsat(c,j) = 0._r8
	if(co2_prod_ace_depth_unsat(c,j) == spval .or. arbinit)		co2_prod_ace_depth_unsat(c,j) = 0._r8
	if(co2_decomp_depth_unsat(c,j) == spval .or. arbinit)		co2_decomp_depth_unsat(c,j) = 0._r8
	if(co2_cons_depth_unsat(c,j) == spval .or. arbinit)			co2_cons_depth_unsat(c,j) = 0._r8
	if(co2_ebul_depth_unsat(c,j) == spval .or. arbinit)			co2_ebul_depth_unsat(c,j) = 0._r8
	if(co2_aere_depth_unsat(c,j) == spval .or. arbinit)			co2_aere_depth_unsat(c,j) = 0._r8
	if(co2_dif_depth_unsat(c,j) == spval .or. arbinit)			co2_dif_depth_unsat(c,j) = 0._r8
	if(o2_cons_depth_unsat(c,j) == spval .or. arbinit)			o2_cons_depth_unsat(c,j) = 0._r8
	if(o2_aere_depth_unsat(c,j) == spval .or. arbinit)			o2_aere_depth_unsat(c,j) = 0._r8
	if(o2_aere_oxid_depth_unsat(c,j) == spval .or. arbinit)		o2_aere_oxid_depth_unsat(c,j) = 0._r8
	if(o2_decomp_depth_unsat(c,j) == spval .or. arbinit)		o2_decomp_depth_unsat(c,j) = 0._r8
	if(o2_cons_depth_unsat(c,j) == spval .or. arbinit)			o2_cons_depth_unsat(c,j) = 0._r8
	if(o2_dif_depth_unsat(c,j) == spval .or. arbinit)			o2_dif_depth_unsat(c,j) = 0._r8
	if(h2_prod_depth_unsat(c,j) == spval .or. arbinit)			h2_prod_depth_unsat(c,j) = 0._r8
	if(h2_cons_depth_unsat(c,j) == spval .or. arbinit)			h2_cons_depth_unsat(c,j) = 0._r8
	if(h2_cons_depth_unsat(c,j) == spval .or. arbinit)			h2_cons_depth_unsat(c,j) = 0._r8
	if(h2_cons_depth_unsat(c,j) == spval .or. arbinit)			h2_cons_depth_unsat(c,j) = 0._r8
	if(h2_aere_depth_unsat(c,j) == spval .or. arbinit)			h2_aere_depth_unsat(c,j) = 0._r8
	if(h2_diff_depth_unsat(c,j) == spval .or. arbinit)			h2_diff_depth_unsat(c,j) = 0._r8
	if(h2_ebul_depth_unsat(c,j) == spval .or. arbinit)			h2_ebul_depth_unsat(c,j) = 0._r8

	if(ch4_prod_ace_depth_sat(c,j) == spval .or. arbinit)		ch4_prod_ace_depth_sat(c,j) = 0._r8
	if(ch4_prod_co2_depth_sat(c,j) == spval .or. arbinit)		ch4_prod_co2_depth_sat(c,j) = 0._r8
	if(ch4_oxid_o2_depth_sat(c,j) == spval .or. arbinit)			ch4_oxid_o2_depth_sat(c,j) = 0._r8
	if(ch4_oxid_aom_depth_sat(c,j) == spval .or. arbinit)		ch4_oxid_aom_depth_sat(c,j) = 0._r8
	if(ch4_aere_depth_sat(c,j) == spval .or. arbinit)			ch4_aere_depth_sat(c,j) = 0._r8
	if(ch4_dif_depth_sat(c,j) == spval .or. arbinit)				ch4_dif_depth_sat(c,j) = 0._r8
	if(ch4_ebul_depth_sat(c,j) == spval .or. arbinit)			ch4_ebul_depth_sat(c,j) = 0._r8
	if(co2_prod_ace_depth_sat(c,j) == spval .or. arbinit)		co2_prod_ace_depth_sat(c,j) = 0._r8
	if(co2_decomp_depth_sat(c,j) == spval .or. arbinit)			co2_decomp_depth_sat(c,j) = 0._r8
	if(co2_cons_depth_sat(c,j) == spval .or. arbinit)			co2_cons_depth_sat(c,j) = 0._r8
	if(co2_ebul_depth_sat(c,j) == spval .or. arbinit)			co2_ebul_depth_sat(c,j) = 0._r8
	if(co2_aere_depth_sat(c,j) == spval .or. arbinit)			co2_aere_depth_sat(c,j) = 0._r8
	if(co2_dif_depth_sat(c,j) == spval .or. arbinit)				co2_dif_depth_sat(c,j) = 0._r8
	if(o2_cons_depth_sat(c,j) == spval .or. arbinit)			o2_cons_depth_sat(c,j) = 0._r8
	if(o2_aere_depth_sat(c,j) == spval .or. arbinit)				o2_aere_depth_sat(c,j) = 0._r8
	if(o2_aere_oxid_depth_sat(c,j) == spval .or. arbinit)			o2_aere_oxid_depth_sat(c,j) = 0._r8
	if(o2_decomp_depth_sat(c,j) == spval .or. arbinit)			o2_decomp_depth_sat(c,j) = 0._r8
	if(o2_cons_depth_sat(c,j) == spval .or. arbinit)			o2_cons_depth_sat(c,j) = 0._r8
	if(o2_dif_depth_sat(c,j) == spval .or. arbinit)				o2_dif_depth_sat(c,j) = 0._r8
	if(h2_prod_depth_sat(c,j) == spval .or. arbinit)			h2_prod_depth_sat(c,j) = 0._r8
	if(h2_cons_depth_sat(c,j) == spval .or. arbinit)			h2_cons_depth_sat(c,j) = 0._r8
	if(h2_cons_depth_sat(c,j) == spval .or. arbinit)			h2_cons_depth_sat(c,j) = 0._r8
	if(h2_cons_depth_sat(c,j) == spval .or. arbinit)			h2_cons_depth_sat(c,j) = 0._r8
	if(h2_aere_depth_sat(c,j) == spval .or. arbinit)			h2_aere_depth_sat(c,j) = 0._r8
	if(h2_diff_depth_sat(c,j) == spval .or. arbinit)				h2_diff_depth_sat(c,j) = 0._r8
	if(h2_ebul_depth_sat(c,j) == spval .or. arbinit)			h2_ebul_depth_sat(c,j) = 0._r8

	if(ch4_prod_ace_depth(c,j) == spval .or. arbinit)			ch4_prod_ace_depth(c,j) = 0._r8
	if(ch4_prod_co2_depth(c,j) == spval .or. arbinit)			ch4_prod_co2_depth(c,j) = 0._r8
	if(ch4_oxid_o2_depth(c,j) == spval .or. arbinit)			ch4_oxid_o2_depth(c,j) = 0._r8
	if(ch4_oxid_aom_depth(c,j) == spval .or. arbinit)			ch4_oxid_aom_depth(c,j) = 0._r8
	if(ch4_aere_depth(c,j) == spval .or. arbinit)				ch4_aere_depth(c,j) = 0._r8
	if(ch4_dif_depth(c,j) == spval .or. arbinit)				ch4_dif_depth(c,j) = 0._r8
	if(ch4_ebul_depth(c,j) == spval .or. arbinit)				ch4_ebul_depth(c,j) = 0._r8
	if(co2_prod_ace_depth(c,j) == spval .or. arbinit)			co2_prod_ace_depth(c,j) = 0._r8
	if(co2_decomp_depth(c,j) == spval .or. arbinit)			co2_decomp_depth(c,j) = 0._r8
	if(co2_cons_depth(c,j) == spval .or. arbinit)				co2_cons_depth(c,j) = 0._r8
	if(co2_ebul_depth(c,j) == spval .or. arbinit)				co2_ebul_depth(c,j) = 0._r8
	if(co2_aere_depth(c,j) == spval .or. arbinit)				co2_aere_depth(c,j) = 0._r8
	if(co2_dif_depth(c,j) == spval .or. arbinit)				co2_dif_depth(c,j) = 0._r8
	if(o2_cons_depth(c,j) == spval .or. arbinit)				o2_cons_depth(c,j) = 0._r8
	if(o2_aere_depth(c,j) == spval .or. arbinit)				o2_aere_depth(c,j) = 0._r8
	if(o2_aere_oxid_depth(c,j) == spval .or. arbinit)			o2_aere_oxid_depth(c,j) = 0._r8
	if(o2_decomp_depth(c,j) == spval .or. arbinit)				o2_decomp_depth(c,j) = 0._r8
	if(o2_cons_depth(c,j) == spval .or. arbinit)				o2_cons_depth(c,j) = 0._r8
	if(o2_dif_depth(c,j) == spval .or. arbinit)					o2_dif_depth(c,j) = 0._r8
	if(h2_prod_depth(c,j) == spval .or. arbinit)				h2_prod_depth(c,j) = 0._r8
	if(h2_cons_depth(c,j) == spval .or. arbinit)				h2_cons_depth(c,j) = 0._r8
	if(h2_cons_depth(c,j) == spval .or. arbinit)				h2_cons_depth(c,j) = 0._r8
	if(h2_cons_depth(c,j) == spval .or. arbinit)				h2_cons_depth(c,j) = 0._r8
	if(h2_aere_depth(c,j) == spval .or. arbinit)				h2_aere_depth(c,j) = 0._r8
	if(h2_diff_depth(c,j) == spval .or. arbinit)				h2_diff_depth(c,j) = 0._r8
	if(h2_ebul_depth(c,j) == spval .or. arbinit)				h2_ebul_depth(c,j) = 0._r8
	
	do k=1,ndecomp_cascade_transitions
	if(micbio_hr_vr(c,j,k) == spval .or. arbinit)				micbio_hr_vr(c,j,k) = 0._r8
	end do
	
        end do
          !if (qflx_surf_lag(c) == spval .or. arbinit) qflx_surf_lag(c) = 0._r8
          !if (finundated_lag(c) == spval .or. arbinit) finundated_lag(c) = 0._r8
          ! finundated will be used to calculate soil decomposition if anoxia is used
        if (fsat_pre(c) == spval .or. arbinit) then
		finundated(c) = 0._r8
        else
		finundated(c) = fsat_pre(c)
	end if
	else if (ltype(l) == istdlak) then
        do j=1,nlevsoi
!  if (conc_ch4_sat(c,j) == spval .or. arbinit)   conc_ch4_sat(c,j)   = 0._r8
!  if (conc_o2_sat(c,j) == spval .or. arbinit)    conc_o2_sat(c,j)    = 0._r8
         ! Need to convert from kg/m^3 organic matter to g C / m^3 (org matter is defined to be 58% C)
	end do
	end if 

	do k=1,ndecomp_cascade_transitions
	if(micbio_hr(c,k) == spval .or. arbinit)				micbio_hr(c,k) = 0._r8
	end do
	
! Set values for all columns equal to zero below nlevsoi
	do j=nlevsoi+1,nlevgrnd
        cmicbiocs(c,j) = 0._r8
	cmicbions(c,j) = 0._r8
	dochr_vr(c,j) = 0._r8
!	micbio_hr(c,j) = 0._r8
	root2doc(c,j) = 0._r8
	froot_r(c,j) = 0._r8
	cn_microbe(c,j) = 10._r8
	  
	cdocs_pre(c,j) = 0._r8
	cdocs(c,j) = 0._r8
	cdocs_unsat(c,j) = 0._r8
	cdocs_sat(c,j) = 0._r8
	cdons(c,j) = 0._r8
	cdons_unsat(c,j) = 0._r8
	cdons_sat(c,j) = 0._r8
	cdons_min(c,j) = 0._r8
	caces(c,j) = 0._r8
	cacebios(c,j) = 0._r8
	cco2bios(c,j) = 0._r8
	caerch4bios(c,j) = 0._r8
	canaerch4bios(c,j) = 0._r8
	caces_unsat(c,j) = 0._r8
	cacebios_unsat(c,j) = 0._r8
	cco2bios_unsat(c,j) = 0._r8
	caerch4bios_unsat(c,j) = 0._r8
	canaerch4bios_unsat(c,j) = 0._r8
	caces_sat(c,j) = 0._r8
	cacebios_sat(c,j) = 0._r8
	cco2bios_sat(c,j) = 0._r8
	caerch4bios_sat(c,j) = 0._r8
	canaerch4bios_sat(c,j) = 0._r8
	caces_prod(c,j) = 0._r8
	caces_unsat_prod(c,j) = 0._r8
	caces_sat_prod(c,j) = 0._r8
	caces_prod_h2(c,j) = 0._r8
	caces_unsat_prod_h2(c,j) = 0._r8
	caces_sat_prod_h2(c,j) = 0._r8	
	ccon_ch4s(c,j) = 0._r8
	ccon_co2s(c,j) = 0._r8
	ccon_o2s(c,j) = 0._r8
	ccon_h2s(c,j) = 0._r8
	ccon_ch4s_unsat(c,j) = 0._r8
	ccon_co2s_unsat(c,j) = 0._r8
	ccon_o2s_unsat(c,j) = 0._r8
	ccon_h2s_unsat(c,j) = 0._r8
	ccon_ch4s_sat(c,j) = 0._r8
	ccon_co2s_sat(c,j) = 0._r8
	ccon_o2s_sat(c,j) = 0._r8
	ccon_h2s_sat(c,j) = 0._r8
	  
	ch4_prod_ace_depth_unsat(c,j) = 0._r8
	ch4_prod_co2_depth_unsat(c,j) = 0._r8
	ch4_oxid_o2_depth_unsat(c,j) = 0._r8
	ch4_oxid_aom_depth_unsat(c,j) = 0._r8
	ch4_aere_depth_unsat(c,j) = 0._r8
	ch4_dif_depth_unsat(c,j) = 0._r8
	ch4_ebul_depth_unsat(c,j) = 0._r8
	co2_prod_ace_depth_unsat(c,j) = 0._r8
	co2_decomp_depth_unsat(c,j) = 0._r8
	co2_cons_depth_unsat(c,j) = 0._r8
	co2_ebul_depth_unsat(c,j) = 0._r8
	co2_aere_depth_unsat(c,j) = 0._r8
	co2_dif_depth_unsat(c,j) = 0._r8
	o2_cons_depth_unsat(c,j) = 0._r8
	o2_aere_depth_unsat(c,j) = 0._r8
	o2_aere_oxid_depth_unsat(c,j) = 0._r8
	o2_decomp_depth_unsat(c,j) = 0._r8
	o2_cons_depth_unsat(c,j) = 0._r8
	o2_dif_depth_unsat(c,j) = 0._r8
	h2_prod_depth_unsat(c,j) = 0._r8
	h2_cons_depth_unsat(c,j) = 0._r8
	h2_cons_depth_unsat(c,j) = 0._r8
	h2_cons_depth_unsat(c,j) = 0._r8
	h2_aere_depth_unsat(c,j) = 0._r8
	h2_diff_depth_unsat(c,j) = 0._r8
	h2_ebul_depth_unsat(c,j) = 0._r8

	ch4_prod_ace_depth_sat(c,j) = 0._r8
	ch4_prod_co2_depth_sat(c,j) = 0._r8
	ch4_oxid_o2_depth_sat(c,j) = 0._r8
	ch4_oxid_aom_depth_sat(c,j) = 0._r8
	ch4_aere_depth_sat(c,j) = 0._r8
	ch4_dif_depth_sat(c,j) = 0._r8
	ch4_ebul_depth_sat(c,j) = 0._r8
	co2_prod_ace_depth_sat(c,j) = 0._r8
	co2_decomp_depth_sat(c,j) = 0._r8
	co2_cons_depth_sat(c,j) = 0._r8
	co2_ebul_depth_sat(c,j) = 0._r8
	co2_aere_depth_sat(c,j) = 0._r8
	co2_dif_depth_sat(c,j) = 0._r8
	o2_cons_depth_sat(c,j) = 0._r8
	o2_aere_depth_sat(c,j) = 0._r8
	o2_aere_oxid_depth_sat(c,j) = 0._r8
	o2_decomp_depth_sat(c,j) = 0._r8
	o2_cons_depth_sat(c,j) = 0._r8
	o2_dif_depth_sat(c,j) = 0._r8
	h2_prod_depth_sat(c,j) = 0._r8
	h2_cons_depth_sat(c,j) = 0._r8
	h2_cons_depth_sat(c,j) = 0._r8
	h2_cons_depth_sat(c,j) = 0._r8
	h2_aere_depth_sat(c,j) = 0._r8
	h2_diff_depth_sat(c,j) = 0._r8
	h2_ebul_depth_sat(c,j) = 0._r8
	  
	ch4_prod_ace_depth(c,j) = 0._r8
	ch4_prod_co2_depth(c,j) = 0._r8
	ch4_oxid_o2_depth(c,j) = 0._r8
	ch4_oxid_aom_depth(c,j) = 0._r8
	ch4_aere_depth(c,j) = 0._r8
	ch4_dif_depth(c,j) = 0._r8
	ch4_ebul_depth(c,j) = 0._r8
	co2_prod_ace_depth(c,j) = 0._r8
	co2_decomp_depth(c,j) = 0._r8
	co2_cons_depth(c,j) = 0._r8
	co2_ebul_depth(c,j) = 0._r8
	co2_aere_depth(c,j) = 0._r8
	co2_dif_depth(c,j) = 0._r8
	o2_cons_depth(c,j) = 0._r8
	o2_aere_depth(c,j) = 0._r8
	o2_aere_oxid_depth(c,j) = 0._r8
	o2_decomp_depth(c,j) = 0._r8
	o2_cons_depth(c,j) = 0._r8
	o2_dif_depth(c,j) = 0._r8
	h2_prod_depth(c,j) = 0._r8
	h2_cons_depth(c,j) = 0._r8
	h2_cons_depth(c,j) = 0._r8
	h2_cons_depth(c,j) = 0._r8
	h2_aere_depth(c,j) = 0._r8
	h2_diff_depth(c,j) = 0._r8
	h2_ebul_depth(c,j) = 0._r8

	end do
	end do

  end subroutine makearbinit_microbe


!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: initTimeConst_microbe
!
! !INTERFACE:
subroutine initTimeConst_microbe
!
! !DESCRIPTION:
! Initialize variables for ch4 code that will not be input from restart/inic file. Also set values for
! inactive CH4 columns to spval so that they will not be averaged in history file.
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use clmtype
  use decompMod   , only : get_proc_bounds, get_proc_global
  use clm_varpar  , only : nlevsoi, ngases, nlevgrnd, nlevdecomp
  use clm_varcon  , only : istsoil, istdlak, spval, istcrop
  use clm_varctl  , only : iulog
  use spmdMod ,     only : masterproc
  !use microbevarcon   , only : allowlakeprod
!
! !ARGUMENTS:
  implicit none
!
! !CALLED FROM:
! subroutine initialize2 in module initializeMod.
!
! !REVISION HISTORY:
! 9/09, Zack Subin
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
  integer , pointer :: clandunit(:)       ! landunit index of column
  integer , pointer :: cgridcell(:)       ! gridcell index of column
  integer , pointer :: ltype(:)           ! landunit type index
  real(r8), pointer :: watsat(:,:)        ! volumetric soil water at saturation (porosity)
  real(r8), pointer :: cellorg(:,:)       ! column 3D organic matter (kg/m^3, 58% by mass carbon) (nlevsoi)

!
! local pointers to implicit out arguments
!
  real(r8), pointer :: cmicbiocs(:,:)
  real(r8), pointer :: cdocs_pre(:,:)
  real(r8), pointer :: cdocs(:,:)
  real(r8), pointer :: cdocs_unsat(:,:)
  real(r8), pointer :: cdocs_sat(:,:)
  real(r8), pointer :: cmicbions(:,:)
  real(r8), pointer :: cdons(:,:)
  real(r8), pointer :: cdons_unsat(:,:)
  real(r8), pointer :: cdons_sat(:,:)  
  real(r8), pointer :: cdons_min(:,:)
  real(r8), pointer :: caces(:,:)
  real(r8), pointer :: cacebios(:,:)
  real(r8), pointer :: cco2bios(:,:)
  real(r8), pointer :: caerch4bios(:,:)
  real(r8), pointer :: canaerch4bios(:,:)
  real(r8), pointer :: caces_unsat(:,:)
  real(r8), pointer :: cacebios_unsat(:,:)
  real(r8), pointer :: cco2bios_unsat(:,:)
  real(r8), pointer :: caerch4bios_unsat(:,:)
  real(r8), pointer :: canaerch4bios_unsat(:,:)
  real(r8), pointer :: caces_sat(:,:)
  real(r8), pointer :: cacebios_sat(:,:)
  real(r8), pointer :: cco2bios_sat(:,:)
  real(r8), pointer :: caerch4bios_sat(:,:)
  real(r8), pointer :: canaerch4bios_sat(:,:)
  real(r8), pointer :: caces_prod(:,:)
  real(r8), pointer :: caces_sat_prod(:,:)
  real(r8), pointer :: caces_unsat_prod(:,:)
  real(r8), pointer :: caces_prod_h2(:,:)
  real(r8), pointer :: caces_sat_prod_h2(:,:)
  real(r8), pointer :: caces_unsat_prod_h2(:,:)
  real(r8), pointer :: dochr_vr(:,:)
  real(r8), pointer :: micbio_hr_vr(:,:,:)
  real(r8), pointer :: micbio_hr(:,:)
  real(r8), pointer :: root2doc(:,:)
  real(r8), pointer :: froot_r(:,:)
  real(r8), pointer :: cn_microbe(:,:)
  
  real(r8), pointer :: ccon_ch4s(:,:)
  real(r8), pointer :: ccon_co2s(:,:)
  real(r8), pointer :: ccon_o2s(:,:)
  real(r8), pointer :: ccon_h2s(:,:)
  real(r8), pointer :: ccon_ch4s_unsat(:,:)
  real(r8), pointer :: ccon_co2s_unsat(:,:)
  real(r8), pointer :: ccon_o2s_unsat(:,:)
  real(r8), pointer :: ccon_h2s_unsat(:,:)
  real(r8), pointer :: ccon_ch4s_sat(:,:)
  real(r8), pointer :: ccon_co2s_sat(:,:)
  real(r8), pointer :: ccon_o2s_sat(:,:)
  real(r8), pointer :: ccon_h2s_sat(:,:)
 
  real(r8), pointer :: finundated(:)               ! inundated gridcell fractional area
  !real(r8), pointer :: fphr(:,:)                   ! fraction of potential HR
  !real(r8), pointer :: sif(:)                      ! (unitless) ratio applied to sat. prod. to account for seasonal inundation
  !real(r8), pointer :: rootfr(:,:)                 ! column-averaged root fraction
  !real(r8), pointer :: o2stress_unsat(:,:)         ! Ratio of oxygen available to that demanded by roots, aerobes, & methanotrophs (nlevsoi)
  !real(r8), pointer :: o2stress_sat(:,:)           ! Ratio of oxygen available to that demanded by roots, aerobes, & methanotrophs (nlevsoi)
  !real(r8), pointer :: ch4stress_unsat(:,:)        ! Ratio of methane available to the total per-timestep methane sinks (nlevsoi)
  !real(r8), pointer :: ch4stress_sat(:,:)          ! Ratio of methane available to the total per-timestep methane sinks (nlevsoi)
  !real(r8), pointer :: totcolch4(:)                ! total methane in soil column (g C / m^2)

! To avoid rare pathologies with allowlakeprod switching between restarts
  !real(r8), pointer :: conc_ch4_sat(:,:)       ! CH4 conc in each soil layer (mol/m3) (nlevsoi)
  !real(r8), pointer :: conc_o2_sat(:,:)        ! O2 conc in each soil layer (mol/m3) (nlevsoi)
  
!
!
!EOP
!
! !OTHER LOCAL VARIABLES:
  integer  :: g,l,c,p,j        ! indices
  integer  :: begp, endp       ! per-proc beginning and ending pft indices
  integer  :: begc, endc       ! per-proc beginning and ending column indices
  integer  :: begl, endl       ! per-proc beginning and ending landunit indices
  integer  :: begg, endg       ! per-proc gridcell ending gridcell indices
	
!  watsat     => cps%watsat
	cmicbiocs          				=> cmic%cmicbiocs
	
	cdocs          				=> cmic%cdocs
	cdocs_pre          				=> cmic%cdocs_pre
	cdocs_unsat          			=> cmic%cdocs_unsat
	cdocs_sat          			=> cmic%cdocs_sat
	cmicbions          			=> cmic%cmicbions
	cdons          				=> cmic%cdons
	cdons_unsat          			=> cmic%cdons_unsat
	cdons_sat          			=> cmic%cdons_sat
	cdons_min          			=> cmic%cdons_min
	
	caces          				=> cmic%caces
	cacebios          				=> cmic%cacebios
	cco2bios          				=> cmic%cco2bios
	caerch4bios          			=> cmic%caerch4bios
	canaerch4bios          			=> cmic%canaerch4bios
	caces_unsat          			=> cmic%caces_unsat
	cacebios_unsat          		=> cmic%cacebios_unsat
	cco2bios_unsat          		=> cmic%cco2bios_unsat
	caerch4bios_unsat       		=> cmic%caerch4bios_unsat
	canaerch4bios_unsat    		=> cmic%canaerch4bios_unsat
	caces_sat          			=> cmic%caces_sat
	cacebios_sat          			=> cmic%cacebios_sat
	cco2bios_sat          			=> cmic%cco2bios_sat
	caerch4bios_sat          		=> cmic%caerch4bios_sat
	canaerch4bios_sat          		=> cmic%canaerch4bios_sat
	caces_prod          			=> cmic%caces_prod
	caces_unsat_prod          		=> cmic%caces_unsat_prod
	caces_sat_prod   			=> cmic%caces_sat_prod
  	caces_prod_h2          			=> cmic%caces_prod_h2
	caces_unsat_prod_h2          		=> cmic%caces_unsat_prod_h2
	caces_sat_prod_h2   			=> cmic%caces_sat_prod_h2
	
	dochr_vr          				=> cmic%dochr_vr
	micbio_hr_vr          			=> cmic%micbio_hr_vr
	micbio_hr          			=> cmic%micbio_hr
	root2doc          				=> cmic%root2doc
	froot_r          				=> cmic%froot_r
	cn_microbe          				=> cmic%cn_microbe
  
	ccon_ch4s          			=> cmic%ccon_ch4s
	ccon_co2s          			=> cmic%ccon_co2s
	ccon_o2s          				=> cmic%ccon_o2s
	ccon_h2s          				=> cmic%ccon_h2s
	ccon_ch4s_unsat          		=> cmic%ccon_ch4s_unsat
	ccon_co2s_unsat          		=> cmic%ccon_co2s_unsat
	ccon_o2s_unsat          		=> cmic%ccon_o2s_unsat
	ccon_h2s_unsat          		=> cmic%ccon_h2s_unsat
	ccon_ch4s_sat          			=> cmic%ccon_ch4s_sat
	ccon_co2s_sat          			=> cmic%ccon_co2s_sat
	ccon_o2s_sat          			=> cmic%ccon_o2s_sat
	ccon_h2s_sat          			=> cmic%ccon_h2s_sat
	finundated    				=> cws%finundated

  !sif           => cmic%sif
  !rootfr        => cps%pps_a%rootfr
  !o2stress_unsat   => cmic%o2stress_unsat
  !o2stress_sat     => cmic%o2stress_sat
  !ch4stress_unsat  => cmic%ch4stress_unsat
  !ch4stress_sat    => cmic%ch4stress_sat
  !totcolch4           => cmic%totcolch4

!------------------------------------------------------------------------
	if (masterproc) write (iulog,*) 'Attempting to initialize non-state variables for CH4 Mod'

	call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

  ! Assign local pointers to derived subtypes components (landunit-level)

	ltype           				=> lun%itype

  ! Assign local pointers to derived subtypes components (column-level)

	clandunit  					=> col%landunit
	cgridcell   					=> col%gridcell

	do c = begc,endc
	l = clandunit(c)

! First set levels from nlevsoi+1 to nlevgrnd = 0
	cmicbiocs(c,nlevsoi+1:nlevgrnd) 			= 0.
	cdocs_pre(c,nlevsoi+1:nlevgrnd) 			= 0.
	cdocs(c,nlevsoi+1:nlevgrnd) 				= 0.
	cdocs_unsat(c,nlevsoi+1:nlevgrnd) 			= 0.
	cdocs_sat(c,nlevsoi+1:nlevgrnd) 			= 0.
	cmicbions(c,nlevsoi+1:nlevgrnd) 			= 0.
	cdons(c,nlevsoi+1:nlevgrnd) 				= 0.
	cdons_unsat(c,nlevsoi+1:nlevgrnd) 			= 0.
	cdons_sat(c,nlevsoi+1:nlevgrnd) 			= 0.	
	cdons_min(c,nlevsoi+1:nlevgrnd) 			= 0.	
	caces(c,nlevsoi+1:nlevgrnd) 				= 0.
	cacebios(c,nlevsoi+1:nlevgrnd) 				= 0.
	cco2bios(c,nlevsoi+1:nlevgrnd) 				= 0.
	caerch4bios(c,nlevsoi+1:nlevgrnd) 			= 0.
	canaerch4bios(c,nlevsoi+1:nlevgrnd) 			= 0.
	caces_unsat(c,nlevsoi+1:nlevgrnd) 			= 0.
	cacebios_unsat(c,nlevsoi+1:nlevgrnd) 		= 0.
	cco2bios_unsat(c,nlevsoi+1:nlevgrnd) 		= 0.
	caerch4bios_unsat(c,nlevsoi+1:nlevgrnd) 		= 0.
	canaerch4bios_unsat(c,nlevsoi+1:nlevgrnd) 	= 0.
	caces_sat(c,nlevsoi+1:nlevgrnd) 			= 0.
	cacebios_sat(c,nlevsoi+1:nlevgrnd) 			= 0.
	cco2bios_sat(c,nlevsoi+1:nlevgrnd) 			= 0.
	caerch4bios_sat(c,nlevsoi+1:nlevgrnd) 		= 0.
	canaerch4bios_sat(c,nlevsoi+1:nlevgrnd) 		= 0.
	
	caces_prod(c,nlevsoi+1:nlevgrnd) 			= 0.
	caces_sat_prod(c,nlevsoi+1:nlevgrnd) 		= 0.
	caces_unsat_prod(c,nlevsoi+1:nlevgrnd) 		= 0.
	caces_prod_h2(c,nlevsoi+1:nlevgrnd) 			= 0.
	caces_sat_prod_h2(c,nlevsoi+1:nlevgrnd) 		= 0.
	caces_unsat_prod_h2(c,nlevsoi+1:nlevgrnd) 		= 0.
      
	dochr_vr(c,nlevsoi+1:nlevgrnd) 				= 0.
!	micbiohr_vr(c,nlevsoi+1:nlevgrnd) 			= 0.
	root2doc(c,nlevsoi+1:nlevgrnd) 				= 0.
	froot_r(c,nlevsoi+1:nlevgrnd) 				= 0.
	cn_microbe(c,nlevsoi+1:nlevgrnd) 			= 10.
	
	ccon_ch4s(c,nlevsoi+1:nlevgrnd) 			= 0.
	ccon_co2s(c,nlevsoi+1:nlevgrnd) 			= 0.
	ccon_o2s(c,nlevsoi+1:nlevgrnd) 			= 0.
	ccon_h2s(c,nlevsoi+1:nlevgrnd) 			= 0.
	ccon_ch4s_unsat(c,nlevsoi+1:nlevgrnd) 		= 0.
	ccon_co2s_unsat(c,nlevsoi+1:nlevgrnd) 		= 0.
	ccon_o2s_unsat(c,nlevsoi+1:nlevgrnd) 		= 0.
	ccon_h2s_unsat(c,nlevsoi+1:nlevgrnd) 		= 0.
	ccon_ch4s_sat(c,nlevsoi+1:nlevgrnd) 		= 0.
	ccon_co2s_sat(c,nlevsoi+1:nlevgrnd) 		= 0.
	ccon_o2s_sat(c,nlevsoi+1:nlevgrnd) 			= 0.
	ccon_h2s_sat(c,nlevsoi+1:nlevgrnd) 			= 0.
	end do
	if (masterproc) write (iulog,*) 'Successfully initialized non-state variables for CH4 component of microbial Mod'
end subroutine initTimeConst_microbe

#endif
!end if microbe
end module initmicrobeMod
