
module CNCStateUpdate4Mod
#ifdef MICROBE

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: CStateUpdate3Mod
!
! !DESCRIPTION:
! Module for carbon state variable update, mortality fluxes.
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    implicit none
    save
    private
! !PUBLIC MEMBER FUNCTIONS:
    public:: CStateUpdate4
!
! !REVISION HISTORY:
! Xiaofeng Xu (4/1/2020) created for methane module
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CStateUpdate4
!
! !INTERFACE:
subroutine CStateUpdate4(num_soilc, filter_soilc, num_soilp, filter_soilp, isotope)
!
! !DESCRIPTION:
! On the radiation time step, update all the prognostic carbon state
! variables affected by fire fluxes
!
! !USES:

! !USES:
	use shr_kind_mod , only : r8 => shr_kind_r8
	use clmtype
	use subgridAveMod , only : p2c, c2g
	use clm_varpar, only : nlevgrnd, nlevdecomp, i_bacteria, i_fungi, i_dom, cn_dom
	use pftvarcon,  only : noveg
	use clm_varcon, only : secspday, istwet, istsoil, istdlak, spval, istcrop
	use microbevarcon 
	use clm_time_manager, only : get_step_size, get_nstep
	use filterMod,  only : filter
	use shr_const_mod, only: SHR_CONST_TKFRZ, SHR_CONST_RGAS
        use abortutils  , only: endrun

! !ARGUMENTS:
   implicit none
   integer, intent(in) :: num_soilc       ! number of soil columns in filter
   integer, intent(in) :: filter_soilc(:) ! filter for soil columns
   integer, intent(in) :: num_soilp       ! number of soil pfts in filter
   integer, intent(in) :: filter_soilp(:) ! filter for soil pfts
   character(len=*), intent(in) :: isotope         ! 'bulk', 'c13' or 'c14'
!
! !LOCAL VARIABLES:
! local pointers to implicit in variables
	real(r8), pointer :: decomp_cpools_vr(:,:,:) 		! bring DOMC to methane module
	real(r8), pointer :: cdocs(:,:)					! column-level concentration of DOC molC/m3 
	
	real(r8), pointer :: cmicbiocs(:,:)				! column-level biomass of all microbes molC/m3
	real(r8), pointer :: caces(:,:)					! column-level concentration of acetate molC/m3
	real(r8), pointer :: cacebios(:,:)				! column-level biomass of methanogen based on acetate molC/m3
	real(r8), pointer :: cco2bios(:,:)				! column-level biomass of methanogen based on CO2/H2 molC/m3
	real(r8), pointer :: caerch4bios(:,:)			! column-level biomass of aerobix methanotrophy molC/m3
	real(r8), pointer :: canaerch4bios(:,:)			! column-level biomass of anaerobic methanotrophy molC/m3

	real(r8), pointer :: caces_prod(:,:)				! acetate acid production
	real(r8), pointer :: caces_prod_h2(:,:)			! acetogenesis
	
	real(r8), pointer :: ccon_ch4s(:,:)				! concentration of CH4 mol C/m3 at each layer
	real(r8), pointer :: ccon_co2s(:,:)				! column-level concentration of CO2 mol C/m3
      
	real(r8), pointer :: ch4_prod_ace_depth(:,:)		! column-level aceteclastic methanogenesis
	real(r8), pointer :: ch4_prod_co2_depth(:,:)		! column-level hydrogenotrophic methanogenesis
	real(r8), pointer :: ch4_oxid_o2_depth(:,:)		! column-level oxic methanotrophy
	real(r8), pointer :: ch4_oxid_aom_depth(:,:)		! column-level anaerobic oxidation of methane
	
	real(r8), pointer :: ch4_aere_depth(:,:)			! column-level 
	real(r8), pointer :: ch4_dif_depth(:,:)
	real(r8), pointer :: ch4_ebul_depth(:,:)
	
	real(r8), pointer :: co2_prod_ace_depth(:,:)
	real(r8), pointer :: co2_decomp_depth(:,:)
	real(r8), pointer :: co2_cons_depth(:,:)
	
	real(r8), pointer :: co2_ebul_depth(:,:)
	real(r8), pointer :: co2_aere_depth(:,:)
	real(r8), pointer :: co2_dif_depth(:,:)
	
	real(r8), pointer :: ch4_surf_aere(:)
	real(r8), pointer :: ch4_surf_ebul(:)
	real(r8), pointer :: ch4_surf_dif(:)
	real(r8), pointer :: ch4_surf_netflux(:)
	
	real(r8), pointer :: co2_surf_aere(:)
	real(r8), pointer :: co2_surf_ebul(:)
	real(r8), pointer :: co2_surf_dif(:)
	real(r8), pointer :: co2_surf_netflux(:)

! !OTHER LOCAL VARIABLES:
   type(column_microbe_type), pointer :: cmiciso
    type(column_cstate_type), pointer :: ccisos
   integer :: c,p,j,l,k,fc      ! indices
   real(r8):: dt       ! radiation time step (seconds)

!EOP
!-----------------------------------------------------------------------
   ! select which isotope
   select case (isotope)
   case ('bulk')
      ccisos =>  ccs
      cmiciso =>  cmic
   case ('c14')
      ccisos =>  cc14s
      cmiciso =>  cmicc14
   case ('c13')
      ccisos =>  cc13s
      cmiciso =>  cmicc13
   case default
      call endrun('CNCIsoStateUpdate3Mod: iso must be bulk, c13 or c14')
   end select

    
! !ARGUMENTS:

   ! Column level pointers
	decomp_cpools_vr			=> ccs%decomp_cpools_vr
	cdocs 					=> cmiciso%cdocs
	
	cmicbiocs 					=> cmiciso%cmicbiocs
	caces 					=> cmiciso%caces
	caces_prod 				=> cmiciso%caces_prod
	caces_prod_h2 				=> cmiciso%caces_prod_h2

	cacebios 					=> cmiciso%cacebios
	cco2bios 					=> cmiciso%cco2bios
	caerch4bios 				=> cmiciso%caerch4bios
	canaerch4bios 				=> cmiciso%canaerch4bios
     
	ccon_ch4s        				=> cmiciso%ccon_ch4s
	ccon_co2s         			=> cmiciso%ccon_co2s
 
	ch4_prod_ace_depth 			=> cmiciso%ch4_prod_ace_depth
	ch4_prod_co2_depth 			=> cmiciso%ch4_prod_co2_depth
	ch4_oxid_o2_depth 			=> cmiciso%ch4_oxid_o2_depth
	ch4_oxid_aom_depth 		=> cmiciso%ch4_oxid_aom_depth
	
	ch4_aere_depth 			=> cmiciso%ch4_aere_depth
	ch4_dif_depth 				=> cmiciso%ch4_dif_depth
	ch4_ebul_depth 			=> cmiciso%ch4_ebul_depth
	
	co2_prod_ace_depth 			=> cmiciso%co2_prod_ace_depth
	co2_decomp_depth 			=> cmiciso%co2_decomp_depth
	co2_cons_depth 			=> cmiciso%co2_cons_depth
	
	co2_ebul_depth 			=> cmiciso%co2_ebul_depth
	co2_aere_depth 			=> cmiciso%co2_aere_depth
	co2_dif_depth 				=> cmiciso%co2_dif_depth
	
	ch4_surf_aere				=> cmiciso%ch4_surf_aere
	ch4_surf_ebul 				=> cmiciso%ch4_surf_ebul
	ch4_surf_dif 				=> cmiciso%ch4_surf_dif
	ch4_surf_netflux 			=> cmiciso%ch4_surf_netflux
	
	co2_surf_aere				=> cmiciso%co2_surf_aere
	co2_surf_ebul 				=> cmiciso%co2_surf_ebul
	co2_surf_dif 				=> cmiciso%co2_surf_dif
	co2_surf_netflux 			=> cmiciso%co2_surf_netflux

    ! set time steps
    dt = real( get_step_size(), r8 )
    
	do j = 1,nlevdecomp
	do fc = 1,num_soilc
		c = filter_soilc(fc)
		cdocs(c,j) = ccisos%decomp_cpools_vr(c,j,i_dom)
		
		caces(c,j) = max(0._r8, (caces(c,j) + caces_prod(c,j)*dt + caces_prod_h2(c,j)*dt - ch4_prod_ace_depth(c,j)*dt))
!		cacebios(c,j) = cacebios(c,j)
!		cco2bios(c,j) = cco2bios(c,j)
!		caerch4bios(c,j) = caerch4bios(c,j)
!		canaerch4bios(c,j) = canaerch4bios(c,j)
		ccon_ch4s(c,j)=max(0._r8, (ccon_ch4s(c,j) + ch4_prod_ace_depth(c,j)*dt + ch4_prod_co2_depth(c,j)*dt - &
		ch4_oxid_o2_depth(c,j)*dt - ch4_oxid_aom_depth(c,j)*dt - &
		ch4_aere_depth(c,j)*dt - ch4_dif_depth(c,j)*dt - ch4_ebul_depth(c,j)*dt))
	end do
	end do

end subroutine CStateUpdate4
!-----------------------------------------------------------------------
#endif

end module CNCStateUpdate4Mod
