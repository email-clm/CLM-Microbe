module CNDecompMod
#ifdef CN

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: CNDecompMod
!
! !DESCRIPTION:
! Module holding routines used in litter and soil decomposition model
! for coupled carbon-nitrogen code.
!
! !USES:
   use shr_kind_mod , only: r8 => shr_kind_r8
   use shr_const_mod, only: SHR_CONST_TKFRZ
    use clm_varctl  , only: iulog
    use clm_varcon, only: dzsoi_decomp, zsoi
#ifdef MICROBE
    use clm_varcon, only: MBC_k

#endif
#ifndef CENTURY_DECOMP
    use CNDecompCascadeMod_BGC, only : decomp_rate_constants
#else
    use CNDecompCascadeMod_CENTURY, only : decomp_rate_constants
#endif
#ifdef NITRIF_DENITRIF
    use CNNitrifDenitrifMod, only: nitrif_denitrif
#endif
    use CNVerticalProfileMod, only: decomp_vertprofiles

   implicit none
   save
   private
! !PUBLIC MEMBER FUNCTIONS:
   public:: CNDecompAlloc

!#ifdef MICROBE
!   real(r8), public :: decomp_depth_efolding = 0.17_r8 !0.5    ! (meters) e-folding depth for reduction in decomposition [set to large number for depth-independance]
!#endif

!
! !REVISION HISTORY:
! 8/15/03: Created by Peter Thornton
! 10/23/03, Peter Thornton: migrated to vector data structures
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNDecompAlloc
!
! !INTERFACE:
subroutine CNDecompAlloc (lbp, ubp, lbc, ubc, num_soilc, filter_soilc, &
   num_soilp, filter_soilp)
!
! !DESCRIPTION:
!
! !USES:
   use clmtype
   use CNAllocationMod , only: CNAllocation
   use clm_time_manager, only: get_step_size
   use clm_varpar   , only: nlevsoi,nlevgrnd,nlevdecomp,ndecomp_cascade_transitions,ndecomp_pools
   use pft2colMod      , only: p2c
 
#ifdef MICROBE
   USE microbevarcon , only: plant2doc, som2doc, micbiocn, micbioMR, CUEref, CUEt, Tcueref 
   USE microbevarcon , only: Tsref, Tmref, Tsmin, Tmmin, Msmin, Msmax, Mmmin, dock 
   use CNDecompCascadeMod_BGC
   use clm_time_manager, only : get_step_size
   use clm_varcon, only: secspday   
   use clm_varpar, only: i_bacteria, i_fungi, i_dom, cn_dom, numpft
   use clm_varpar      , only: max_pft_per_col

#endif
!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: lbp, ubp        ! pft-index bounds
   integer, intent(in) :: lbc, ubc        ! column-index bounds
   integer, intent(in) :: num_soilc       ! number of soil columns in filter
   integer, intent(in) :: filter_soilc(ubc-lbc+1) ! filter for soil columns
   integer, intent(in) :: num_soilp       ! number of soil pfts in filter
   integer, intent(in) :: filter_soilp(ubp-lbp+1) ! filter for soil pfts
!
! !CALLED FROM:
! subroutine CNEcosystemDyn in module CNEcosystemDynMod.F90
!
! !REVISION HISTORY:
! 8/15/03: Created by Peter Thornton
!
! !LOCAL VARIABLES:
! local pointers to implicit in scalars
!
   ! all c pools involved in decomposition
   real(r8), pointer :: decomp_cpools_vr(:,:,:)    ! (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) c pools
   ! all n pools involved in decomposition
   real(r8), pointer :: decomp_npools_vr(:,:,:)    ! (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools

   integer, pointer :: clandunit(:)      			! index into landunit level quantities
   integer , pointer :: itypelun(:)      			! landunit type
   ! pft level
   real(r8), pointer :: rootfr(:,:)      			! fraction of roots in each soil layer  (nlevgrnd)
!
! local pointers to implicit in/out scalars
!
   real(r8), pointer :: fpi_vr(:,:)                            ! fraction of potential immobilization (no units)
   real(r8), pointer :: decomp_cascade_hr_vr(:,:,:)            ! vertically-resolved het. resp. from decomposing C pools (gC/m3/s)
   real(r8), pointer :: decomp_cascade_ctransfer_vr(:,:,:)     ! vertically-resolved het. resp. from decomposing C pools (gC/m3/s)
   real(r8), pointer :: decomp_pools_hr(:,:)                   ! het. resp. from decomposing C pools (gC/m2/s)
   real(r8), pointer :: decomp_cascade_ctransfer(:,:)          ! vertically-resolved C transferred along deomposition cascade (gC/m3/s)
   real(r8), pointer :: decomp_cascade_ntransfer_vr(:,:,:)     ! vert-res transfer of N from donor to receiver pool along decomp. cascade (gN/m3/s)
   real(r8), pointer :: decomp_cascade_sminn_flux_vr(:,:,:)    ! vert-res mineral N flux for transition along decomposition cascade (gN/m3/s)

   real(r8), pointer :: potential_immob_vr(:,:)
#ifndef NITRIF_DENITRIF
   real(r8), pointer :: sminn_to_denit_decomp_cascade_vr(:,:,:)
   real(r8), pointer :: sminn_to_denit_excess_vr(:,:)
#endif
   real(r8), pointer :: gross_nmin_vr(:,:)
   real(r8), pointer :: net_nmin_vr(:,:)
   real(r8), pointer :: gross_nmin(:)            ! gross rate of N mineralization (gN/m2/s)
   real(r8), pointer :: net_nmin(:)              ! net rate of N mineralization (gN/m2/s)
   ! For methane code
#ifdef LCH4
   real(r8), pointer :: fphr(:,:)                ! fraction of potential SOM + LITTER heterotrophic respiration
   real(r8), pointer :: w_scalar(:,:)            ! fraction by which decomposition is limited by moisture availability
#endif

#ifdef MICROBE
!   real(r8), pointer :: fphr(:,:)                ! fraction of potential SOM + LITTER heterotrophic respiration
   real(r8), pointer :: w_scalar(:,:)            ! fraction by which decomposition is limited by moisture availability
   real(r8), pointer :: t_scalar(:,:)            ! fraction by which decomposition is limited by temperature
   real(r8), pointer :: o_scalar(:,:)            ! fraction by which decomposition is limited by oxygen
!   real(r8), pointer :: depth_scalar(:,:)            ! fraction by which decomposition is reduced along soil profile
   real(r8), pointer :: cdocs_pre(:,:)                 ! the available carbon for methane processes
   real(r8), pointer :: cdocs(:,:)                 ! the available carbon for methane processes
   real(r8), pointer :: cdocs_unsat(:,:)             ! the available carbon for methane processes in unsaturated fraction
   real(r8), pointer :: cdocs_sat(:,:)                 ! the available carbon for methane processes in saturated fraction
   real(r8), pointer :: cmicbiocs(:,:)                 ! microbial biomass carbon at columne-level
   real(r8), pointer :: cdons(:,:)                 	! availabe nitrogen
   real(r8), pointer :: cdons_unsat(:,:)             ! nitrogen
   real(r8), pointer :: cdons_sat(:,:)                 ! nitrogen
   real(r8), pointer :: cdons_min(:,:)                 ! nitrogen
   real(r8), pointer :: cmicbions(:,:)                 ! microbial biomass nitrogen at columne-level
   real(r8), pointer :: h2osoi_vol(:,:)                 ! water volume for each column
   real(r8), pointer :: watsat(:,:)                 ! the available carbon for methane processes
   real(r8), pointer :: t_soisno(:,:)                 ! the available carbon for methane processes   
   real(r8), pointer :: dochr_vr(:,:)                 ! vertically resolved decomposition of dissolved organic carbon 
   real(r8), pointer :: micbio_hr_vr(:,:,:)                 ! vertically resolved decomposition of dissolved organic carbon 
   real(r8), pointer :: micbio_hr(:,:)                 ! vertically resolved decomposition of dissolved organic carbon 
   real(r8), pointer :: finundated(:)			!inundated fractin of soil column
   real(r8), Pointer :: froot_r(:,:)
   real(r8), Pointer :: root2doc(:,:)			! part of fine root exudate carbon to soil to develop dissolved organic carbon
									! we set it as part of fine root respiration
   real(r8), Pointer :: cn_microbe(:,:)			! C:N ration of the microbe as a combination of bacteria and fungi
   logical, pointer :: is_microbe(:)                            ! TRUE => pool is a microbe pool

   real(r8),allocatable :: decomp_depth_efolding_in(:) ! define the active decomposing soil depth
   real(r8), pointer :: decomp_depth_efolding(:) ! define the active decomposing soil depth

   integer , pointer :: ivt(:)             ! pft vegetation type
   integer :: pi, p

    real(r8), pointer :: wtcol(:)                ! pft weight relative to column (0-1)
    integer, allocatable :: pft_index(:)
    integer , pointer :: pfti(:)        ! pft index array
    integer , pointer :: npfts(:)       ! number of pfts on the column

#endif

   logical, pointer :: is_litter(:)                         ! TRUE => pool is a litter pool
   logical, pointer :: is_soil(:)                           ! TRUE => pool is a soil pool

   real(r8), pointer :: decomp_k(:,:,:)                       ! rate constant for decomposition (1./sec)
   real(r8), pointer :: rf_decomp_cascade(:,:,:)              ! respired fraction in decomposition step (frac)
   integer,  pointer :: cascade_donor_pool(:)                 ! which pool is C taken from for a given decomposition step
   integer,  pointer :: cascade_receiver_pool(:)              ! which pool is C added to for a given decomposition step
   real(r8), pointer :: pathfrac_decomp_cascade(:,:,:)        ! what fraction of C leaving a given pool passes through a given transition (frac)
   logical,  pointer :: floating_cn_ratio_decomp_pools(:)     ! TRUE => pool has fixed C:N ratio

!
! local pointers to implicit out scalars
!
! !OTHER LOCAL VARIABLES:
   integer :: c,j,k,l,m          !indices
   integer :: fc           !lake filter column index
   real(r8):: p_decomp_cpool_loss(lbc:ubc,1:nlevdecomp,1:ndecomp_cascade_transitions) !potential C loss from one pool to another
   real(r8):: pmnf_decomp_cascade(lbc:ubc,1:nlevdecomp,1:ndecomp_cascade_transitions) !potential mineral N flux, from one pool to another
   real(r8):: immob(lbc:ubc,1:nlevdecomp)        !potential N immobilization
   real(r8):: ratio        !temporary variable
   real(r8):: dnp          !denitrification proportion
   real(r8):: cn_decomp_pools(lbc:ubc,1:nlevdecomp,1:ndecomp_pools)
   real(r8), pointer :: initial_cn_ratio(:)               ! c:n ratio for initialization of pools
   integer, parameter :: i_atm = 0
   integer, pointer :: altmax_indx(:)                  ! maximum annual depth of thaw
   integer, pointer :: altmax_lastyear_indx(:)         ! prior year maximum annual depth of thaw

#ifdef MICROBE
   real(r8) :: CUE, fm_t,  fm_m
   real(r8) :: depth_scalar(lbc:ubc,1:nlevdecomp) 
   real(r8) :: microbeMR
   real(r8) :: doc_k
   integer :: dt, dtd
#endif

   ! For methane code
#ifndef NITRIF_DENITRIF
   real(r8):: phr_vr(lbc:ubc,1:nlevdecomp)       !potential HR (gC/m3/s)
#else
   real(r8), pointer :: phr_vr(:,:)              !potential HR (gC/m3/s)
#endif
   real(r8):: hrsum(lbc:ubc,1:nlevdecomp)        !sum of HR (gC/m2/s)
   
   !EOP
   !-----------------------------------------------------------------------
   
   decomp_cpools_vr              		=> ccs%decomp_cpools_vr
   decomp_cascade_hr_vr         		=> ccf%decomp_cascade_hr_vr
   decomp_cascade_ctransfer_vr   	=> ccf%decomp_cascade_ctransfer_vr
   decomp_npools_vr              		=> cns%decomp_npools_vr
   decomp_cascade_ntransfer_vr   	=> cnf%decomp_cascade_ntransfer_vr
   decomp_cascade_sminn_flux_vr  	=> cnf%decomp_cascade_sminn_flux_vr
   fpi_vr                				=> cps%fpi_vr
   potential_immob_vr    			=> cnf%potential_immob_vr
   
   decomp_k                        		=> ccf%decomp_k
   rf_decomp_cascade               		=> cps%rf_decomp_cascade
   cascade_donor_pool              		=> decomp_cascade_con%cascade_donor_pool
   cascade_receiver_pool           		=> decomp_cascade_con%cascade_receiver_pool
   pathfrac_decomp_cascade         	=> cps%pathfrac_decomp_cascade
   floating_cn_ratio_decomp_pools  	=> decomp_cascade_con%floating_cn_ratio_decomp_pools
   initial_cn_ratio                		=> decomp_cascade_con%initial_cn_ratio
   altmax_indx                     		=> cps%altmax_indx
   altmax_lastyear_indx            		=> cps%altmax_lastyear_indx
   
#ifndef NITRIF_DENITRIF
   sminn_to_denit_decomp_cascade_vr => cnf%sminn_to_denit_decomp_cascade_vr
   sminn_to_denit_excess_vr 		=> cnf%sminn_to_denit_excess_vr
#else
   phr_vr                   				=> ccf%phr_vr
#endif
   gross_nmin_vr            			=> cnf%gross_nmin_vr
   net_nmin_vr              			=> cnf%net_nmin_vr
   gross_nmin               			=> cnf%gross_nmin
   net_nmin                 			=> cnf%net_nmin
   ! For methane code
#ifdef LCH4
   fphr               				=> cch4%fphr
   w_scalar           				=> ccf%w_scalar
#endif
   
   rootfr                				=> pps%rootfr
   clandunit             				=>col%landunit
   itypelun              				=> lun%itype
   
#ifdef MICROBE
   w_scalar           				=> ccf%w_scalar
   t_scalar           				=> ccf%t_scalar
   o_scalar           				=> ccf%o_scalar
!   depth_scalar           => ccf%depth_scalar
   cdocs          					=> cmic%cdocs
   cdocs_unsat       				=> cmic%cdocs_unsat
   cdocs_sat           				=> cmic%cdocs_sat
   cmicbiocs           				=> cmic%cmicbiocs
   cdons           					=> cmic%cdons
   cdons_unsat       				=> cmic%cdons_unsat
   cdons_sat           				=> cmic%cdons_sat
   cdons_min           				=> cmic%cdons_min
   cmicbions           				=> cmic%cmicbions
   h2osoi_vol           				=> cws%h2osoi_vol
   watsat           					=> cps%watsat
   t_soisno           				=> ces%t_soisno
   dochr_vr           				=> cmic%dochr_vr
   micbio_hr_vr           			=> cmic%micbio_hr_vr
   micbio_hr           				=> cmic%micbio_hr
   root2doc           				=> cmic%root2doc
   froot_r						=> cmic%froot_r
   cn_microbe           				=> cmic%cn_microbe
   finundated					=> cws%finundated   
   is_microbe                                  	=> decomp_cascade_con%is_microbe

   decomp_depth_efolding =>pftcon%decomp_depth_efolding 
   ivt                           =>pft%itype
  
   wtcol          =>pft%wtcol
   pfti             =>col%pfti
   npfts            =>col%npfts

  if (num_soilc .gt. 0) then

  allocate(pft_index(0))

  do fc = 1,num_soilc
  c = filter_soilc(fc)

  do pi = 1,max_pft_per_col

        if (pi <=  npfts(c)) then
           p = pfti(c) + pi - 1
        end if
 
  end do

   pft_index = [pft_index, (MAXLOC(wtcol(pfti(c):p), DIM=1, mask = wtcol(pfti(c):p) .gt. 0) - 1)]

  end do
   
  end if

decomp_depth_efolding_in = decomp_depth_efolding(pft_index(:))


#endif  

   is_litter                               		=> decomp_cascade_con%is_litter
   is_soil                                 		=> decomp_cascade_con%is_soil

   call decomp_rate_constants(lbc, ubc, num_soilc, filter_soilc)
   
   ! set initial values for potential C and N fluxes
   p_decomp_cpool_loss(:,:,:) = 0._r8
   pmnf_decomp_cascade(:,:,:) = 0._r8
   
   ! column loop to calculate potential decomp rates and total immobilization
   ! demand.
   
   !!! calculate c:n ratios of applicable pools
   do l = 1, ndecomp_pools
      if ( floating_cn_ratio_decomp_pools(l) ) then
         do j = 1,nlevdecomp
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               if ( decomp_npools_vr(c,j,l) .gt. 0._r8 ) then
                  cn_decomp_pools(c,j,l) = decomp_cpools_vr(c,j,l) / decomp_npools_vr(c,j,l)
               end if
            end do
         end do
      else
         do j = 1,nlevdecomp
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               cn_decomp_pools(c,j,l) = initial_cn_ratio(l)
            end do
         end do
      end if
   end do
   
#ifdef MICROBE
      ! add a term to reduce decomposition rate at depth
      ! for now used a fixed e-folding depth
      do j = 1, nlevdecomp
         do fc = 1,num_soilc
            c = filter_soilc(fc)
	    if((decomp_npools_vr(c,j,i_bacteria) + decomp_npools_vr(c,j,i_fungi)) < 1e-15) then
	    cn_microbe(c,j) = 10._r8
	    else
	    cn_microbe(c,j) = (decomp_cpools_vr(c,j,i_bacteria) + decomp_cpools_vr(c,j,i_fungi)) / (decomp_npools_vr(c,j,i_bacteria) + decomp_npools_vr(c,j,i_fungi))
	    end if
            depth_scalar(c,j) = exp(-zsoi(j)/decomp_depth_efolding_in(fc))
	    
	    cdocs(c,j) = max(cdocs(c,j),1e-20)
!write(*,*) "here debugging1: ", c,j, cdocs(c,j), cdons(c,j), decomp_cpools_vr(c,j,i_dom),decomp_npools_vr(c,j,i_dom)	    
	    !decomp_cpools_vr(c,j,i_dom) = cdocs(c,j)
	    !decomp_npools_vr(c,j,i_dom) = cdocs(c,j) / cn_dom
!write(*,*) "here debugging: ", c,j, cdocs(c,j), cdons(c,j)
         end do
      end do

   dt = real( get_step_size(), r8 )
   dtd = dt/secspday
   
   !~ microbeMR = -log(1.0_r8-micbioMR)   
   !~ doc_k = -log(1.0_r8-dock)  
   !~ ! calculate the new discrete-time decay rate for model timestep
   !~ microbeMR = 1.0_r8-exp(-microbeMR*dtd)   
   !~ doc_k = 1.0_r8-exp(-doc_k*dtd)  
#endif   
   
   ! calculate the non-nitrogen-limited fluxes
   ! these fluxes include the  "/ dt" term to put them on a
   ! per second basis, since the rate constants have been
   ! calculated on a per timestep basis.
   
   do k = 1, ndecomp_cascade_transitions
      do j = 1,nlevdecomp
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            
            if (decomp_cpools_vr(c,j,cascade_donor_pool(k)) > 0._r8 .and. decomp_k(c,j,cascade_donor_pool(k)) .gt. 0._r8 ) then
#ifdef MICROBE
if ((is_litter(cascade_donor_pool(k)) .or. is_soil(cascade_donor_pool(k))) .and. (is_microbe(cascade_donor_pool(k)) .eqv. .false.)) then
		p_decomp_cpool_loss(c,j,k) = decomp_cpools_vr(c,j,cascade_donor_pool(k)) * decomp_k(c,j,cascade_donor_pool(k))  * pathfrac_decomp_cascade(c,j,k) !* (cmicbiocs(c,j) + 1e-5) / (cmicbiocs(c,j) + MBC_k(cascade_donor_pool(k)))
else
		p_decomp_cpool_loss(c,j,k) = decomp_cpools_vr(c,j,cascade_donor_pool(k)) * decomp_k(c,j,cascade_donor_pool(k))  * pathfrac_decomp_cascade(c,j,k)
end if
#else
                p_decomp_cpool_loss(c,j,k) = decomp_cpools_vr(c,j,cascade_donor_pool(k)) * decomp_k(c,j,cascade_donor_pool(k))  * pathfrac_decomp_cascade(c,j,k)
#endif
               if ( .not. floating_cn_ratio_decomp_pools(cascade_receiver_pool(k)) ) then  !! not transition of cwd to litter
                  
                  if (cascade_receiver_pool(k) .ne. i_atm ) then  ! not 100% respiration
                     ratio = 0._r8
                     
                     if (decomp_npools_vr(c,j,cascade_donor_pool(k)) > 0._r8) then
                        ratio = cn_decomp_pools(c,j,cascade_receiver_pool(k))/cn_decomp_pools(c,j,cascade_donor_pool(k))
                     endif
                     
                     pmnf_decomp_cascade(c,j,k) = (p_decomp_cpool_loss(c,j,k) * (1.0_r8 - rf_decomp_cascade(c,j,k) - ratio) &
                          / cn_decomp_pools(c,j,cascade_receiver_pool(k)) )
                     
                  else   ! 100% respiration
                     pmnf_decomp_cascade(c,j,k) = - p_decomp_cpool_loss(c,j,k) / cn_decomp_pools(c,j,cascade_donor_pool(k))
                  endif
                  
               else   ! CWD -> litter
                  pmnf_decomp_cascade(c,j,k) = 0._r8
               end if
            end if
         end do
         
      end do
   end do
   
   ! Sum up all the potential immobilization fluxes (positive pmnf flux)
   ! and all the mineralization fluxes (negative pmnf flux)
   do j = 1,nlevdecomp
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         immob(c,j) = 0._r8
      end do
   end do
   do k = 1, ndecomp_cascade_transitions
      do j = 1,nlevdecomp
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            if (pmnf_decomp_cascade(c,j,k) > 0._r8) then
               immob(c,j) = immob(c,j) + pmnf_decomp_cascade(c,j,k)
            else
               gross_nmin_vr(c,j) = gross_nmin_vr(c,j) - pmnf_decomp_cascade(c,j,k)
            end if
         end do
      end do
   end do
   
   do j = 1,nlevdecomp
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         potential_immob_vr(c,j) = immob(c,j)
      end do
   end do
   
   ! Add up potential hr for methane calculations
   do j = 1,nlevdecomp
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         phr_vr(c,j) = 0._r8
      end do
   end do
   do k = 1, ndecomp_cascade_transitions
      do j = 1,nlevdecomp
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            phr_vr(c,j) = phr_vr(c,j) + rf_decomp_cascade(c,j,k) * p_decomp_cpool_loss(c,j,k)
         end do
      end do
   end do
   
   call decomp_vertprofiles(lbp, ubp, lbc,ubc,num_soilc,filter_soilc,num_soilp,filter_soilp)
   
#ifdef NITRIF_DENITRIF
   ! calculate nitrification and denitrification rates
   call nitrif_denitrif(lbc, ubc, num_soilc, filter_soilc)
#endif
   
   
   ! now that potential N immobilization is known, call allocation
   ! to resolve the competition between plants and soil heterotrophs
   ! for available soil mineral N resource.
   
   call CNAllocation(lbp, ubp, lbc,ubc,num_soilc,filter_soilc,num_soilp, &
        filter_soilp)
   
   ! column loop to calculate actual immobilization and decomp rates, following
   ! resolution of plant/heterotroph  competition for mineral N
   
   dnp = 0.01_r8
   
   ! calculate c:n ratios of applicable pools
   do l = 1, ndecomp_pools
      if ( floating_cn_ratio_decomp_pools(l) ) then
         do j = 1,nlevdecomp
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               if ( decomp_npools_vr(c,j,l) .gt. 0._r8 ) then
                  cn_decomp_pools(c,j,l) = decomp_cpools_vr(c,j,l) / decomp_npools_vr(c,j,l)
               end if
            end do
         end do
      else
         do j = 1,nlevdecomp
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               cn_decomp_pools(c,j,l) = initial_cn_ratio(l)
            end do
         end do
      end if
   end do
   
   ! upon return from CNAllocation, the fraction of potential immobilization
   ! has been set (cps%fpi_vr). now finish the decomp calculations.
   ! Only the immobilization steps are limited by fpi_vr (pmnf > 0)
   ! Also calculate denitrification losses as a simple proportion
   ! of mineralization flux.
   
   do k = 1, ndecomp_cascade_transitions
      do j = 1,nlevdecomp
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            
            if (decomp_cpools_vr(c,j,cascade_donor_pool(k)) .gt. 0._r8) then
               if ( pmnf_decomp_cascade(c,j,k) .gt. 0._r8 ) then
                  p_decomp_cpool_loss(c,j,k) = p_decomp_cpool_loss(c,j,k) * fpi_vr(c,j)
                  pmnf_decomp_cascade(c,j,k) = pmnf_decomp_cascade(c,j,k) * fpi_vr(c,j)
#ifndef NITRIF_DENITRIF
                  sminn_to_denit_decomp_cascade_vr(c,j,k) = 0._r8
               else
                  sminn_to_denit_decomp_cascade_vr(c,j,k) = -dnp * pmnf_decomp_cascade(c,j,k)
#endif
               end if
               decomp_cascade_hr_vr(c,j,k) = rf_decomp_cascade(c,j,k) * p_decomp_cpool_loss(c,j,k)
               decomp_cascade_ctransfer_vr(c,j,k) = (1._r8 - rf_decomp_cascade(c,j,k)) * p_decomp_cpool_loss(c,j,k)
                              
               if (decomp_npools_vr(c,j,cascade_donor_pool(k)) > 0._r8 .and. cascade_receiver_pool(k) .ne. i_atm) then
                  decomp_cascade_ntransfer_vr(c,j,k) = p_decomp_cpool_loss(c,j,k) / cn_decomp_pools(c,j,cascade_donor_pool(k))
               else
                  decomp_cascade_ntransfer_vr(c,j,k) = 0._r8
               endif
               if ( cascade_receiver_pool(k) .ne. 0 ) then
                  decomp_cascade_sminn_flux_vr(c,j,k) = pmnf_decomp_cascade(c,j,k)
               else  ! keep sign convention negative for terminal pools
                  decomp_cascade_sminn_flux_vr(c,j,k) = - pmnf_decomp_cascade(c,j,k)
               endif
               net_nmin_vr(c,j) = net_nmin_vr(c,j) - pmnf_decomp_cascade(c,j,k)
            else
               decomp_cascade_ntransfer_vr(c,j,k) = 0._r8
#ifndef NITRIF_DENITRIF
               sminn_to_denit_decomp_cascade_vr(c,j,k) = 0._r8
#endif
               decomp_cascade_sminn_flux_vr(c,j,k) = 0._r8
            end if
            
!~ #ifdef MICROBE    
	!~ if(decomp_cpools_vr(c,j,cascade_donor_pool(k)) .gt. 0._r8)
	!~ decomp_cascade_hr_vr(c,j,k) = rf_decomp_cascade(c,j,k) * p_decomp_cpool_loss(c,j,k)
!~ #endif

         end do
      end do
   end do
   
   
#ifdef MICROBE
!   do l = 1, ndecomp_pools
       do j = 1,nlevdecomp
            do fc = 1,num_soilc
               c = filter_soilc(fc)
      !~ decomp_cpools_vr(c,j,i_dom) = decomp_cpools_vr(c,j,i_dom) !- cdons_min(c,j) * cn_dom
!write(iulog, *) "xiaofeng here1 ", cdocs(c,j), decomp_cpools_vr(c,j,i_dom)  
!write(*,*) "here debugging2: ", c,j, cdocs(c,j), cdons(c,j), decomp_cpools_vr(c,j,i_dom),decomp_npools_vr(c,j,i_dom)	    
      !~ cdocs(c,j) = decomp_cpools_vr(c,j,i_dom)   ! gC/m3
!write(iulog, *) "xiaofeng here2 ", cdocs(c,j), decomp_cpools_vr(c,j,i_dom)     
      !~ cdocs_unsat(c,j) = cdocs(c,j) * (1. - finundated(c))
      !~ cdocs_sat(c,j) = cdocs(c,j) * finundated(c)
      
      !~ decomp_npools_vr(c,j,i_dom) = decomp_npools_vr(c,j,i_dom) !- cdons_min(c,j)
      !~ cdons(c,j) = decomp_npools_vr(c,j,i_dom)
      !~ cdons_unsat(c,j) = cdons(c,j) * (1. - finundated(c))
      !~ cdons_sat(c,j) = cdons(c,j) * finundated(c)
      
      cmicbiocs(c,j) = decomp_cpools_vr(c,j,i_bacteria) +  decomp_cpools_vr(c,j,i_fungi)
      cmicbions(c,j) = decomp_npools_vr(c,j,i_bacteria) +  decomp_npools_vr(c,j,i_fungi)
           end do
      end do
!end do
   !~ do j = 1,nlevsoi
      !~ do fc = 1,num_soilc
         !~ c = filter_soilc(fc)
!~ !write(iulog, *) "here2", cdocs(c,j), " doc_k ", doc_k
	!~ cdocs(c,j) = cdocs_unsat(c,j) * (1. - finundated(c)) + cdocs(c,j) * finundated(c)
        !~ dochr_vr(c,j) = cdocs(c,j) * doc_k * 1.5**((t_soisno(c,j) - 273.15 - 25.) / 10) ! * t_scalar(c,j) * w_scalar(c,j)
	!~ cdocs_unsat(c,j) = cdocs_unsat(c,j) * (1. - finundated(c)) * (1. - doc_k * 1.5**((t_soisno(c,j) - 273.15 - 25.) / 10)) ! * t_scalar(c,j) * w_scalar(c,j)
	!~ cdocs_sat(c,j) = cdocs_sat(c,j) * finundated(c) * (1.0 - doc_k * 1.5**((t_soisno(c,j) - 273.15 - 25.) / 10)) ! * t_scalar(c,j) * w_scalar(c,j)
	!~ cdocs_unsat(c,j) = max(0._r8, cdocs_unsat(c,j))
	!~ cdocs_sat(c,j) = max(0._r8, cdocs_sat(c,j))
	!~ cdocs(c,j) = cdocs_unsat(c,j) * (1. - finundated(c)) + cdocs(c,j) * finundated(c)
!~ !write(iulog, *) dochr_vr(c,j), " dochr doc ", cdocs(c,j)
!~ !  calculation of death of the microbes

!~ if(t_soisno(c,j) < Tmmin) then
!~ fm_t = (2. ** ((Tmmin - Tmref)/10)) ** 2
!~ else
!~ fm_t = 2. ** ((t_soisno(c,j) - Tmref)/10)
!~ endif

!~ if((h2osoi_vol(c,j) / watsat(c,j)) < Mmmin) then
!~ fm_m = 0._r8
!~ else
!~ fm_m = 1._r8
!~ endif

!~ micbiohr_vr(c,j) = cmicbios(c,j) * microbeMR / dt * fm_t * fm_m * depth_scalar(c,j) * o_scalar(c,j)
!~ if(cmicbios(c,j) - micbiohr_vr(c,j) > 1e-15) then
!~ cmicbios(c,j) =  cmicbios(c,j) - micbiohr_vr(c,j)
!~ else
!~ micbiohr_vr(c,j) = cmicbios(c,j) - 1e-15
!~ cmicbios(c,j) =  15e-15
!~ endif
      !~ end do
   !~ end do
#endif   
   
#ifdef LCH4
   ! Calculate total fraction of potential HR, for methane code
   do j = 1,nlevdecomp
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         hrsum(c,j) = 0._r8
      end do
   end do
   do k = 1, ndecomp_cascade_transitions
      do j = 1,nlevdecomp
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            hrsum(c,j) = hrsum(c,j) + rf_decomp_cascade(c,j,k) * p_decomp_cpool_loss(c,j,k)
         end do
      end do
   end do
   
   ! Nitrogen limitation / (low)-moisture limitation
   do j = 1,nlevdecomp
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         if (phr_vr(c,j) > 0._r8) then
            fphr(c,j) = hrsum(c,j) / phr_vr(c,j) * w_scalar(c,j)
            fphr(c,j) = max(fphr(c,j), 0.01_r8) ! Prevent overflow errors for 0 respiration
         else
            fphr(c,j) = 1._r8
         end if
      end do
   end do
#endif
   
   ! vertically integrate net and gross mineralization fluxes for diagnostic output
   do j = 1,nlevdecomp
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         net_nmin(c) = net_nmin(c) + net_nmin_vr(c,j) * dzsoi_decomp(j)
         gross_nmin(c) = gross_nmin(c) + gross_nmin_vr(c,j) * dzsoi_decomp(j)   
      end do
   end do
   
 end subroutine CNDecompAlloc
 
 
#endif
 
end module CNDecompMod