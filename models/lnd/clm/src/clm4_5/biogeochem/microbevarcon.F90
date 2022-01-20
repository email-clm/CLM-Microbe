module microbevarcon

#ifdef MICROBE
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: microbevarcon
!
! !DESCRIPTION:
! Module containing microbe parameters and logical switches and routine to read constants from CLM namelist.
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils  , only : endrun
  use clm_varctl  , only : iulog
!
! !PUBLIC TYPES:
  implicit none
  save
!
! Microbe Model Parameters
!
  real(r8), PARAMETER :: micbiocn = 8  
  real(r8), PARAMETER :: micbioMR = 0.123
  real(r8), PARAMETER :: dock = 0.031
  real(r8), PARAMETER :: CUEref = 0.43 
  real(r8), PARAMETER :: CUEt = -0.012  
  real(r8), PARAMETER :: Tcueref = 288.15  
  real(r8), PARAMETER :: Tsref = 283.15  
  real(r8), PARAMETER :: Tmref = 285.15  
  real(r8), PARAMETER :: Tsmin = 272.15 
  real(r8), PARAMETER :: Tmmin = 271.15 
  real(r8), PARAMETER :: Msmin = 0.01  
  real(r8), PARAMETER :: Msmax = 1.0  
  real(r8), PARAMETER :: Mmmin = 0.1
  real(r8), PARAMETER :: MFGbiomin = 1e-15

  integer, parameter :: nummicrobepar = 110
  real(r8) :: q10ch4base = 295._r8  ! Rough estimate from comparison between Walter and previous CLM-CH4 data
  ! Uses Michigan bog data from Shannon & White
  ! This is the temperature at which the effective f_ch4 actually equals the constant f_ch4.

  real(r8) :: q10ch4 = 1.33_r8 ! Production Q10
  ! Note that this is the additional Q10 for methane production ABOVE the soil decomposition temperature relationship.
  ! Corresponds to a methanogenesis Q10 of 2 when SOM HR has Q10 of 1.5.
  ! (No assumption is made in CH4 code about SOM HR temperature relationship.)
  ! Note that this formulation should be improved by making anaerobic decomposition in general a stronger function of
  ! temperature, not just methane production.

  real(r8) :: vmax_ch4_oxid = 45.e-6_r8 * 1000._r8 / 3600._r8 ! [mol/m3-w/s];
  ! 45 uM/h from Walter and Heimann for the Mich. site (2000)
  ! oxidation rate constant (Walter and Heimann 2000)

  real(r8) :: k_m = 5.e-6_r8 * 1000._r8 ! [mol/m3-w]
  ! Michaelis-Menten oxidation rate constant for CH4 concentration (Walter and Heimann 2000)

  real(r8) :: q10_ch4oxid = 1.9_r8 ! Segers, 1998
  ! Q10 oxidation constant (Walter and Heimann 2000)

  real(r8) :: smp_crit = -2.4e5_r8 ! mm. From Schnell & King, 1996.
  ! Critical soil moisture potential to reduce oxidation (mm) due to dessication of methanotrophs above the water table.
  ! To turn off limitation, set to very large negative value.

  real(r8) :: aereoxid = -1._r8 ! fraction of methane flux entering aerenchyma rhizosphere that will be
  ! oxidized rather than emitted.  From Wania.
  ! Note, this has been replaced by prognostic O2 diffusion into aerenchyma and is set to -1 by default.
  ! Set to value between 0 & 1 (inclusive) for sensitivity tests.  In particular, to calculate the 
  ! prognostic fraction analogous to the Wania parameter, set to 0. and compare to a normal run.

  real(r8) :: mino2lim = 0.2_r8 ! minimum anaerobic decomposition rate as a fraction of potential aerobic rate
  ! for soil decomposition (or diagnostic O2-limitation / seasonal inundation factor)

  real(r8) :: rootlitfrac = 0.50_r8 ! Fraction of soil organic matter associated with roots
  ! Used to partition the production between rootfr and the top 5 layers (ifndef VERTSOILC)

  real(r8) :: scale_factor_aere = 1.0_r8 ! scale factor on the aerenchyma area for sensitivity tests

  real(r8) :: vgc_max = 0.15_r8 ! ratio of saturation pressure triggering ebullition

  real(r8) :: organic_max  = 130._r8 ! organic matter content (kg/m3) where soil is assumed to act like peat
  ! for diffusion. Very large values will lead to all soil being treated as mineral. Negative values will lead   ! to all soil being treated as peat.

  real(r8) :: satpow = 2._r8 ! exponent on watsat for saturated soil solute diffusion
  ! (2 = Buckingham / Moldrup; 4/3 = Millington-Quirk)

  real(r8) :: cnscalefactor = 1._r8 ! scale factor on CN decomposition for assigning methane flux
                                    ! This should equal 1 except for sensitivity studies.

  real(r8) :: f_ch4 = 0.2_r8
                          ! originally 25% / (100% + 25%) from Wania. This is the ratio of CH4 production to total C
                          ! mineralization.
                          ! fraction of total decomposition that comes off as CH4 rather than CO2
                          ! Effective value will depend on temperature, redox, & pH but cannot exceed 50%
                          ! based on stoichiometry.
                          ! Note this is a crude parameterization: values in the field vary widely.

  logical :: transpirationloss = .true. ! switch for activating CH4 loss from transpiration
                                      ! Transpiration loss assumes that the methane concentration in dissolved soil
                                      ! water remains constant through the plant and is released when the water evaporates
                                      ! from the stomata.
                                      ! Currently hard-wired to true; impact is < 1 Tg CH4/yr

  real(r8) :: k_m_o2 = 20.e-6_r8 * 1000._r8 ! [mol/m3-w]
  ! Michaelis-Menten oxidation rate constant for O2 concentration (Segers 1998, Lidstrom and Somers 1984)

  real(r8) :: nongrassporosratio = 0.33_r8 ! Ratio of root porosity in non-grass to grass, used for aerenchyma transport
                                           ! Some values in Colmer 2003

 ! logical :: allowlakeprod = .false. ! Switch to allow production under lakes based on soil carbon dataset
                                     ! (Methane can be produced, and CO2 produced from methane oxidation,
                                     ! which will slowly reduce the available carbon stock, if ! replenishlakec, but no other biogeochem is done.)
                                     ! Note: switching this off turns off ALL lake methane biogeochem. However, 0 values
                                     ! will still be averaged into the concentration _sat history fields.

  logical :: usephfact = .false. ! Switch to use pH factor in methane production

  real(r8) :: k_m_unsat = 5.e-6_r8 * 1000._r8 / 10._r8 ! [mol/m3-w]
  ! Michaelis-Menten oxidation rate constant for CH4 concentration: literature suggests that methanotrophs
  ! in upland areas have higher affinity for methane in order to access the low ambient concentrations above the water
  ! table. (See Bender & Conrad, 1992, etc.)

  real(r8) :: vmax_oxid_unsat = 45.e-6_r8 * 1000._r8 / 3600._r8 / 10._r8 ! [mol/m3-w/s]
  ! Literature suggests that while k_m is lower, vmax is also lower in upland areas.

  real(r8) :: scale_factor_gasdiff = 1.0_r8 ! For sensitivity tests; convection would allow this to be > 1

  real(r8) :: scale_factor_liqdiff = 1.0_r8 ! For sensitivity tests; convection would allow this to be > 1

  real(r8) :: redoxlag = 30.0_r8 ! Number of days to lag in the calculation of finundated_lag, which will
                                 ! be used to assess the availability of alternative electron acceptors in recently
                                 ! inundated area, reducing production.
                                 ! Set to 0 to turn off this feature.

  ! New namelists added 6/12/11

  logical :: fin_use_fsat = .false. ! Use fsat rather than the inversion to Prigent satellite inundation obs. (applied to
                                    ! CLM water table depth and surface runoff) to calculated finundated which is
                                    ! used in methane code and potentially soil code
                                    !!!! Attn EK: Set this to true when Sean Swenson's prognostic, tested
                                       ! fsat is integrated. (CLM4 fsat is bad for these purposes.)

 ! real(r8) :: unsat_aere_ratio = 0.05_r8 / 0.3_r8 ! Ratio to multiply upland vegetation aerenchyma porosity by compared to
                                        ! inundated systems. Note: porosity will be kept at above a minimum residual value
                                        ! porosmin set in subroutine ch4_aere.

  logical :: usefrootc = .false.    ! Use CLMCN fine root C rather than ann NPP & LAI based parameterization to
                                    ! calculate tiller C for aerenchyma area calculation.
                                    ! The NPP & LAI param. was based on Wania for Arctic sedges and may not be
                                    ! appropriate for woody PFTs, although nongrassporosratio above partly adjusts
                                    ! for this.  However, using fine root C reduces the aerenchyma area by a large
                                    ! factor.

  logical :: ch4offline = .true.    ! true --> Methane is not passed between the land & atmosphere.
                                    ! NEM is not added to NEE flux to atm. to correct for methane production,
                                    ! and ambient CH4 is set to constant 2009 value.

  logical :: ch4rmcnlim = .false.   ! Remove the N and low moisture limitations on SOM HR when calculating
                                    ! methanogenesis.
                                    ! Note: this option has not been extensively tested.
                                    ! Currently hardwired off.

  logical :: anoxicmicrosites = .false. ! Use Arah & Stephen 1998 expression to allow production above the water table
                                        ! Currently hardwired off; expression is crude.

  logical :: ch4frzout = .false.    ! Exclude CH4 from frozen fraction of soil pore H2O, to simulate "freeze-out" pulse
                                    ! as in Mastepanov 2008.
                                    ! Causes slight increase in emissions in the fall and decrease in the spring.
                                    ! Currently hardwired off; small impact.

  real(r8) :: redoxlag_vertical = 0._r8   ! time lag (days) to inhibit production for newly unsaturated layers
                                          ! when decreasing WT depth for unsat. zone. See the description for redoxlag.

  real(r8) :: atmch4 = 1.79e-6_r8    ! Atmospheric CH4 mixing ratio to prescribe if not provided by the atmospheric model  ! (mol/mol)
  real(r8) :: atmo2 = 209460e-6_r8    ! Atmospheric CH4 mixing ratio to prescribe if not provided by the atmospheric model
  real(r8) :: atmco2 = 397e-6_r8    ! Atmospheric CH4 mixing ratio to prescribe if not provided by the atmospheric model
  real(r8) :: atmh2 = 0.55e-6_r8    ! Atmospheric CH4 mixing ratio to prescribe if not provided by the atmospheric model

  real(r8) :: plant2doc = 0.2             ! fraction of litter decomposition to available carbon
  real(r8) :: som2doc = 0.05          ! fraction of soil organic matter decomposition to available carbon

!xiaofeng xu creared new mechanisms and  added the following new parameters
  	real(r8) :: m_dKAce = 15 ! //3000	//36*1000/12/50/4    /50 for the microbe c in soc, // mMol C m-3 : convert from 36 gC m-3
	real(r8) :: m_dAceProdACmax = 0.005 !;//14.4; //14.4E3; // mMol/m3/d 0.999995; // 1-pow(0.6, 23)  0.4h-1  // Grant R.F, 1998
	real(r8) :: m_dKAceProdO2 = 10 !;  // mMol/m3
	real(r8) :: m_dH2ProdAcemax = 1.35e-5 !;  // mMol/m3   1.35; // 0.05625 * 24mMol / mMol Acetate-C h-1 : convert from 6.75e-3 mMol acetate g-1 h-1
	real(r8) :: m_dKH2ProdAce = 1.65e-3 !;	 // mMol/m3   //16.5e-3;
	real(r8) :: m_dKCO2ProdAce = 8.25e-3 !;		// mMol/m3   8.25e-3;
	real(r8) :: m_dGrowRH2Methanogens = 0.023 !;		// d-1  0.023*24=0.552, (1-0.023)^24 = 0.5589
	real(r8) :: m_dDeadRH2Methanogens = 0.007 !;		//0.0735 d-1   0.016*24  //0.384;  0.007 * 24 = 0.168; 0.002 - 0.012 h-1
	real(r8) :: m_dYH2Methanogens = 0.20 !;
	real(r8) :: m_dGrowRAceMethanogens = 0.02 !;		// d-1
	real(r8) :: m_dDeadRAceMethanogens = 0.001 !;		 //0.01 d-1
	real(r8) :: m_dYAceMethanogens = 0.04 !;   // 0.04; molC/(molAcetate -C) : convert from 0.08 mol C (mol acetate)-1
	real(r8) :: m_dGrowRMethanotrophs = 0.024 !;		// d-1
	real(r8) :: m_dDeadRMethanotrophs = 0.002 !;		//d-1
	real(r8) :: m_dYMethanotrophs = 0.40 !;	// molC/(mol CH4 -C) : convert from 0.4 mol C (mol CH4)-1
	real(r8) :: m_dGrowRAOMMethanotrophs = 0.024 !;		// d-1
	real(r8) :: m_dDeadRAOMMethanotrophs = 0.002 !;		//d-1
	real(r8) :: m_dYAOMMethanotrophs = 0.40 !;	// molC/(mol CH4 -C) : convert from 0.4 mol C (mol CH4)-1
	real(r8) :: m_dAceProdQ10 = 2.0 !;
	real(r8) :: m_dACProdQ10 = 1.5 !;
	real(r8) :: m_dACMinQ10 = 1.5 !;
	real(r8) :: m_dAceH2min = 0.48e-3 !;		 //0.48e-3;   mMol/m3
	real(r8) :: m_dCH4H2min = 0.27e-4 !;	 // 0.27e-4;   mMol/m3
	real(r8) :: m_dKH2ProdCH4 = 7.75e-3 !;	// mMol/m3
	real(r8) :: m_dKCO2ProdCH4 = 1.98e-3 !;	 // mMol/m3
	real(r8) :: m_dH2CH4ProdQ10 = 3.5 !;
	real(r8) :: m_dH2AceProdQ10 = 3.5 !;
	real(r8) :: m_dKCH4ProdAce = 5 !; //0.005 !;		 // mMol/m3
	real(r8) :: m_dKCH4ProdO2 = 10 !; // 0.01;	// mMol/m3
	real(r8) :: m_dCH4ProdQ10 = 4 !;
	real(r8) :: m_drCH4Prod = 0.5 !;		// MolCH4(Mol acetate -C)-1 : convert from 2 MolCH4(Mol acetate)-1
	real(r8) :: m_dKCH4OxidCH4 = 5 !; //0.005;		// mMol/m3
	real(r8) :: m_dKAOMCH4OxidCH4 = 2
	real(r8) :: m_dKCH4OxidO2 = 10 !; // 0.01;		// mMol/m3
	real(r8) :: m_dCH4OxidQ10 = 2.0 !;
	real(r8) :: m_dAOMCH4OxidQ10 = 2.0 !
	real(r8) :: m_drAer = 2 !;					// MolO2(mol C)-1
	real(r8) :: m_dKAerO2 = 10 !;				 //0.01 mMol/m3
	real(r8) :: m_dAerDecomQ10 = 2.0 !;
	real(r8) :: m_drCH4Oxid = 2 !;			// MolO2(molCH4)-1
	real(r8) :: m_dKe = 0.03 !;			 // 0.05 h-1
	real(r8) :: m_dCH4min = 0.5 !;			 //0.5 mMol/m3
	real(r8) :: m_dAirCH4 = 0.0893 !;  // mM/m3
	real(r8) :: m_dAirH2 = 0.0257 !; // mM/m3
	real(r8) :: m_dAirO2 = 9.375e3 !; // mM/m3
	real(r8) :: m_dAirCO2 = 16.295 !; // mM/m3
	
	! not in the parameter file
	real(r8) :: m_dPlantTrans = 0.68 / 48.0
	real(r8) :: g_dMaxH2inWater = 4.73e-4
	
	!Penning, H., P. Claus, P. Casper, and R. Conrad. 2006. Carbon isotope fractionation during acetoclastic methanogenesis by Methanosaeta concilii in culture and a lake sediment. 
	! Applied and Environmental Microbiology 72:5648-5652.
	real(r8) :: frac_doc = 1.0			! fractionation factor for DOC production
	real(r8) :: frac_ace = 1.0			! fractionation factor for acetate production
	real(r8) :: frac_acch4 = 1.0		! fractionation factor for aceclatic methanogenesis
	real(r8) :: frac_hych4 = 1.0		! fractionation factor for hydrogenotrophic methanogenesis
 	real(r8) :: frac_acetogenesis = 1.0	! fractionation factor for acetogenesis
	real(r8) :: frac_ch4ox = 1.0		! fractionation factor for methanotrophy
 	real(r8) :: frac_ch4aom = 1.0		! fractionation factor for anaerobic oxidation of methane
  
	real(r8) :: pHmin = 4.0
	real(r8) :: pHmax = 10.
	real(r8) :: pHopt = 7.
  
	real(r8) :: k_dom = 0.042
	real(r8) :: k_bacteria = 0.56
	real(r8) :: k_fungi = 0.56
	real(r8) :: dom_diffus = 10. / 3600. / 365.
	real(r8) :: m_Fick_ad = 0.75
	
	real(r8) :: m_rf_s1m = 0.28
	real(r8) :: m_rf_s2m = 0.46
	real(r8) :: m_rf_s3m = 0.55
	real(r8) :: m_rf_s4m = 0.75
  
	real(r8) :: m_batm_f = 0.05
	real(r8) :: m_bdom_f = 0.25
	real(r8) :: m_bs1_f = 0.1
	real(r8) :: m_bs2_f = 0.12
	real(r8) :: m_bs3_f = 0.18
	real(r8) :: m_fatm_f = 0.05
	real(r8) :: m_fdom_f = 0.25
	real(r8) :: m_fs1_f = 0.1
	real(r8) :: m_fs2_f = 0.12
	real(r8) :: m_fs3_f = 0.18
	real(r8) :: m_domb_f = 0.3
	real(r8) :: m_domf_f = 0.3
	real(r8) :: m_doms1_f = 0.2
	real(r8) :: m_doms2_f = 0.15
	real(r8) :: m_doms3_f = 0.05
  
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: ch4conrd ! Read and initialize CH4 constants
!
! !REVISION HISTORY:
! Created by Zack Subin
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ch4conrd
!
! !INTERFACE:
  subroutine ch4conrd ()
!
! !DESCRIPTION:
! Read and initialize CH4 constants
!
! !USES:
    use fileutils , only : opnfil, getfil, relavu, getavu
    use spmdMod   , only : masterproc, mpicom, MPI_REAL8, MPI_LOGICAL
    !use controlMod  , only: NLFilename

!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
! routine initch4
!
! !REVISION HISTORY:
! Created by Zack Subin
! Modified 6/12/2011 to use namelist rather than ASCII file.
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: i,j,n                ! loop indices
    integer :: ierr                 ! error code
    integer :: unitn                ! unit for namelist file
    character(len=40) :: ch4parname(nummicrobepar)  ! subroutine name
    character(LEN=256)::locfn='./microbepar_in'
!-----------------------------------------------------------------------

real(r8)::dummy(nummicrobepar)
    if (masterproc) then

       write(iulog,*) 'Attempting to read CH4 parameters .....'
       unitn = getavu()
       call opnfil(locfn, unitn, 'f')
       do i = 1, nummicrobepar
       read(unitn, *, IOSTAT=ierr) ch4parname(i), dummy(i)
       if(ierr/=0) then 
       write(iulog, *) 'error in reading in microbepar_in'
       call endrun()
       end if
       end do
       call relavu( unitn )

i=1

q10ch4base                         = dummy(i); i=i+1
    q10ch4 = dummy(i); i=i+1
    vmax_ch4_oxid = dummy(i); i=i+1
    k_m = dummy(i); i=i+1
    q10_ch4oxid = dummy(i); i=i+1
    smp_crit = dummy(i); i=i+1
    aereoxid = dummy(i); i=i+1
    mino2lim = dummy(i); i=i+1
    rootlitfrac = dummy(i); i=i+1
    scale_factor_aere = dummy(i); i=i+1
    vgc_max = dummy(i); i=i+1
    organic_max = dummy(i); i=i+1
    satpow = dummy(i); i=i+1
    cnscalefactor = dummy(i); i=i+1
    f_ch4 = dummy(i); i=i+1
    k_m_o2 = dummy(i); i=i+1
    nongrassporosratio = dummy(i); i=i+1
 !   allowlakeprod = dummy(i); i=i+1
 !   lake_decomp_fact = dummy(i); i=i+1
    if (dummy(i) .eq. 1) usephfact=.true.
    i=i+1
!    usephfact = int(dummy(i)); i=i+1
    k_m_unsat = dummy(i); i=i+1
    vmax_oxid_unsat = dummy(i); i=i+1
 !   replenishlakec = dummy(i); i=i+1
    scale_factor_gasdiff = dummy(i); i=i+1
    scale_factor_liqdiff = dummy(i); i=i+1
    redoxlag = dummy(i); i=i+1
    if (dummy(i) .eq. 1) usefrootc=.true.
    i=i+1
!    unsat_aere_ratio = dummy(i); i=i+1
    if (dummy(i) .eq. 1) ch4offline=.true.
    i=i+1
!    ch4offline = dummy(i); i=i+1
    redoxlag_vertical = dummy(i); i=i+1
    if (dummy(i) .eq. 1) fin_use_fsat=.true.
    i=i+1
    atmch4 = dummy(i); i=i+1
    plant2doc = dummy(i); i=i+1
    som2doc = dummy(i); i=i+1
    
	m_dKAce = dummy(i); i=i+1
	m_dAceProdACmax = dummy(i); i=i+1
	m_dKAceProdO2 = dummy(i); i=i+1
	m_dH2ProdAcemax = dummy(i); i=i+1
	m_dKH2ProdAce = dummy(i); i=i+1
	m_dKCO2ProdAce = dummy(i); i=i+1
	m_dGrowRH2Methanogens = dummy(i); i=i+1
	m_dDeadRH2Methanogens = dummy(i); i=i+1
	m_dYH2Methanogens = dummy(i); i=i+1
	m_dGrowRAceMethanogens = dummy(i); i=i+1
	m_dDeadRAceMethanogens = dummy(i); i=i+1
	m_dYAceMethanogens = dummy(i); i=i+1
	m_dGrowRMethanotrophs = dummy(i); i=i+1
	m_dDeadRMethanotrophs = dummy(i); i=i+1
	m_dYMethanotrophs = dummy(i); i=i+1
	m_dGrowRAOMMethanotrophs = dummy(i); i=i+1
	m_dDeadRAOMMethanotrophs = dummy(i); i=i+1
	m_dYAOMMethanotrophs = dummy(i); i=i+1
	m_dAceProdQ10 = dummy(i); i=i+1
	m_dACProdQ10 = dummy(i); i=i+1
	m_dACMinQ10 = dummy(i); i=i+1
	m_dAceH2min = dummy(i); i=i+1
	m_dCH4H2min = dummy(i); i=i+1
	m_dKH2ProdCH4 = dummy(i); i=i+1
	m_dKCO2ProdCH4 = dummy(i); i=i+1
	m_dH2CH4ProdQ10 = dummy(i); i=i+1
	m_dH2AceProdQ10 = dummy(i); i=i+1
	m_dKCH4ProdAce = dummy(i); i=i+1
	m_dKCH4ProdO2 = dummy(i); i=i+1
	m_dCH4ProdQ10 = dummy(i); i=i+1
	m_drCH4Prod = dummy(i); i=i+1
	m_dKCH4OxidCH4 = dummy(i); i=i+1
	m_dKAOMCH4OxidCH4 = dummy(i); i=i+1
	m_dKCH4OxidO2 = dummy(i); i=i+1
	m_dCH4OxidQ10 = dummy(i); i=i+1
	m_dAOMCH4OxidQ10 = dummy(i); i=i+1
	m_drAer = dummy(i); i=i+1
	m_dKAerO2 = dummy(i); i=i+1
	m_dAerDecomQ10 = dummy(i); i=i+1
	m_drCH4Oxid = dummy(i); i=i+1
	m_dKe = dummy(i); i=i+1
	m_dCH4min = dummy(i); i=i+1
	m_dAirCH4 = dummy(i); i=i+1
	m_dAirH2 = dummy(i); i=i+1
	m_dAirO2 = dummy(i); i=i+1
	m_dAirCO2 = dummy(i); i=i+1
	
	frac_doc = dummy(i); i=i+1	
	frac_ace = dummy(i); i=i+1	
	frac_acch4 = dummy(i); i=i+1	
	frac_hych4 = dummy(i); i=i+1		
	frac_acetogenesis = dummy(i); i=i+1		
	frac_ch4ox = dummy(i); i=i+1		
	frac_ch4aom = dummy(i); i=i+1		
	
	k_dom = dummy(i); i=i+1
	k_bacteria = dummy(i); i=i+1
	k_fungi = dummy(i); i=i+1
	dom_diffus = dummy(i); i=i+1
	m_Fick_ad = dummy(i); i=i+1
	m_rf_s1m = dummy(i); i=i+1
	m_rf_s2m = dummy(i); i=i+1
	m_rf_s3m = dummy(i); i=i+1
	m_rf_s4m = dummy(i); i=i+1  
	m_batm_f = dummy(i); i=i+1
	m_bdom_f = dummy(i); i=i+1
	m_bs1_f = dummy(i); i=i+1
	m_bs2_f = dummy(i); i=i+1
	m_bs3_f = dummy(i); i=i+1
	m_fatm_f = dummy(i); i=i+1
	m_fdom_f = dummy(i); i=i+1
	m_fs1_f = dummy(i); i=i+1
	m_fs2_f = dummy(i); i=i+1
	m_fs3_f = dummy(i); i=i+1
	m_domb_f = dummy(i); i=i+1
	m_domf_f = dummy(i); i=i+1
	m_doms1_f = dummy(i); i=i+1
	m_doms2_f = dummy(i); i=i+1
	m_doms3_f = dummy(i); i=i+1
	m_dPlantTrans = dummy(i); i=i+1
    
    !xiaofeng xu creared new mechanisms and the new parameters         
end if

    call mpi_bcast (q10ch4base, 1 , MPI_REAL8, 0, mpicom, ierr)         
    call mpi_bcast (q10ch4, 1 , MPI_REAL8, 0, mpicom, ierr)             
    call mpi_bcast (vmax_ch4_oxid, 1 , MPI_REAL8, 0, mpicom, ierr)     
    call mpi_bcast (k_m, 1 , MPI_REAL8, 0, mpicom, ierr)               
    call mpi_bcast (q10_ch4oxid, 1 , MPI_REAL8, 0, mpicom, ierr)       
    call mpi_bcast (smp_crit, 1 , MPI_REAL8, 0, mpicom, ierr)          
    call mpi_bcast (aereoxid, 1 , MPI_REAL8, 0, mpicom, ierr)          
    call mpi_bcast (mino2lim, 1 , MPI_REAL8, 0, mpicom, ierr)          
    call mpi_bcast (rootlitfrac, 1 , MPI_REAL8, 0, mpicom, ierr)       
    call mpi_bcast (scale_factor_aere, 1 , MPI_REAL8, 0, mpicom, ierr) 
    call mpi_bcast (vgc_max, 1 , MPI_REAL8, 0, mpicom, ierr)           
    call mpi_bcast (organic_max, 1 , MPI_REAL8, 0, mpicom, ierr)       
    call mpi_bcast (satpow, 1 , MPI_REAL8, 0, mpicom, ierr)            
    call mpi_bcast (cnscalefactor, 1 , MPI_REAL8, 0, mpicom, ierr)            
    call mpi_bcast (f_ch4, 1 , MPI_REAL8, 0, mpicom, ierr)            
    !call mpi_bcast (transpirationloss, 1 , MPI_LOGICAL, 0, mpicom, ierr)            
    call mpi_bcast (k_m_o2, 1 , MPI_REAL8, 0, mpicom, ierr)            
    call mpi_bcast (nongrassporosratio, 1 , MPI_REAL8, 0, mpicom, ierr)            
!    call mpi_bcast (allowlakeprod, 1 , MPI_LOGICAL, 0, mpicom, ierr)            
!    call mpi_bcast (lake_decomp_fact, 1 , MPI_REAL8, 0, mpicom, ierr)            
    call mpi_bcast (usephfact, 1 , MPI_LOGICAL, 0, mpicom, ierr)            
    call mpi_bcast (k_m_unsat, 1 , MPI_REAL8, 0, mpicom, ierr)            
    call mpi_bcast (vmax_oxid_unsat, 1 , MPI_REAL8, 0, mpicom, ierr)            
 !   call mpi_bcast (replenishlakec, 1 , MPI_LOGICAL, 0, mpicom, ierr)            
    call mpi_bcast (scale_factor_gasdiff, 1 , MPI_REAL8, 0, mpicom, ierr)            
    call mpi_bcast (scale_factor_liqdiff, 1 , MPI_REAL8, 0, mpicom, ierr)            
    call mpi_bcast (redoxlag, 1 , MPI_REAL8, 0, mpicom, ierr)            
    call mpi_bcast (fin_use_fsat, 1 , MPI_LOGICAL, 0, mpicom, ierr)            
!    call mpi_bcast (unsat_aere_ratio, 1 , MPI_REAL8, 0, mpicom, ierr)            
    call mpi_bcast (usefrootc, 1 , MPI_LOGICAL, 0, mpicom, ierr)            
    call mpi_bcast (ch4offline, 1 , MPI_LOGICAL, 0, mpicom, ierr)            
    !call mpi_bcast (ch4rmcnlim, 1 , MPI_LOGICAL, 0, mpicom, ierr)            
    !call mpi_bcast (anoxicmicrosites, 1 , MPI_LOGICAL, 0, mpicom, ierr)            
    !call mpi_bcast (ch4frzout, 1 , MPI_LOGICAL, 0, mpicom, ierr)            
    call mpi_bcast (redoxlag_vertical, 1 , MPI_REAL8, 0, mpicom, ierr)            
    call mpi_bcast (atmch4, 1 , MPI_REAL8, 0, mpicom, ierr)
    call mpi_bcast (plant2doc, 1 , MPI_REAL8, 0, mpicom, ierr)
    call mpi_bcast (som2doc, 1 , MPI_REAL8, 0, mpicom, ierr)
    
    call mpi_bcast (m_dKAce, 1 , MPI_REAL8, 0, mpicom, ierr)
    call mpi_bcast (m_dAceProdACmax, 1 , MPI_REAL8, 0, mpicom, ierr)
    call mpi_bcast (m_dKAceProdO2, 1 , MPI_REAL8, 0, mpicom, ierr)
    call mpi_bcast (m_dH2ProdAcemax, 1 , MPI_REAL8, 0, mpicom, ierr)
    call mpi_bcast (m_dKH2ProdAce, 1 , MPI_REAL8, 0, mpicom, ierr)
    call mpi_bcast (m_dKCO2ProdAce, 1 , MPI_REAL8, 0, mpicom, ierr)
    call mpi_bcast (m_dGrowRH2Methanogens, 1 , MPI_REAL8, 0, mpicom, ierr)
    call mpi_bcast (m_dDeadRH2Methanogens, 1 , MPI_REAL8, 0, mpicom, ierr)
    call mpi_bcast (m_dYH2Methanogens, 1 , MPI_REAL8, 0, mpicom, ierr)
    call mpi_bcast (m_dGrowRAceMethanogens, 1 , MPI_REAL8, 0, mpicom, ierr)
    call mpi_bcast (m_dDeadRAceMethanogens, 1 , MPI_REAL8, 0, mpicom, ierr)    
    call mpi_bcast (m_dYAceMethanogens, 1 , MPI_REAL8, 0, mpicom, ierr)
    call mpi_bcast (m_dGrowRMethanotrophs, 1 , MPI_REAL8, 0, mpicom, ierr)
    call mpi_bcast (m_dDeadRMethanotrophs, 1 , MPI_REAL8, 0, mpicom, ierr)
    call mpi_bcast (m_dYMethanotrophs, 1 , MPI_REAL8, 0, mpicom, ierr)
    call mpi_bcast (m_dGrowRAOMMethanotrophs, 1 , MPI_REAL8, 0, mpicom, ierr)    
    call mpi_bcast (m_dDeadRAOMMethanotrophs, 1 , MPI_REAL8, 0, mpicom, ierr)
    call mpi_bcast (m_dYAOMMethanotrophs, 1 , MPI_REAL8, 0, mpicom, ierr)
    call mpi_bcast (m_dAceProdQ10, 1 , MPI_REAL8, 0, mpicom, ierr)
    call mpi_bcast (m_dACProdQ10, 1 , MPI_REAL8, 0, mpicom, ierr)
    call mpi_bcast (m_dACMinQ10, 1 , MPI_REAL8, 0, mpicom, ierr)    
    call mpi_bcast (m_dAceH2min, 1 , MPI_REAL8, 0, mpicom, ierr)
    call mpi_bcast (m_dCH4H2min, 1 , MPI_REAL8, 0, mpicom, ierr)
    call mpi_bcast (m_dKH2ProdCH4, 1 , MPI_REAL8, 0, mpicom, ierr)
    call mpi_bcast (m_dKCO2ProdCH4, 1 , MPI_REAL8, 0, mpicom, ierr)
    call mpi_bcast (m_dH2CH4ProdQ10, 1 , MPI_REAL8, 0, mpicom, ierr)    
    call mpi_bcast (m_dH2AceProdQ10, 1 , MPI_REAL8, 0, mpicom, ierr)
    call mpi_bcast (m_dKCH4ProdAce, 1 , MPI_REAL8, 0, mpicom, ierr)
    call mpi_bcast (m_dKCH4ProdO2, 1 , MPI_REAL8, 0, mpicom, ierr)
    call mpi_bcast (m_dCH4ProdQ10, 1 , MPI_REAL8, 0, mpicom, ierr)
    call mpi_bcast (m_drCH4Prod, 1 , MPI_REAL8, 0, mpicom, ierr)    
    call mpi_bcast (m_dKCH4OxidCH4, 1 , MPI_REAL8, 0, mpicom, ierr)
    call mpi_bcast (m_dKAOMCH4OxidCH4, 1 , MPI_REAL8, 0, mpicom, ierr)
    call mpi_bcast (m_dKCH4OxidO2, 1 , MPI_REAL8, 0, mpicom, ierr)
    call mpi_bcast (m_dCH4OxidQ10, 1 , MPI_REAL8, 0, mpicom, ierr)
    call mpi_bcast (m_dAOMCH4OxidQ10, 1 , MPI_REAL8, 0, mpicom, ierr)    
    call mpi_bcast (m_drAer, 1 , MPI_REAL8, 0, mpicom, ierr)
    call mpi_bcast (m_dKAerO2, 1 , MPI_REAL8, 0, mpicom, ierr)
    call mpi_bcast (m_dAerDecomQ10, 1 , MPI_REAL8, 0, mpicom, ierr)
    call mpi_bcast (m_drCH4Oxid, 1 , MPI_REAL8, 0, mpicom, ierr)
    call mpi_bcast (m_dKe, 1 , MPI_REAL8, 0, mpicom, ierr)    
    call mpi_bcast (m_dCH4min, 1 , MPI_REAL8, 0, mpicom, ierr)
    call mpi_bcast (m_dAirCH4, 1 , MPI_REAL8, 0, mpicom, ierr)
    call mpi_bcast (m_dAirH2, 1 , MPI_REAL8, 0, mpicom, ierr)
    call mpi_bcast (m_dAirO2, 1 , MPI_REAL8, 0, mpicom, ierr)
    call mpi_bcast (m_dAirCO2, 1 , MPI_REAL8, 0, mpicom, ierr)
 
    call mpi_bcast (frac_doc, 1 , MPI_REAL8, 0, mpicom, ierr)
    call mpi_bcast (frac_ace, 1 , MPI_REAL8, 0, mpicom, ierr)
    call mpi_bcast (frac_acch4, 1 , MPI_REAL8, 0, mpicom, ierr)
    call mpi_bcast (frac_hych4, 1 , MPI_REAL8, 0, mpicom, ierr)
    call mpi_bcast (frac_acetogenesis, 1 , MPI_REAL8, 0, mpicom, ierr)
    call mpi_bcast (frac_ch4ox, 1 , MPI_REAL8, 0, mpicom, ierr)
    call mpi_bcast (frac_ch4aom, 1 , MPI_REAL8, 0, mpicom, ierr)
    	
    call mpi_bcast (k_dom, 1 , MPI_REAL8, 0, mpicom, ierr)
    call mpi_bcast (k_bacteria, 1 , MPI_REAL8, 0, mpicom, ierr)
    call mpi_bcast (k_fungi, 1 , MPI_REAL8, 0, mpicom, ierr)
    call mpi_bcast (dom_diffus, 1 , MPI_REAL8, 0, mpicom, ierr)
    call mpi_bcast (m_Fick_ad, 1 , MPI_REAL8, 0, mpicom, ierr)
    
    call mpi_bcast (m_rf_s1m, 1 , MPI_REAL8, 0, mpicom, ierr)
    call mpi_bcast (m_rf_s2m, 1 , MPI_REAL8, 0, mpicom, ierr)
    call mpi_bcast (m_rf_s3m, 1 , MPI_REAL8, 0, mpicom, ierr)
    call mpi_bcast (m_rf_s4m, 1 , MPI_REAL8, 0, mpicom, ierr)
    call mpi_bcast (m_batm_f, 1 , MPI_REAL8, 0, mpicom, ierr)
    call mpi_bcast (m_bdom_f, 1 , MPI_REAL8, 0, mpicom, ierr)
    call mpi_bcast (m_bs1_f, 1 , MPI_REAL8, 0, mpicom, ierr)
    call mpi_bcast (m_bs2_f, 1 , MPI_REAL8, 0, mpicom, ierr)
    call mpi_bcast (m_bs3_f, 1 , MPI_REAL8, 0, mpicom, ierr)
    call mpi_bcast (m_fatm_f, 1 , MPI_REAL8, 0, mpicom, ierr)
    call mpi_bcast (m_fdom_f, 1 , MPI_REAL8, 0, mpicom, ierr)
    call mpi_bcast (m_fs1_f, 1 , MPI_REAL8, 0, mpicom, ierr)
    call mpi_bcast (m_fs2_f, 1 , MPI_REAL8, 0, mpicom, ierr)
    call mpi_bcast (m_fs3_f, 1 , MPI_REAL8, 0, mpicom, ierr)   
    call mpi_bcast (m_domb_f, 1 , MPI_REAL8, 0, mpicom, ierr)
    call mpi_bcast (m_domf_f, 1 , MPI_REAL8, 0, mpicom, ierr)
    call mpi_bcast (m_doms1_f, 1 , MPI_REAL8, 0, mpicom, ierr)
    call mpi_bcast (m_doms2_f, 1 , MPI_REAL8, 0, mpicom, ierr)
    call mpi_bcast (m_doms3_f, 1 , MPI_REAL8, 0, mpicom, ierr)
    call mpi_bcast (m_dPlantTrans, 1 , MPI_REAL8, 0, mpicom, ierr)
    
    if (masterproc) then
       write(iulog,*) 'Successfully read CH4 namelist'
       write(iulog,*)' '
       write(iulog,*)'q10ch4base = ', q10ch4base
       write(iulog,*)'q10ch4 = ', q10ch4
       write(iulog,*)'vmax_ch4_oxid = ', vmax_ch4_oxid
       write(iulog,*)'k_m = ', k_m
       write(iulog,*)'q10_ch4oxid = ', q10_ch4oxid
       write(iulog,*)'smp_crit = ', smp_crit
       write(iulog,*)'aereoxid = ', aereoxid
       write(iulog,*)'mino2lim = ', mino2lim
       write(iulog,*)'rootlitfrac = ', rootlitfrac
       write(iulog,*)'scale_factor_aere = ', scale_factor_aere
       write(iulog,*)'vgc_max = ', vgc_max
       write(iulog,*)'organic_max = ', organic_max
       write(iulog,*)'satpow = ', satpow
       write(iulog,*)'cnscalefactor = ', cnscalefactor
       write(iulog,*)'f_ch4 = ', f_ch4
       !write(iulog,*)'transpirationloss = ', transpirationloss
       write(iulog,*)'k_m_o2 = ', k_m_o2
       write(iulog,*)'nongrassporosratio = ', nongrassporosratio
!       write(iulog,*)'allowlakeprod = ', allowlakeprod
!       write(iulog,*)'lake_decomp_fact = ', lake_decomp_fact
       write(iulog,*)'usephfact = ', usephfact
       write(iulog,*)'k_m_unsat = ', k_m_unsat
       write(iulog,*)'vmax_oxid_unsat = ', vmax_oxid_unsat
!       write(iulog,*)'replenishlakec = ', replenishlakec
       write(iulog,*)'scale_factor_gasdiff = ', scale_factor_gasdiff
       write(iulog,*)'scale_factor_liqdiff = ', scale_factor_liqdiff
       write(iulog,*)'redoxlag = ', redoxlag
       write(iulog,*)'fin_use_fsat = ', fin_use_fsat
 !      write(iulog,*)'unsat_aere_ratio = ', unsat_aere_ratio
       write(iulog,*)'usefrootc = ', usefrootc
       write(iulog,*)'ch4offline = ', ch4offline
       !write(iulog,*)'ch4rmcnlim = ', ch4rmcnlim
       !write(iulog,*)'anoxicmicrosites = ', anoxicmicrosites
       !write(iulog,*)'ch4frzout = ', ch4frzout
       write(iulog,*)'redoxlag_vertical = ', redoxlag_vertical
       write(iulog,*)'atmch4 = ', atmch4
       write(iulog,*)'plant2doc = ', plant2doc
       write(iulog,*)'som2doc = ', som2doc

	write(iulog,*)'m_dKAce = ', m_dKAce
	write(iulog,*)'m_dAceProdACmax = ', m_dAceProdACmax
	write(iulog,*)'m_dKAceProdO2 = ', m_dKAceProdO2
	write(iulog,*)'m_dH2ProdAcemax = ', m_dH2ProdAcemax
	write(iulog,*)'m_dKH2ProdAce = ', m_dKH2ProdAce
	write(iulog,*)'m_dKCO2ProdAce = ', m_dKCO2ProdAce
	write(iulog,*)'m_dGrowRH2Methanogens = ', m_dGrowRH2Methanogens
	write(iulog,*)'m_dDeadRH2Methanogens = ', m_dDeadRH2Methanogens
	write(iulog,*)'m_dYH2Methanogens = ', m_dYH2Methanogens
	write(iulog,*)'m_dGrowRAceMethanogens = ', m_dGrowRAceMethanogens
	write(iulog,*)'m_dDeadRAceMethanogens = ', m_dDeadRAceMethanogens
	write(iulog,*)'m_dYAceMethanogens = ', m_dYAceMethanogens
	write(iulog,*)'m_dGrowRMethanotrophs = ', m_dGrowRMethanotrophs
	write(iulog,*)'m_dDeadRMethanotrophs = ', m_dDeadRMethanotrophs
	write(iulog,*)'m_dYMethanotrophs = ', m_dYMethanotrophs
	write(iulog,*)'m_dGrowRAOMMethanotrophs = ', m_dGrowRAOMMethanotrophs
	write(iulog,*)'m_dDeadRAOMMethanotrophs = ', m_dDeadRAOMMethanotrophs
	write(iulog,*)'m_dYAOMMethanotrophs = ', m_dYAOMMethanotrophs
	write(iulog,*)'m_dAceProdQ10 = ', m_dAceProdQ10
	write(iulog,*)'m_dACProdQ10 = ', m_dACProdQ10
	write(iulog,*)'m_dACMinQ10 = ', m_dACMinQ10
	write(iulog,*)'m_dAceH2min = ', m_dAceH2min
	write(iulog,*)'m_dCH4H2min = ', m_dCH4H2min
	write(iulog,*)'m_dKH2ProdCH4 = ', m_dKH2ProdCH4
	write(iulog,*)'m_dKCO2ProdCH4 = ', m_dKCO2ProdCH4
	write(iulog,*)'m_dH2CH4ProdQ10 = ', m_dH2CH4ProdQ10
	write(iulog,*)'m_dH2AceProdQ10 = ', m_dH2AceProdQ10
	write(iulog,*)'m_dKCH4ProdAce = ', m_dKCH4ProdAce
	write(iulog,*)'m_dKCH4ProdO2 = ', m_dKCH4ProdO2
	write(iulog,*)'m_dCH4ProdQ10 = ', m_dCH4ProdQ10
	write(iulog,*)'m_drCH4Prod = ', m_drCH4Prod
	write(iulog,*)'m_dKCH4OxidCH4 = ', m_dKCH4OxidCH4
	write(iulog,*)'m_dKAOMCH4OxidCH4 = ', m_dKAOMCH4OxidCH4
	write(iulog,*)'m_dKCH4OxidO2 = ', m_dKCH4OxidO2
	write(iulog,*)'m_dCH4OxidQ10 = ', m_dCH4OxidQ10
	write(iulog,*)'m_dAOMCH4OxidQ10 = ', m_dAOMCH4OxidQ10
	write(iulog,*)'m_drAer = ', m_drAer
	write(iulog,*)'m_dKAerO2 = ', m_dKAerO2
	write(iulog,*)'m_dAerDecomQ10 = ', m_dAerDecomQ10
	write(iulog,*)'m_drCH4Oxid = ', m_drCH4Oxid
	write(iulog,*)'m_dKe = ', m_dKe
	write(iulog,*)'m_dCH4min = ', m_dCH4min
	write(iulog,*)'m_dAirCH4 = ', m_dAirCH4
	write(iulog,*)'m_dAirH2 = ', m_dAirH2
	write(iulog,*)'m_dAirO2 = ', m_dAirO2
	write(iulog,*)'m_dAirCO2 = ', m_dAirCO2
	
	write(iulog,*)'frac_doc = ', frac_doc
	write(iulog,*)'frac_ace = ', frac_ace
	write(iulog,*)'frac_acch4 = ', frac_acch4
	write(iulog,*)'frac_hych4 = ', frac_hych4
	write(iulog,*)'frac_acetogenesis = ', frac_acetogenesis
	write(iulog,*)'frac_ch4ox = ', frac_ch4ox
	write(iulog,*)'frac_ch4aom = ', frac_ch4aom
	
	write(iulog,*)'k_dom = ', k_dom
	write(iulog,*)'k_bacteria = ', k_bacteria
	write(iulog,*)'k_fungi = ', k_fungi
	write(iulog,*)'dom_diffus = ', dom_diffus
	
	write(iulog,*)'m_rf_s1m = ', m_rf_s1m
	write(iulog,*)'m_rf_s2m = ', m_rf_s2m
	write(iulog,*)'m_rf_s3m = ', m_rf_s3m
	write(iulog,*)'m_rf_s4m = ', m_rf_s4m
	write(iulog,*)'m_batm_f = ', m_batm_f
	write(iulog,*)'m_bdom_f = ', m_bdom_f
	write(iulog,*)'m_bs1_f = ', m_bs1_f
	write(iulog,*)'m_bs2_f = ', m_bs2_f
	write(iulog,*)'m_bs3_f = ', m_bs3_f
	write(iulog,*)'m_fatm_f = ', m_fatm_f
	write(iulog,*)'m_fdom_f = ', m_fdom_f
	write(iulog,*)'m_fs1_f = ', m_fs1_f
	write(iulog,*)'m_fs2_f = ', m_fs2_f
	write(iulog,*)'m_fs3_f = ', m_fs3_f
	write(iulog,*)'m_domb_f = ', m_domb_f
	write(iulog,*)'m_domf_f = ', m_domf_f
	write(iulog,*)'m_doms1_f = ', m_doms1_f
	write(iulog,*)'m_doms2_f = ', m_doms2_f
	write(iulog,*)'m_doms3_f = ', m_doms3_f
	write(iulog,*)'m_dPlantTrans = ', m_dPlantTrans
	
       if (ch4offline) write(iulog,*)'CH4 Model will be running offline and not affect fluxes to atmosphere.'
       if (aereoxid >= 0._r8) write(iulog,*) 'Fixed aerenchyma oxidation has been selected.'
!       if (.not. allowlakeprod) write(iulog,*) 'Lake production has been disabled.  Lakes will not factor into CH4 BGC.  "Sat" history fields will not average over lakes except for concentrations, which will average zero from lakes.'
       write(iulog,*)'Successfully initialized CH4 parameters'
       write(iulog,*)

    end if

  end subroutine ch4conrd

#endif
! defined MICROBE
end module microbevarcon

