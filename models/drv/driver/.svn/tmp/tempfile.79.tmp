module mrg_x2o_mct

  use shr_kind_mod, only: r8 => shr_kind_r8, cl => shr_kind_cl
  use mct_mod
  use seq_flds_mod
  use seq_flds_indices
  use seq_comm_mct
  use seq_cdata_mod
  use seq_infodata_mod

  implicit none
  save
  private ! except

!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------

  public :: mrg_x2o_init_mct
  public :: mrg_x2o_run_mct
  public :: mrg_x2o_final_mct

!--------------------------------------------------------------------------
! Private data
!--------------------------------------------------------------------------

!===========================================================================================
contains
!===========================================================================================

  subroutine mrg_x2o_init_mct( cdata_o, a2x_o, i2x_o, r2x_o )

    !----------------------------------------------------------------------- 
    type(seq_cdata), intent(in)    :: cdata_o
    type(mct_aVect), intent(inout) :: a2x_o
    type(mct_aVect), intent(inout) :: i2x_o
    type(mct_aVect), intent(inout) :: r2x_o

    type(mct_gsMap), pointer       :: gsMap_ocn
    integer                        :: mpicom
    !----------------------------------------------------------------------- 

    ! Set gsMap
    call seq_cdata_setptrs(cdata_o, gsMap=gsMap_ocn, mpicom=mpicom)

    ! Initialize av for atmosphere export state on ocn decomp

    call mct_aVect_init(a2x_o, rList=seq_flds_a2x_fields, &
               lsize=MCT_GSMap_lsize(GSMap_ocn, mpicom))
    call mct_aVect_zero(a2x_o)

    ! Initialize av for ice export state on ocn decomp

    call mct_aVect_init(i2x_o, rList=seq_flds_i2x_fields, &
               lsize=mct_GSMap_lsize(GSMap_ocn, mpicom))
    call mct_aVect_zero(i2x_o)

    ! Initialize av for rof export state on ocn decomp

    call mct_aVect_init(r2x_o, rList=seq_flds_r2x_fields, &
               lsize=mct_GSMap_lsize(GSMap_ocn, mpicom))
    call mct_aVect_zero(r2x_o)

  end subroutine mrg_x2o_init_mct

!===========================================================================================

!!  subroutine mrg_x2o_run_mct( cdata_o, a2x_o, i2x_o, r2x_o, xao_o, fractions_o, x2o_o )
  subroutine mrg_x2o_run_mct( cdata_o, a2x_o, i2x_o, xao_o, fractions_o, x2o_o )

    !----------------------------------------------------------------------- 
    !
    ! Arguments
    !
    
    type(seq_cdata), intent(in)    :: cdata_o
    type(mct_aVect), intent(in)    :: a2x_o
    type(mct_aVect), intent(in)    :: i2x_o
!!    type(mct_aVect), intent(in)    :: r2x_o
    type(mct_aVect), intent(in)    :: xao_o
    type(mct_aVect), intent(in)    :: fractions_o
    type(mct_aVect), intent(inout) :: x2o_o
    !
    ! Local variables
    !
    integer  :: n, ki, ko, kir, kor
    integer  :: lsize
    real(r8) :: ifrac,ifracr
    real(r8) :: afrac,afracr
    logical  :: usevector    ! use vector-friendly mct_copy
    real(r8) :: flux_epbalfact
    character(len=cl) :: flux_epbal
    type(seq_infodata_type),pointer :: infodata
    real(r8) :: frac_sum
    real(r8) :: avsdr, anidr, avsdf, anidf   ! albedos
    real(r8) :: fswabsv, fswabsi             ! sw

    !----- formats -----
    character(*),parameter :: F01 =   "('(mrg_x2o_run_mct) ',a,3e11.3,a,f9.6)"
    character(*),parameter :: subName = '(mrg_x2o_run_mct) '
    !----------------------------------------------------------------------- 
    ! 
    ! Copy runoff vector directly
    !
#ifdef CPP_VECTOR
    usevector = .true.
#else
    usevector = .false.
#endif

    call seq_cdata_setptrs(cdata_o, infodata=infodata)
    call mct_aVect_zero(x2o_o)

    call mct_aVect_copy(aVin=a2x_o, aVout=x2o_o, vector=usevector)
    call mct_aVect_copy(aVin=i2x_o, aVout=x2o_o, vector=usevector)
! tcx moved out to a separate accumulate
!!    call mct_aVect_copy(aVin=r2x_o, aVout=x2o_o, vector=usevector)
    call mct_aVect_copy(aVin=xao_o, aVout=x2o_o, vector=usevector)

    call seq_infodata_GetData(infodata, flux_epbalfact = flux_epbalfact)

    ! 
    ! Compute input ocn state (note that this only applies to non-land portion of gridcell)
    !
    ki  = mct_aVect_indexRa(fractions_o,"ifrac",perrWith=subName)
    ko  = mct_aVect_indexRa(fractions_o,"ofrac",perrWith=subName)
    kir = mct_aVect_indexRa(fractions_o,"ifrad",perrWith=subName)
    kor = mct_aVect_indexRa(fractions_o,"ofrad",perrWith=subName)
    lsize = mct_aVect_lsize(x2o_o)
    do n = 1,lsize

       ifrac = fractions_o%rAttr(ki,n)
       afrac = fractions_o%rAttr(ko,n)
       frac_sum = ifrac + afrac
       if ((frac_sum) /= 0._r8) then
          ifrac = ifrac / (frac_sum)
          afrac = afrac / (frac_sum)
       endif

       ifracr = fractions_o%rAttr(kir,n)
       afracr = fractions_o%rAttr(kor,n)
       frac_sum = ifracr + afracr
       if ((frac_sum) /= 0._r8) then
          ifracr = ifracr / (frac_sum)
          afracr = afracr / (frac_sum)
       endif

       x2o_o%rAttr(index_x2o_Foxx_taux ,n) = xao_o%rAttr(index_xao_Faox_taux ,n) * afrac + &
                                             i2x_o%rAttr(index_i2x_Fioi_taux ,n) * ifrac

       x2o_o%rAttr(index_x2o_Foxx_tauy ,n) = xao_o%rAttr(index_xao_Faox_tauy ,n) * afrac + &
                                             i2x_o%rAttr(index_i2x_Fioi_tauy ,n) * ifrac 

       ! --- was flux_solar:
       avsdr = xao_o%rAttr(index_xao_So_avsdr,n)  
       anidr = xao_o%rAttr(index_xao_So_anidr,n)  
       avsdf = xao_o%rAttr(index_xao_So_avsdf,n)  
       anidf = xao_o%rAttr(index_xao_So_anidf,n)  
       fswabsv  =  a2x_o%rAttr(index_a2x_Faxa_swvdr,n) * (1.0_R8 - avsdr) &
                 + a2x_o%rAttr(index_a2x_Faxa_swvdf,n) * (1.0_R8 - avsdf)
       fswabsi  =  a2x_o%rAttr(index_a2x_Faxa_swndr,n) * (1.0_R8 - anidr) &
                 + a2x_o%rAttr(index_a2x_Faxa_swndf,n) * (1.0_R8 - anidf)

       x2o_o%rAttr(index_x2o_Foxx_swnet,n) = (fswabsv + fswabsi)                 * afracr + &
                                             i2x_o%rAttr(index_i2x_Fioi_swpen,n) * ifrac

       x2o_o%rAttr(index_x2o_Foxx_lat  ,n) = xao_o%rAttr(index_xao_Faox_lat  ,n) * afrac

       x2o_o%rAttr(index_x2o_Foxx_sen  ,n) = xao_o%rAttr(index_xao_Faox_sen  ,n) * afrac

       x2o_o%rAttr(index_x2o_Foxx_evap ,n) = xao_o%rAttr(index_xao_Faox_evap ,n) * afrac

       x2o_o%rAttr(index_x2o_Foxx_lwup ,n) = xao_o%rAttr(index_xao_Faox_lwup ,n) * afrac

       x2o_o%rAttr(index_x2o_Foxx_lwdn ,n) = a2x_o%rAttr(index_a2x_Faxa_lwdn ,n) * afrac

       x2o_o%rAttr(index_x2o_Foxx_snow ,n) = a2x_o%rAttr(index_a2x_Faxa_snowc,n) * afrac + &
                                             a2x_o%rAttr(index_a2x_Faxa_snowl,n) * afrac 

       x2o_o%rAttr(index_x2o_Foxx_rain ,n) = a2x_o%rAttr(index_a2x_Faxa_rainc,n) * afrac + &
                                             a2x_o%rAttr(index_a2x_Faxa_rainl,n) * afrac

       x2o_o%rAttr(index_x2o_Foxx_melth,n) = i2x_o%rAttr(index_i2x_Fioi_melth,n) * ifrac

       x2o_o%rAttr(index_x2o_Foxx_meltw,n) = i2x_o%rAttr(index_i2x_Fioi_meltw,n) * ifrac

       x2o_o%rAttr(index_x2o_Foxx_salt ,n) = i2x_o%rAttr(index_i2x_Fioi_salt ,n) * ifrac

       ! scale total precip and runoff by flux_epbalfact (TODO: note in cpl6 this was always 
       ! done over points where imask was > 0 - how does this translate here?)

       x2o_o%rAttr(index_x2o_Foxx_rain ,n) = x2o_o%rAttr(index_x2o_Foxx_rain ,n) * flux_epbalfact
       x2o_o%rAttr(index_x2o_Foxx_snow ,n) = x2o_o%rAttr(index_x2o_Foxx_snow ,n) * flux_epbalfact
! this has been moved into r2xacc_rx
!       x2o_o%rAttr(index_x2o_Forr_roff ,n) = x2o_o%rAttr(index_x2o_Forr_roff ,n) * flux_epbalfact
!       x2o_o%rAttr(index_x2o_Forr_ioff ,n) = x2o_o%rAttr(index_x2o_Forr_ioff ,n) * flux_epbalfact

       x2o_o%rAttr(index_x2o_Foxx_prec ,n) = x2o_o%rAttr(index_x2o_Foxx_rain ,n) + &
                                             x2o_o%rAttr(index_x2o_Foxx_snow ,n)

       x2o_o%rAttr(index_x2o_Foxx_bcphidry,n) = a2x_o%rAttr(index_a2x_Faxa_bcphidry,n) * afrac
       x2o_o%rAttr(index_x2o_Foxx_bcphodry,n) = a2x_o%rAttr(index_a2x_Faxa_bcphodry,n) * afrac
       x2o_o%rAttr(index_x2o_Foxx_bcphiwet,n) = a2x_o%rAttr(index_a2x_Faxa_bcphiwet,n) * afrac
       x2o_o%rAttr(index_x2o_Foxx_ocphidry,n) = a2x_o%rAttr(index_a2x_Faxa_ocphidry,n) * afrac
       x2o_o%rAttr(index_x2o_Foxx_ocphodry,n) = a2x_o%rAttr(index_a2x_Faxa_ocphodry,n) * afrac
       x2o_o%rAttr(index_x2o_Foxx_ocphiwet,n) = a2x_o%rAttr(index_a2x_Faxa_ocphiwet,n) * afrac
       x2o_o%rAttr(index_x2o_Foxx_dstwet1,n)  = a2x_o%rAttr(index_a2x_Faxa_dstwet1,n)  * afrac
       x2o_o%rAttr(index_x2o_Foxx_dstwet2,n)  = a2x_o%rAttr(index_a2x_Faxa_dstwet2,n)  * afrac
       x2o_o%rAttr(index_x2o_Foxx_dstwet3,n)  = a2x_o%rAttr(index_a2x_Faxa_dstwet3,n)  * afrac
       x2o_o%rAttr(index_x2o_Foxx_dstwet4,n)  = a2x_o%rAttr(index_a2x_Faxa_dstwet4,n)  * afrac
       x2o_o%rAttr(index_x2o_Foxx_dstdry1,n)  = a2x_o%rAttr(index_a2x_Faxa_dstdry1,n)  * afrac
       x2o_o%rAttr(index_x2o_Foxx_dstdry2,n)  = a2x_o%rAttr(index_a2x_Faxa_dstdry2,n)  * afrac
       x2o_o%rAttr(index_x2o_Foxx_dstdry3,n)  = a2x_o%rAttr(index_a2x_Faxa_dstdry3,n)  * afrac
       x2o_o%rAttr(index_x2o_Foxx_dstdry4,n)  = a2x_o%rAttr(index_a2x_Faxa_dstdry4,n)  * afrac

    end do
  end subroutine mrg_x2o_run_mct
!
!===========================================================================================
!
  subroutine mrg_x2o_final_mct
    ! ******************
    ! Do nothing for now
    ! ******************
  end subroutine mrg_x2o_final_mct

end module mrg_x2o_mct


