module mrg_x2i_mct

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

  public :: mrg_x2i_init_mct
  public :: mrg_x2i_run_mct
  public :: mrg_x2i_final_mct

!--------------------------------------------------------------------------
! Private data
!--------------------------------------------------------------------------

!===========================================================================================
contains
!===========================================================================================

  subroutine mrg_x2i_init_mct( cdata_i, a2x_i, o2x_i )

    type(seq_cdata), intent(in)     :: cdata_i
    type(mct_aVect), intent(inout)  :: a2x_i
    type(mct_aVect), intent(inout)  :: o2x_i

    type(mct_gsMap), pointer        :: gsMap_ice
    integer                         :: mpicom

    ! Set gsMap
    call seq_cdata_setptrs(cdata_i,gsMap=gsMap_ice,mpicom=mpicom)

    ! Initialize av for atmosphere export state on ice decomp

    call MCT_aVect_init(a2x_i, rList=seq_flds_a2x_fields, &
               lsize=mct_gsMap_lsize(gsMap_ice, mpicom))
    call MCT_aVect_zero(a2x_i)

    ! Initialize av for ocn export state on ice decomp

    call MCT_aVect_init(o2x_i, rList=seq_flds_o2x_fields, &
               lsize=mct_gsMap_lsize(gsMap_ice, mpicom))
    call MCT_aVect_zero(o2x_i)

  end subroutine mrg_x2i_init_mct
        
!===========================================================================================

  subroutine mrg_x2i_run_mct( cdata_i, a2x_i, o2x_i, x2i_i )

    !----------------------------------------------------------------------- 
    !
    ! Arguments
    !
    type(seq_cdata),intent(in) :: cdata_i
    type(mct_aVect),intent(in) :: a2x_i
    type(mct_aVect),intent(in) :: o2x_i
    type(mct_aVect),intent(out):: x2i_i
    !
    ! Local variables
    !
    logical :: usevector    ! use vector-friendly mct_copy
    integer :: i
    real(r8):: flux_epbalfact
    character(len=cl) :: flux_epbal
    type(seq_infodata_type),pointer :: infodata
    !
    !----- formats -----
    character(*),parameter :: F01 =   "('(mrg_x2i_run_mct) ',a,3e11.3,a,f9.6)"
    character(*),parameter :: subName = '(mrg_x2i_run_mct) '
    !----------------------------------------------------------------------- 
    ! 
    ! Combine atm/ocn states to compute input ice state
    !
#ifdef CPP_VECTOR
    usevector = .true.
#else
    usevector = .false.
#endif

    call seq_cdata_setptrs(cdata_i,infodata=infodata)

    call mct_aVect_copy(aVin=o2x_i, aVout=x2i_i, vector=usevector)
    call mct_aVect_copy(aVin=a2x_i, aVout=x2i_i, vector=usevector)

    ! Apply correction to precipitation of requested driver namelist

    call seq_infodata_GetData(infodata, flux_epbalfact = flux_epbalfact)

    ! Merge total snow and precip for ice input

    do i = 1,mct_aVect_lsize(x2i_i)
       x2i_i%rAttr(index_x2i_Faxa_rain,i) = a2x_i%rAttr(index_a2x_Faxa_rainc,i) + &
	                                    a2x_i%rAttr(index_a2x_Faxa_rainl,i)
       x2i_i%rAttr(index_x2i_Faxa_snow,i) = a2x_i%rAttr(index_a2x_Faxa_snowc,i) + &
	                                    a2x_i%rAttr(index_a2x_Faxa_snowl,i) 

       ! scale total precip and runoff by flux_epbalfact (TODO: note in cpl6 this was always 
       ! done over points where imask was > 0 - how does this translate here?)

       x2i_i%rAttr(index_x2i_Faxa_rain,i) = x2i_i%rAttr(index_x2i_Faxa_rain,i) * flux_epbalfact
       x2i_i%rAttr(index_x2i_Faxa_snow,i) = x2i_i%rAttr(index_x2i_Faxa_snow,i) * flux_epbalfact
       	
    end do

  end subroutine mrg_x2i_run_mct
!
!===========================================================================================
!
  subroutine mrg_x2i_final_mct
    ! ******************
    ! Do nothing for now
    ! ******************
  end subroutine mrg_x2i_final_mct

end module mrg_x2i_mct


