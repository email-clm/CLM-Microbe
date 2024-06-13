module mrg_x2g_mct

  use shr_kind_mod
  use mct_mod
  use seq_flds_mod
  use seq_flds_indices
  use seq_comm_mct
  use seq_cdata_mod

  implicit none
  save
  private ! except

!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------

  public :: mrg_x2g_init_mct
  public :: mrg_x2g_run_mct
  public :: mrg_x2g_final_mct

!--------------------------------------------------------------------------
! Private data
!--------------------------------------------------------------------------

!===========================================================================================
contains
!===========================================================================================

  subroutine mrg_x2g_init_mct( cdata_g, s2x_g)

    type(seq_cdata) ,intent(in)     :: cdata_g
    type(mct_aVect), intent(inout)  :: s2x_g

    type(mct_GsMap), pointer        :: GSMap_glc
    integer                         :: mpicom

    ! Set gsMap
    call seq_cdata_setptrs(cdata_g, gsMap=gsMap_glc, mpicom=mpicom)

    ! Initialize av for atmosphere export state on lnd decomp

    call mct_aVect_init(s2x_g, rList=seq_flds_s2x_fields, &
         lsize=mct_gsMap_lsize(gsMap_glc, mpicom))
    call mct_aVect_zero(s2x_g)

  end subroutine mrg_x2g_init_mct

!===========================================================================================

  subroutine mrg_x2g_run_mct( cdata_g, s2x_g, x2g_g )

    !----------------------------------------------------------------------- 
    !
    ! Arguments
    !
    type(seq_cdata), intent(in)     :: cdata_g
    type(mct_aVect), intent(inout)  :: s2x_g  ! input
    type(mct_aVect), intent(inout)  :: x2g_g  ! output
    !
    ! Local variables
    !
    logical :: usevector    ! use vector-friendly mct_copy
    !----------------------------------------------------------------------- 
    ! 
    ! Create input land state directly from atm output state
    !
#ifdef CPP_VECTOR
   usevector = .true.
#else
   usevector = .false.
#endif

    call mct_aVect_copy(aVin=s2x_g, aVout=x2g_g, vector=usevector)

  end subroutine mrg_x2g_run_mct
!
!===========================================================================================
!
  subroutine mrg_x2g_final_mct
    ! ******************
    ! Do nothing for now
    ! ******************
  end subroutine mrg_x2g_final_mct

end module mrg_x2g_mct


