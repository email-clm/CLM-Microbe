module mrg_x2s_mct

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

  public :: mrg_x2s_init_mct
  public :: mrg_x2s_run_mct
  public :: mrg_x2s_final_mct

!--------------------------------------------------------------------------
! Private data
!--------------------------------------------------------------------------

!===========================================================================================
contains
!===========================================================================================

  subroutine mrg_x2s_init_mct( cdata_s, g2x_s)

    type(seq_cdata) ,intent(in)     :: cdata_s
    type(mct_aVect), intent(inout)  :: g2x_s

    type(mct_GsMap), pointer        :: GSMap_sno
    integer                         :: mpicom

    ! Set gsMap
    call seq_cdata_setptrs(cdata_s, gsMap=gsMap_sno, mpicom=mpicom)

    ! Initialize av for atmosphere export state on lnd decomp

    call mct_aVect_init(g2x_s, rList=seq_flds_g2x_fields, &
         lsize=mct_gsMap_lsize(gsMap_sno, mpicom))
    call mct_aVect_zero(g2x_s)

  end subroutine mrg_x2s_init_mct

!===========================================================================================

  subroutine mrg_x2s_run_mct( cdata_s, g2x_s, x2s_s )

    !----------------------------------------------------------------------- 
    !
    ! Arguments
    !
    type(seq_cdata), intent(in)     :: cdata_s
    type(mct_aVect), intent(inout)  :: g2x_s  ! input
    type(mct_aVect), intent(inout)  :: x2s_s  ! output
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

    call mct_aVect_copy(aVin=g2x_s, aVout=x2s_s, vector=usevector)

  end subroutine mrg_x2s_run_mct
!
!===========================================================================================
!
  subroutine mrg_x2s_final_mct
    ! ******************
    ! Do nothing for now
    ! ******************
  end subroutine mrg_x2s_final_mct

end module mrg_x2s_mct


