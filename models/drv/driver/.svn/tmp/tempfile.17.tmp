module mrg_x2l_mct

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

  public :: mrg_x2l_init_mct
  public :: mrg_x2l_run_mct
  public :: mrg_x2l_final_mct

!--------------------------------------------------------------------------
! Private data
!--------------------------------------------------------------------------

!===========================================================================================
contains
!===========================================================================================

  subroutine mrg_x2l_init_mct( cdata_l, a2x_l )

    type(seq_cdata) ,intent(in)     :: cdata_l
    type(mct_aVect), intent(inout)  :: a2x_l

    type(mct_GsMap), pointer        :: GSMap_lnd
    integer                         :: mpicom

    ! Set gsMap
    call seq_cdata_setptrs(cdata_l, gsMap=gsMap_lnd, mpicom=mpicom)

    ! Initialize av for atmosphere export state on lnd decomp

    call mct_aVect_init(a2x_l, rList=seq_flds_a2x_fields, &
         lsize=mct_gsMap_lsize(gsMap_lnd, mpicom))
    call mct_aVect_zero(a2x_l)

  end subroutine mrg_x2l_init_mct

!===========================================================================================

  subroutine mrg_x2l_run_mct( cdata_l, a2x_l, x2l_l )

    !----------------------------------------------------------------------- 
    !
    ! Arguments
    !
    type(seq_cdata), intent(in)     :: cdata_l
    type(mct_aVect), intent(inout)  :: a2x_l  ! input
    type(mct_aVect), intent(inout)  :: x2l_l  ! output
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

    call mct_aVect_copy(aVin=a2x_l, aVout=x2l_l, vector=usevector)

  end subroutine mrg_x2l_run_mct
!
!===========================================================================================
!
  subroutine mrg_x2l_final_mct
    ! ******************
    ! Do nothing for now
    ! ******************
  end subroutine mrg_x2l_final_mct

end module mrg_x2l_mct


