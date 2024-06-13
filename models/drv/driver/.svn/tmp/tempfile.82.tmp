module map_atmatm_mct

!---------------------------------------------------------------------
!
! Purpose:
!
! Collect coupling routines for sequential coupling of atmx - atma
!       
! Author: T. Craig
!
!---------------------------------------------------------------------

  use shr_sys_mod
  use mct_mod
  use seq_comm_mct
  use seq_cdata_mod
  use seq_rearr_mod

  implicit none
  save
  private  ! except

!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------

  public :: map_atm2atm_init_mct
  public :: map_atmx2atma_mct
  public :: map_atma2atmx_mct

  ! Note ccc is component, xxx is driver/coupler

  interface map_atm2atm_init_mct; module procedure &
    map_ccc2ccc_init_mct
  end interface
  interface map_atmx2atma_mct; module procedure &
    map_xxx2ccc_mct
  end interface
  interface map_atma2atmx_mct; module procedure &
    map_ccc2xxx_mct
  end interface

!--------------------------------------------------------------------------
! Public data
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
! Private data
!--------------------------------------------------------------------------

  type(mct_rearr), private :: Re_xxx2ccc
  type(mct_rearr), private :: Re_ccc2xxx

#ifdef CPP_VECTOR
    logical :: usevector = .true.
#else
    logical :: usevector = .false.
#endif

#ifdef SYSUNICOS
    logical :: usealltoall = .true.
#else
    logical :: usealltoall = .false.
#endif

  logical :: copyoption = .true.    ! copy AV if possible
  logical :: gsmapSame  = .false.    ! if gsmaps are same
  character(*),parameter :: subName = '(map_atmatm_mct)'

!=======================================================================
   contains
!=======================================================================

  subroutine map_ccc2ccc_init_mct( cdata_cc, x2c_cc, c2x_cc, ID_cc, &
                                   cdata_cx, x2c_cx, c2x_cx, ID_cx, ID_join)

    implicit none
    !-----------------------------------------------------
    ! 
    ! Arguments
    !
    type(seq_cdata),intent(in)    :: cdata_cc
    type(mct_aVect),intent(inout) :: x2c_cc
    type(mct_aVect),intent(inout) :: c2x_cc
    integer        ,intent(in)    :: ID_cc
    type(seq_cdata),intent(inout) :: cdata_cx
    type(mct_aVect),intent(inout) :: x2c_cx
    type(mct_aVect),intent(inout) :: c2x_cx
    integer        ,intent(in)    :: ID_cx
    integer        ,intent(in)    :: ID_join
    !
    ! Local Variables
    !
    character(len=*),parameter :: subname = "(map_atm2atm_init_mct) "
    type(mct_gsmap),pointer :: gsmap_cc
    type(mct_gsmap),pointer :: gsmap_cx

    !-----------------------------------------------------

    call seq_rearr_init( cdata_cc, x2c_cc, c2x_cc, ID_cc, &
                         cdata_cx, x2c_cx, c2x_cx, ID_cx, ID_join, &
                         Re_ccc2xxx, Re_xxx2ccc)

    call seq_cdata_setptrs(cdata_cc,gsmap=gsmap_cc)
    call seq_cdata_setptrs(cdata_cx,gsmap=gsmap_cx)

    if (seq_rearr_gsmapIdentical(gsmap_cc,gsmap_cx)) then
       gsmapsame = .true.
       if (seq_comm_iamroot(ID_join)) then
          write(logunit,'(2A,L2)') subname,' gsmaps ARE IDENTICAL, copyoption = ',copyoption
       endif
    else
       if (seq_comm_iamroot(ID_join)) write(logunit,'(2A)') subname,' gsmaps are not identical'
       gsmapsame = .false.
    endif

  end subroutine map_ccc2ccc_init_mct

!=======================================================================

  subroutine map_ccc2xxx_mct( cdata_cc, av_cc, cdata_cx, av_cx, fldlist )

    !-----------------------------------------------------
    ! 
    ! Arguments
    !
    type(seq_cdata),intent(in) :: cdata_cc
    type(mct_aVect),intent(in) :: av_cc
    type(seq_cdata),intent(in) :: cdata_cx
    type(mct_aVect),intent(out):: av_cx
    character(len=*),intent(in), optional :: fldlist  ! this is an rList
    !-----------------------------------------------------
    character(len=*),parameter :: subname = "(map_atma2atmx_mct) "
    type(mct_aVect) :: av_test
    integer :: k,n

    if (copyoption .and. gsmapsame) then
       if (present(fldlist)) then
          call mct_aVect_copy(aVin=av_cc,aVout=av_cx,rList=fldlist,vector=usevector)
       else
          call mct_aVect_copy(aVin=av_cc,aVout=av_cx,vector=usevector)
       endif
    else
       if (present(fldlist)) then
          call mct_rearr_rearrange_fldlist(av_cc, av_cx, Re_ccc2xxx, VECTOR=usevector, &
               ALLTOALL=usealltoall, fldlist=fldlist)
       else
          call mct_rearr_rearrange(av_cc, av_cx, Re_ccc2xxx, VECTOR=usevector, ALLTOALL=usealltoall)
       endif
    end if

!tcx verifies rearranging is working properly
#if (1 == 0)
      call mct_avect_init(av_test,av_cc,mct_avect_lsize(av_cc))
      call mct_rearr_rearrange(av_cx, av_test, Re_xxx2ccc, VECTOR=usevector, ALLTOALL=usealltoall)
      do k = 1,mct_avect_nRattr(av_cc)
      do n = 1,mct_avect_lsize(av_cc)
         if (av_cc%rAttr(k,n) /= av_test%rAttr(k,n)) then
            write(6,*) 'tcz r1 ',mct_avect_nRattr(av_cc),mct_avect_nRattr(av_test),mct_avect_lsize(av_cc),mct_avect_lsize(av_test),mct_avect_lsize(av_cx)
            write(6,*) 'tcz diff ',k,n,av_cc%rAttr(k,n),av_test%rAttr(k,n)
            call shr_sys_flush(6)
            call shr_sys_abort()
         endif
      enddo
      enddo
      call mct_avect_clean(av_test)
#endif

  end subroutine map_ccc2xxx_mct

!=======================================================================

  subroutine map_xxx2ccc_mct( cdata_cx, av_cx, cdata_cc, av_cc, fldlist)

    !-----------------------------------------------------
    ! 
    ! Arguments
    !
    type(seq_cdata),intent(in) :: cdata_cx
    type(mct_aVect),intent(in) :: av_cx
    type(seq_cdata),intent(in) :: cdata_cc
    type(mct_aVect),intent(out):: av_cc
    character(len=*),intent(in), optional :: fldlist
    !-----------------------------------------------------
    character(len=*),parameter :: subname = "(map_atmx2atma_mct) "

    if (copyoption .and. gsmapsame) then
       if (present(fldlist)) then
          call mct_aVect_copy(aVin=av_cx,aVout=av_cc,rList=fldlist,vector=usevector)
       else
          call mct_aVect_copy(aVin=av_cx,aVout=av_cc,vector=usevector)
       endif
    else
       if (present(fldlist)) then
          call mct_rearr_rearrange_fldlist(av_cx, av_cc, Re_xxx2ccc, VECTOR=usevector, &
               ALLTOALL=usealltoall, fldlist=fldlist)
       else
          call mct_rearr_rearrange(av_cx, av_cc, Re_xxx2ccc, VECTOR=usevector, ALLTOALL=usealltoall)
       endif
    end if

  end subroutine map_xxx2ccc_mct

!=======================================================================


end module map_atmatm_mct
