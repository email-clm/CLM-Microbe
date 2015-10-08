module map_snoglc_mct

!---------------------------------------------------------------------
!
! Purpose:
!
! Collect coupling routines for sequential coupling of SNO-GLC.
!       
! Author: R. Jacob, M. Vertenstein
!
!---------------------------------------------------------------------

  use shr_sys_mod
  use mct_mod
  use seq_cdata_mod
  use seq_comm_mct

  implicit none
  save
  private  ! except

!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------

  public :: map_sno2glc_init_mct
  public :: map_glc2sno_init_mct
  public :: map_sno2glc_mct
  public :: map_glc2sno_mct

!--------------------------------------------------------------------------
! Public data
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
! Private data
!--------------------------------------------------------------------------

  type(mct_rearr), private :: Re_sno2glc
  type(mct_rearr), private :: Re_glc2sno

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

  character(*),parameter :: subName = '(map_snoglc_mct)'

!=======================================================================
   contains
!=======================================================================

  subroutine map_glc2sno_init_mct( cdata_g, cdata_s)

    !-----------------------------------------------------
    ! 
    ! Arguments
    !
    type(seq_cdata),intent(in) :: cdata_g
    type(seq_cdata),intent(in) :: cdata_s
    !
    ! Local Variables
    !
    integer                  :: glcsize, snosize  ! global grid sizes      
    type(mct_gsMap), pointer :: gsMap_s           ! sno gsMap
    type(mct_gsMap), pointer :: gsMap_g           ! glc gsMap
    integer                  :: mpicom            ! communicator spanning glc and sno
    !-----------------------------------------------------

    call seq_cdata_setptrs(cdata_g, gsMap=gsMap_g)
    call seq_cdata_setptrs(cdata_s, gsMap=gsMap_s)
    call seq_cdata_setptrs(cdata_g, mpicom=mpicom)

    ! Sanity check

    glcsize = mct_gsMap_gsize(gsMap_g)
    snosize = mct_gsMap_gsize(gsMap_s)
    if (glcsize /= snosize) then
      write(logunit,*) "(map_glc2sno_init_mct) sno and glc are different."
      write(logunit,*) "(map_glc2sno_init_mct) Must be same size. Exiting."
      call shr_sys_abort(subName // "different size")
    endif

    ! Initialize rearranger

    call mct_rearr_init(gsMap_g, gsMap_s, mpicom, Re_glc2sno)

  end subroutine map_glc2sno_init_mct

!=============================================================

  subroutine map_sno2glc_init_mct( cdata_s, cdata_g)

    !---------------------------------------------
    ! 
    ! Arguments
    !
    type(seq_cdata),intent(in) :: cdata_s
    type(seq_cdata),intent(in) :: cdata_g
    !
    ! Local Variables
    !
    integer                  :: glcsize, snosize  ! global grid sizes
    type(mct_gsMap), pointer :: gsMap_s           ! sno gsMap
    type(mct_gsMap), pointer :: gsMap_g           ! glc gsMap
    integer                  :: mpicom            ! communicator spanning glc and sno
    integer                  :: ka, km            ! indsnos
    type(mct_gGrid), pointer :: dom_s             ! sno domain
    type(mct_gGrid), pointer :: dom_g             ! glc domain
    integer                  :: lsize             ! size of attribute vector
    type(mct_aVect)          :: areasrc           ! glc areas from mapping file
    type(mct_aVect)          :: areadst           ! sno areas from mapping file
    !---------------------------------------------

    call seq_cdata_setptrs(cdata_s, gsMap=gsMap_s, dom=dom_s) 
    call seq_cdata_setptrs(cdata_g, gsMap=gsMap_g, dom=dom_g)
    call seq_cdata_setptrs(cdata_g, mpicom=mpicom)

    ! Sanity check

    glcsize = mct_gsMap_gsize(gsMap_g)
    snosize = mct_gsMap_gsize(gsMap_s)
    if (glcsize /= snosize) then
      write(logunit,*) "(map_sno2glc_init_mct) glc and sno grids are different."
      write(logunit,*) "(map_sno2glc_init_mct) Must be the same....Exiting."
      call shr_sys_abort(subName // "different glc,sno grids")
    endif

    ! Initialize rearranger

    call mct_rearr_init(gsMap_s, gsMap_g, mpicom, Re_sno2glc)

    ! Set glc aream to sno aream
    ! Note that "aream" attribute of dom_s%data is set in routine map_atm2sno_init_mct

    lsize = mct_gsMap_lsize(gsMap_s, mpicom)
    call mct_aVect_init( areasrc, rList="aream", lsize=lsize )

    lsize = mct_gsMap_lsize(gsMap_g, mpicom)
    call mct_aVect_init( areadst, rList="aream", lsize=lsize )
       
    km = mct_aVect_indexRA(dom_s%data,"aream")
    ka = mct_aVect_indexRA(areasrc   ,"aream")
    areasrc%rAttr(ka,:) = dom_s%data%rAttr(km,:)

    ! tcraig, 6/2010, initialize areadst with dom_g km because areasrc is only over land
    ka = mct_aVect_indexRA(areadst   ,"aream")
    km = mct_aVect_indexRA(dom_g%data,"aream")
    areadst%rAttr(ka,:) = dom_g%data%rAttr(km,:)

    call mct_rearr_rearrange(areasrc, areadst, Re_sno2glc, VECTOR=usevector, ALLTOALL=usealltoall)

    ka = mct_aVect_indexRA(areadst   ,"aream")
    km = mct_aVect_indexRA(dom_g%data,"aream")
    dom_g%data%rAttr(km,:) = areadst%rAttr(ka,:)

    call mct_aVect_clean(areasrc)      
    call mct_aVect_clean(areadst)      

  end subroutine map_sno2glc_init_mct  

!=======================================================================

  subroutine map_glc2sno_mct( cdata_g, av_g, cdata_s, av_s, fluxlist, statelist )

    !-----------------------------------------------------
    ! 
    ! Arguments
    !
    type(seq_cdata),intent(in) :: cdata_g
    type(mct_aVect),intent(in) :: av_g
    type(seq_cdata),intent(in) :: cdata_s
    type(mct_aVect),intent(out):: av_s
    character(len=*),intent(in), optional :: fluxlist
    character(len=*),intent(in), optional :: statelist
    !-----------------------------------------------------

    if (present(fluxlist) .or. present(statelist)) then
       if (present(fluxlist)) then
          call mct_rearr_rearrange_fldlist(av_g, av_s, Re_glc2sno, VECTOR=usevector, &
               ALLTOALL=usealltoall, fldlist=fluxlist)
       endif
       if (present(statelist)) then
          call mct_rearr_rearrange_fldlist(av_g, av_s, Re_glc2sno, VECTOR=usevector, &
               ALLTOALL=usealltoall, fldlist=statelist)
       endif
    else
       call mct_rearr_rearrange(av_g, av_s, Re_glc2sno, VECTOR=usevector, ALLTOALL=usealltoall)
    end if

  end subroutine map_glc2sno_mct

!=======================================================================

  subroutine map_sno2glc_mct( cdata_s, av_s, cdata_g, av_g, fluxlist, statelist)

    !-----------------------------------------------------
    ! 
    ! Arguments
    !
    type(seq_cdata),intent(in) :: cdata_s
    type(mct_aVect),intent(in) :: av_s
    type(seq_cdata),intent(in) :: cdata_g
    type(mct_aVect),intent(out):: av_g
    character(len=*),intent(in), optional :: fluxlist
    character(len=*),intent(in), optional :: statelist
    !-----------------------------------------------------

    if (present(fluxlist) .or. present(statelist)) then
       if (present(fluxlist)) then
          call mct_rearr_rearrange_fldlist(av_s, av_g, Re_sno2glc, VECTOR=usevector, &
               ALLTOALL=usealltoall, fldlist=fluxlist)
       endif
       if (present(statelist)) then
          call mct_rearr_rearrange_fldlist(av_s, av_g, Re_sno2glc, VECTOR=usevector, &
               ALLTOALL=usealltoall, fldlist=statelist)
       endif
    else
       call mct_rearr_rearrange(av_s, av_g, Re_sno2glc, VECTOR=usevector, ALLTOALL=usealltoall)
    end if

  end subroutine map_sno2glc_mct

end module map_snoglc_mct
