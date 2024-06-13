module map_atmlnd_mct

!---------------------------------------------------------------------
!
! Purpose:
!
! Collect coupling routines for sequential coupling of LND-ATM.
!       
!
! Author: R. Jacob, M. Vertenstein
! Revision History:
! 30Mar06 - P. Worley - added optional arguments to MCT_Rearrange call
! 13Apr06 - M. Vertenstein - cleaned up interfaces 
!
!---------------------------------------------------------------------

  use shr_kind_mod      ,only: R8 => SHR_KIND_R8, IN=>SHR_KIND_IN
  use shr_sys_mod
  use shr_const_mod
  use shr_mct_mod, only: shr_mct_sMatPInitnc
  use mct_mod

  use seq_comm_mct, only : logunit, loglevel
  use seq_cdata_mod
  use seq_flds_indices
  use seq_infodata_mod
  use m_die

  implicit none
  save
  private  ! except

!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------

  public :: map_lnd2atm_init_mct
  public :: map_atm2lnd_init_mct
  public :: map_lnd2atm_mct
  public :: map_atm2lnd_mct

!--------------------------------------------------------------------------
! Public data
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
! Private data
!--------------------------------------------------------------------------

  type(mct_rearr), private :: Re_lnd2atm
  type(mct_rearr), private :: Re_atm2lnd
  type(mct_sMatp), private :: sMatp_Fa2l
  type(mct_sMatp), private :: sMatp_Sa2l
  type(mct_sMatp), private :: sMatp_Fl2a
  type(mct_sMatp), private :: sMatp_Sl2a

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
  logical, private :: samegrid_mapa2l

!=======================================================================
   contains
!=======================================================================

  subroutine map_atm2lnd_init_mct( cdata_a, cdata_l)

    !--------------------------------------------------
    ! 
    ! Arguments
    !
    type(seq_cdata),intent(in) :: cdata_a
    type(seq_cdata),intent(in) :: cdata_l
    ! 
    ! Local variables
    !
    type(seq_infodata_type), pointer :: infodata
    type(mct_gsMap), pointer :: gsMap_a           ! atm gsMap
    type(mct_gsMap), pointer :: gsMap_l           ! lnd gsMap
    type(mct_gGrid), pointer :: dom_a             ! atm domain
    type(mct_gGrid), pointer :: dom_l             ! lnd domain
    integer                  :: mpicom            ! communicator spanning atm and lnd
    integer                  :: ka, km            ! indices
    integer                  :: lsize             ! size of attribute vector
    type(mct_aVect)          :: areasrc           ! atm areas from mapping file
    type(mct_aVect)          :: areadst           ! lnd areas set to atm areas
    character(*),parameter :: subName = '(map_atm2lnd_init_mct) '
    !--------------------------------------------------

    call seq_cdata_setptrs(cdata_a, gsMap=gsMap_a, dom=dom_a)
    call seq_cdata_setptrs(cdata_l, gsMap=gsMap_l, dom=dom_l)
    call seq_cdata_setptrs(cdata_a, mpicom=mpicom, infodata=infodata)

    call seq_infodata_GetData( infodata, samegrid_al=samegrid_mapa2l)

    lsize = mct_gsMap_lsize(gsMap_a, mpicom)
    call mct_aVect_init( areasrc, rList="aream", lsize=lsize )
       
    lsize = mct_gsMap_lsize(gsMap_l, mpicom)
    call mct_aVect_init( areadst, rList="aream", lsize=lsize )

    if (samegrid_mapa2l) then

       call mct_rearr_init(gsMap_a, gsMap_l, mpicom, Re_atm2lnd)

       ! --- copy atm area into land aream

       ka = mct_aVect_indexRA(areasrc   , "aream")
       km = mct_aVect_indexRa(dom_a%data, "aream" )
       areasrc%rAttr(ka,:) = dom_a%data%rAttr(km,:)
       call mct_rearr_rearrange(areasrc, areadst, Re_atm2lnd, VECTOR=usevector, &
          ALLTOALL=usealltoall)

    else

       call shr_mct_sMatPInitnc(sMatp_Fa2l,gsMap_a,gsMap_l,"seq_maps.rc", &
          "atm2lndFmapname:","atm2lndFmaptype:",mpicom, &
          areasrc=areasrc, areadst=areadst)

       call shr_mct_sMatPInitnc(sMatp_Sa2l,gsMap_a,gsMap_l,"seq_maps.rc", &
          "atm2lndSmapname:","atm2lndSmaptype:",mpicom)

    endif

    ka = mct_aVect_indexRA(areadst   ,"aream")
    km = mct_aVect_indexRA(dom_l%data,"aream")
    dom_l%data%rAttr(km,:) = areadst%rAttr(ka,:)

    call mct_aVect_clean(areasrc)      
    call mct_aVect_clean(areadst)      

  end subroutine map_atm2lnd_init_mct

!=======================================================================

  subroutine map_lnd2atm_init_mct( cdata_l, cdata_a )

    !--------------------------------------------------
    ! 
    ! Arguments
    !
    type(seq_cdata),intent(in) :: cdata_l
    type(seq_cdata),intent(in) :: cdata_a
    ! 
    ! Local variables
    !
    type(seq_infodata_type), pointer :: infodata
    type(mct_gsMap), pointer :: gsMap_a           ! atm gsMap
    type(mct_gsMap), pointer :: gsMap_l           ! lnd gsMap
    type(mct_gGrid), pointer :: dom_l             ! lnd domain
    integer                  :: kf,iam,ierr,lsize
    integer                  :: mpicom            ! communicator spanning atm and lnd
    character(*),parameter :: subName = '(map_lnd2atm_init_mct) '
    !--------------------------------------------------

    call seq_cdata_setptrs(cdata_l, gsMap=gsMap_l, dom=dom_l)
    call seq_cdata_setptrs(cdata_a, gsMap=gsMap_a)
    call seq_cdata_setptrs(cdata_a, mpicom=mpicom,infodata=infodata)
    call mpi_comm_rank(mpicom,iam,ierr)
    
    call seq_infodata_GetData( infodata, samegrid_al=samegrid_mapa2l)

    ! Initialize lnd -> atm mapping or rearranger

    if (samegrid_mapa2l) then

       call mct_rearr_init(gsMap_l, gsMap_a, mpicom, Re_lnd2atm)

    else

       call shr_mct_sMatPInitnc(sMatp_Fl2a, gsMap_l, gsMap_a, "seq_maps.rc", &
            "lnd2atmFmapname:", "lnd2atmFmaptype:", mpicom)

       call shr_mct_sMatPInitnc(sMatp_Sl2a, gsMap_l, gsMap_a, "seq_maps.rc", &
            "lnd2atmSmapname:", "lnd2atmSmaptype:", mpicom)

    endif

  end subroutine map_lnd2atm_init_mct

!=======================================================================

  subroutine map_atm2lnd_mct( cdata_a, av_a, cdata_l, av_l, fluxlist, statelist )

    !--------------------------------------------------
    ! 
    ! Arguments
    !
    type(seq_cdata), intent(in)  :: cdata_a
    type(mct_aVect), intent(in)  :: av_a
    type(seq_cdata), intent(in)  :: cdata_l
    type(mct_aVect), intent(out) :: av_l
    character(len=*),intent(in), optional :: statelist
    character(len=*),intent(in), optional :: fluxlist
    !
    ! Local
    ! 
    integer                :: lsize
    type(mct_aVect)        :: av_l_f     ! temporary flux attribute vector
    type(mct_aVect)        :: av_l_s     ! temporary state attribute vector
    character(*),parameter :: subName = '(map_atm2lnd_mct) '
    !--------------------------------------------------

    if (samegrid_mapa2l) then

       if (present(fluxlist) .or. present(statelist)) then
          if (present(fluxlist)) then
             call mct_rearr_rearrange_fldlist(av_a, av_l, Re_atm2lnd, VECTOR=usevector, &
                  ALLTOALL=usealltoall, fldlist=fluxlist)
          endif
          if (present(statelist)) then
              call mct_rearr_rearrange_fldlist(av_a, av_l, Re_atm2lnd, VECTOR=usevector, &
                  ALLTOALL=usealltoall, fldlist=statelist)
          endif
       else
          call mct_rearr_rearrange(av_a, av_l, Re_atm2lnd, VECTOR=usevector, ALLTOALL=usealltoall)
       end if

    else

       if (present(fluxlist) .or. present(statelist)) then
          if (present(fluxlist)) then
             lsize = mct_aVect_lsize(av_l)
             call mct_aVect_init (av_l_f, rlist=fluxlist , lsize=lsize)
             call mct_sMat_avMult(av_a, sMatp_Fa2l, av_l_f, VECTOR=usevector, rList=fluxlist)
             call mct_aVect_copy (aVin=av_l_f, aVout=av_l, vector=usevector)
             call mct_aVect_clean(av_l_f)
          end if
          if (present(statelist)) then
             lsize = mct_aVect_lsize(av_l)
             call mct_aVect_init (av_l_s, rlist=statelist, lsize=lsize)
             call mct_sMat_avMult(av_a, sMatp_Sa2l, av_l_s, VECTOR=usevector, rList=statelist)
             call mct_aVect_copy (aVin=av_l_s, aVout=av_l, vector=usevector)
             call mct_aVect_clean(av_l_s)
          end if
       else
          !--- default is flux mapping
          call mct_sMat_avMult(av_a, sMatp_Fa2l, av_l, VECTOR=usevector)
       endif
    endif
       
  end subroutine map_atm2lnd_mct

!=======================================================================

  subroutine map_lnd2atm_mct( cdata_l, av_l, cdata_a, av_a, &
                              fractions_l, fractions_a, &
                              fluxlist, statelist )

    !--------------------------------------------------
    ! 
    ! Arguments
    !
    type(seq_cdata) ,intent(in)  :: cdata_l
    type(mct_aVect) ,intent(in)  :: av_l
    type(seq_cdata) ,intent(in)  :: cdata_a
    type(mct_aVect) ,intent(out) :: av_a
    type(mct_aVect) ,intent(in), optional :: fractions_l
    type(mct_aVect) ,intent(in), optional :: fractions_a
    character(len=*),intent(in), optional :: fluxlist
    character(len=*),intent(in), optional :: statelist
    !
    ! Local
    !
    integer  :: i,j,kl,lsize,numats,ier
    real(r8) :: recip
    type(mct_aVect)          :: av_a_f     ! temporary flux attribute vector
    type(mct_aVect)          :: av_a_s     ! temporary state attribute vector
    type(mct_aVect)          :: temp       ! temporary attribute vector
    character(*),parameter :: subName = '(map_lnd2atm_mct) '
    !--------------------------------------------------

    if (samegrid_mapa2l) then

       if (present(fluxlist) .or. present(statelist)) then
          if (present(fluxlist)) then
             call mct_rearr_rearrange_fldlist(av_l, av_a, Re_lnd2atm, VECTOR=usevector, &
                  ALLTOALL=usealltoall, fldlist=fluxlist)
          endif
          if (present(statelist)) then
             call mct_rearr_rearrange_fldlist(av_l, av_a, Re_lnd2atm, VECTOR=usevector, &
                  ALLTOALL=usealltoall, fldlist=statelist)
          endif
       else
          call mct_rearr_rearrange(av_l, av_a, Re_lnd2atm, VECTOR=usevector, ALLTOALL=usealltoall)
       end if

    else

       ! Normalize input data with fraction of land on land grid

       if (present(fractions_l) .and. present(fractions_a)) then

          lsize = mct_aVect_lsize(av_l)
          numats = mct_aVect_nRAttr(av_l)
          kl = mct_aVect_indexRA(fractions_l,"lfrin")
          call mct_aVect_init(temp, av_l, lsize=lsize)
          do j=1,lsize
             do i=1,numats
                temp%rAttr(i,j) = av_l%rAttr(i,j)* fractions_l%rAttr(kl,j)
             end do
          end do

          ! Perform remapping

          if (present(fluxlist) .or. present(statelist)) then
             if (present (fluxlist)) then
                lsize = mct_aVect_lsize(av_a)
                call mct_aVect_init (av_a_f, rlist=fluxlist , lsize=lsize)
                call mct_sMat_avMult(temp, sMatp_Fl2a, av_a_f, VECTOR=usevector, rList=fluxlist )
                call mct_aVect_copy (aVin=av_a_f, aVout=av_a, vector=usevector)
                call mct_aVect_clean(av_a_f)
             end if
             if (present(statelist)) then
                lsize = mct_aVect_lsize(av_a)
                call mct_aVect_init (av_a_s, rlist=statelist, lsize=lsize)
                call mct_sMat_avMult(temp, sMatp_Sl2a, av_a_s, VECTOR=usevector, rList=statelist)
                call mct_aVect_copy (aVin=av_a_s, aVout=av_a, vector=usevector)
                call mct_aVect_clean(av_a_s)
             end if
          else
             ! --- default is flux mapping
             call mct_sMat_avMult(temp, sMatp_Fl2a, av_a_f, VECTOR=usevector)
          endif

          ! Clean up temporary vector
          call mct_aVect_clean(temp)

          ! Denormalize output data with fraction of land on atmosphere grid 
          lsize = mct_aVect_lsize(av_a)
          numats = mct_aVect_nRAttr(av_a)

          kl    = mct_aVect_indexRA(fractions_a, "lfrin")
          do j=1,lsize
             if (fractions_a%rAttr(kl,j) /= 0.0_R8) then
                recip = 1.0_R8 / fractions_a%rAttr(kl,j)
             else
                recip = 0.0_R8
             end if
             do i =1,numats
                   av_a%rAttr(i,j) = av_a%rAttr(i,j) * recip
             end do
          end do

       else

          if (present(fluxlist) .or. present(statelist)) then
             if (present (fluxlist)) then
                lsize = mct_aVect_lsize(av_a)
                call mct_aVect_init (av_a_f, rlist=fluxlist , lsize=lsize)
                call mct_sMat_avMult(av_l, sMatp_Fl2a, av_a_f, VECTOR=usevector, rList=fluxlist )
                call mct_aVect_copy (aVin=av_a_f, aVout=av_a, vector=usevector)
                call mct_aVect_clean(av_a_f)
             end if
             if (present(statelist)) then
                lsize = mct_aVect_lsize(av_a)
                call mct_aVect_init (av_a_s, rlist=statelist, lsize=lsize)
                call mct_sMat_avMult(av_l, sMatp_Sl2a, av_a_s, VECTOR=usevector, rList=statelist)
                call mct_aVect_copy (aVin=av_a_s, aVout=av_a, vector=usevector)
                call mct_aVect_clean(av_a_s)
             end if
          else
             ! --- default is flux mapping
             call mct_sMat_avMult(av_l, sMatp_Fl2a, av_a, VECTOR=usevector)
          endif

       endif

    endif

  end subroutine map_lnd2atm_mct

end module map_atmlnd_mct
