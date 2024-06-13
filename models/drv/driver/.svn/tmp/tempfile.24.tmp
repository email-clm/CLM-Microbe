module map_atmocn_mct

!---------------------------------------------------------------------
!
! Purpose:
!
! Collect coupling routines for sequential coupling of OCN-ATM.
!       
! Author: R. Jacob, M. Vertenstein
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

  public :: map_ocn2atm_init_mct
  public :: map_atm2ocn_init_mct
  public :: map_ocn2atm_mct
  public :: map_atm2ocn_mct

!--------------------------------------------------------------------------
! Public data
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
! Private data
!--------------------------------------------------------------------------

  type(mct_rearr), private :: Re_ocn2atm
  type(mct_rearr), private :: Re_atm2ocn
  type(mct_sMatp), private :: sMatp_Fa2o
  type(mct_sMatp), private :: sMatp_Sa2o
  type(mct_sMatp), private :: sMatp_Fo2a
  type(mct_sMatp), private :: sMatp_So2a

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
  logical, private :: samegrid_mapa2o
  
  integer         :: ni_a   ! number of longitudes on input grid
  integer         :: nj_a   ! number of latitudes  on input grid
  integer         :: ni_o   ! number of longitudes on output grid
  integer         :: nj_o   ! number of latitudes  on output grid
   
!=======================================================================
contains
!=======================================================================

  subroutine map_atm2ocn_init_mct( cdata_a, cdata_o)

    !-----------------------------------------------------
    ! 
    ! Arguments
    !
    type(seq_cdata),intent(in) :: cdata_a
    type(seq_cdata),intent(in) :: cdata_o
    !
    ! Local Variables
    !
    integer                  :: ka, km            ! indices
    integer                  :: lsize             ! size of attribute vector
    type(seq_infodata_type), pointer :: infodata
    type(mct_gsMap), pointer :: gsMap_a
    type(mct_gsMap), pointer :: gsMap_o
    type(mct_ggrid), pointer :: dom_a
    type(mct_ggrid), pointer :: dom_o
    type(mct_aVect)          :: areasrc           ! ocn areas from mapping file
    type(mct_aVect)          :: areadst           ! atm areas from mapping file
    integer                  :: mpicom
    character(*),parameter :: subName = '(map_atm2ocn_init_mct) '
    !-----------------------------------------------------

    call seq_cdata_setptrs(cdata_a, gsMap=gsMap_a, dom=dom_a)
    call seq_cdata_setptrs(cdata_o, gsMap=gsMap_o, dom=dom_o)
    call seq_cdata_setptrs(cdata_a, mpicom=mpicom, infodata=infodata)

    call seq_infodata_GetData( infodata, samegrid_ao=samegrid_mapa2o)

    if (samegrid_mapa2o) then

       call mct_rearr_init(gsMap_a, gsMap_o, mpicom, Re_atm2ocn)

       !--- want to copy atm area into atm and ocn aream ---

       lsize = mct_gsMap_lsize(gsMap_a, mpicom)
       call mct_aVect_init( areasrc, rList="aream", lsize=lsize )
       
       lsize = mct_gsMap_lsize(gsMap_o, mpicom)
       call mct_aVect_init( areadst, rList="aream", lsize=lsize )
       
       ka = mct_aVect_indexRa(dom_a%data, "area" )
       km = mct_aVect_indexRA(areasrc   , "aream")
       areasrc%rAttr(km,:) = dom_a%data%rAttr(ka,:)

       call mct_rearr_rearrange(areasrc, areadst, Re_atm2ocn, VECTOR=usevector, &
          ALLTOALL=usealltoall)

       ka = mct_aVect_indexRA(areasrc   ,"aream")
       km = mct_aVect_indexRA(dom_a%data,"aream")
       dom_a%data%rAttr(km,:) = areasrc%rAttr(ka,:)

       ka = mct_aVect_indexRA(areadst   ,"aream")
       km = mct_aVect_indexRA(dom_o%data,"aream")
       dom_o%data%rAttr(km,:) = areadst%rAttr(ka,:)

       call mct_aVect_clean(areasrc)
       call mct_aVect_clean(areadst)      

    else      

       call shr_mct_sMatPInitnc(sMatp_Fa2o,gsMap_a,gsMap_o,"seq_maps.rc", &
          "atm2ocnFmapname:","atm2ocnFmaptype:",mpicom,&
 	  ni_i=ni_a,nj_i=nj_a,ni_o=ni_o,nj_o=nj_o)

       call shr_mct_sMatPInitnc(sMatp_Sa2o,gsMap_a,gsMap_o,"seq_maps.rc", &
          "atm2ocnSmapname:","atm2ocnSmaptype:",mpicom)

    endif

  end subroutine map_atm2ocn_init_mct

!=======================================================================

  subroutine map_ocn2atm_init_mct( cdata_o, cdata_a)

    !-----------------------------------------------------
    !
    ! Arguments
    !
    type(seq_cdata),intent(in) :: cdata_o
    type(seq_cdata),intent(in) :: cdata_a
    !
    ! Local Variables
    !
    type(seq_infodata_type), pointer :: infodata
    integer                  :: ka, km            ! indices
    type(mct_gsMap), pointer :: gsMap_o           ! ocn gsMap
    type(mct_gsMap), pointer :: gsMap_a           ! atm gsMap
    type(mct_ggrid), pointer :: dom_o             ! ocn domain
    type(mct_ggrid), pointer :: dom_a             ! atm domain
    integer                  :: mpicom            ! communicator spanning atm and ocn
    integer                  :: lsize             ! size of attribute vector
    type(mct_aVect)          :: areasrc           ! ocn areas from mapping file
    type(mct_aVect)          :: areadst           ! atm areas from mapping file
    character(*),parameter :: subName = '(map_ocn2atm_init_mct) '
    !-----------------------------------------------------

    call seq_cdata_setptrs(cdata_o, gsMap=gsMap_o, dom=dom_o)
    call seq_cdata_setptrs(cdata_a, gsMap=gsMap_a, dom=dom_a)
    call seq_cdata_setptrs(cdata_o, mpicom=mpicom, infodata=infodata)

    call seq_infodata_GetData( infodata, samegrid_ao=samegrid_mapa2o)

    ! Initialize ocn->atm mapping or rearranging

    if (samegrid_mapa2o) then

       call mct_rearr_init(gsMap_o, gsMap_a, mpicom, Re_ocn2atm)

    else

       lsize = mct_gsMap_lsize(gsMap_o, mpicom)
       call mct_aVect_init( areasrc, rList="aream", lsize=lsize )
       
       lsize = mct_gsMap_lsize(gsMap_a, mpicom)
       call mct_aVect_init( areadst, rList="aream", lsize=lsize )
       
       call shr_mct_sMatPInitnc(sMatp_Fo2a, gsMap_o, gsMap_a, "seq_maps.rc", &
            "ocn2atmFmapname:", "ocn2atmFmaptype:", mpicom, &
            areasrc=areasrc, areadst=areadst)

       call shr_mct_sMatPInitnc(sMatp_So2a, gsMap_o, gsMap_a, "seq_maps.rc", &
            "ocn2atmSmapname:", "ocn2atmSmaptype:", mpicom)

       !--- copy aream from mapping files

       km = mct_aVect_indexRA(dom_o%data,"aream")
       ka = mct_aVect_indexRA(areasrc   ,"aream")
       dom_o%data%rAttr(km,:) = areasrc%rAttr(ka,:)

       km = mct_aVect_indexRA(dom_a%data,"aream")
       ka = mct_aVect_indexRA(areadst   ,"aream")
       dom_a%data%rAttr(km,:) = areadst%rAttr(ka,:)

       call mct_aVect_clean(areasrc)
       call mct_aVect_clean(areadst)      

    endif

 end subroutine map_ocn2atm_init_mct

!=======================================================================

 subroutine map_atm2ocn_mct( cdata_a, av_a, cdata_o, av_o, &
                             fluxlist, statelist)

    !-----------------------------------------------------
    !
    ! Arguments
    !
    type(seq_cdata) ,intent(in)          :: cdata_a
    type(mct_aVect) ,intent(in)          :: av_a
    type(seq_cdata) ,intent(in)          :: cdata_o
    type(mct_aVect) ,intent(out)         :: av_o
    character(len=*),intent(in),optional :: fluxlist
    character(len=*),intent(in),optional :: statelist
    !
    ! Local Variables
    !
    type(seq_infodata_type), pointer :: infodata
    type(mct_ggrid), pointer :: dom_o
    type(mct_ggrid), pointer :: dom_a
    type(mct_gsMap), pointer :: gsmap_a
    integer                  :: mpicom
    integer                  :: lsize
    integer                  :: ku,kv
    type(mct_aVect)          :: av_o_f     ! temporary flux attribute vector
    type(mct_aVect)          :: av_o_s     ! temporary state attribute vector
    logical                  :: do_npfix          ! npfix on or off
    character(*),parameter :: subName = '(map_atm2ocn_mct) '
    !-----------------------------------------------------

    call seq_cdata_setptrs(cdata_o, infodata=infodata)
    call seq_infodata_GetData( infodata, npfix=do_npfix)

    if (samegrid_mapa2o) then

       if (present(fluxlist) .or. present(statelist)) then
          if (present(fluxlist)) then
	     call mct_rearr_rearrange_fldlist(av_a, av_o, Re_atm2ocn, VECTOR=usevector, &
                  ALLTOALL=usealltoall, fldlist=fluxlist)
          endif
          if (present(statelist)) then
	     call mct_rearr_rearrange_fldlist(av_a, av_o, Re_atm2ocn, VECTOR=usevector, &
                  ALLTOALL=usealltoall, fldlist=statelist)
          endif
       else
          call mct_rearr_rearrange(av_a, av_o, Re_atm2ocn, VECTOR=usevector, ALLTOALL=usealltoall)
       endif

    else
       
      if (present(fluxlist) .or. present(statelist)) then
         if (present(fluxlist)) then
            lsize = mct_aVect_lsize(av_o)
            call mct_aVect_init (av_o_f, rlist=fluxlist , lsize=lsize)
            call mct_sMat_avMult(av_a, sMatp_Fa2o, av_o_f, VECTOR=usevector, rList=fluxlist)
            call mct_aVect_copy (aVin=av_o_f, aVout=av_o, vector=usevector)
            call mct_aVect_clean(av_o_f)
         end if
         if (present(statelist)) then
            lsize = mct_aVect_lsize(av_o)
            call mct_aVect_init (av_o_s, rlist=statelist, lsize=lsize)
            call mct_sMat_avMult(av_a, sMatp_Sa2o, av_o_s, VECTOR=usevector, rList=statelist)
            call mct_aVect_copy (aVin=av_o_s, aVout=av_o, vector=usevector)
            call mct_aVect_clean(av_o_s)
         end if
      else
         !--- default is flux mapping
         call mct_sMat_avMult(av_a, sMatp_Fa2o, av_o, VECTOR=usevector)
      endif
         
      ! Correct a->o vector mapping near NP if appropriate
      
      if (do_npfix) then
      if (present(statelist)) then 
         ku = mct_aVect_indexRA(av_a, 'Sa_u', perrwith='quiet')
         kv = mct_aVect_indexRA(av_a, 'Sa_v', perrwith='quiet')
         if (ku /= 0 .and. kv /= 0) then
            call seq_cdata_setptrs(cdata_o, dom=dom_o)
            call seq_cdata_setptrs(cdata_a, dom=dom_a, gsmap=gsmap_a)
            call seq_cdata_setptrs(cdata_a, mpicom=mpicom)
            call map_npfixNew4R(av_a, av_o, 'Sa_u', 'Sa_v', gsmap_a, dom_a, dom_o, ni_a, nj_a, mpicom)
         end if
      end if
      end if

    endif

  end subroutine map_atm2ocn_mct

!=======================================================================

  subroutine map_ocn2atm_mct( cdata_o, av_o, cdata_a, av_a, &
                              fractions_o, fractions_a, &
	                      fluxlist, statelist )

    !-----------------------------------------------------
    !
    ! Arguments
    !
    type(seq_cdata) , intent(in)           :: cdata_o
    type(mct_aVect) , intent(in)           :: av_o
    type(seq_cdata) , intent(in)           :: cdata_a
    type(mct_aVect) , intent(out)          :: av_a
    type(mct_aVect) , intent(in), optional :: fractions_o
    type(mct_aVect) , intent(in), optional :: fractions_a
    character(len=*), intent(in), optional :: fluxlist
    character(len=*), intent(in), optional :: statelist
    !
    ! Local variables
    !
    type(mct_aVect)          :: temp
    type(mct_aVect)          :: av_a_f, av_a_s
    integer                  :: ka, ko, kSo_t
    integer                  :: numats,i,j,ier
    integer                  :: lsize
    real(R8),allocatable     :: recip(:)
    integer                  :: rcnt
    real(R8)                 :: rmax,rsum,rval
    character(*),parameter :: subName = '(map_ocn2atm_mct) '
    !-----------------------------------------------------

    if (samegrid_mapa2o) then

       if (present(fluxlist) .or. present(statelist)) then
          if (present(fluxlist)) then
             call mct_rearr_rearrange_fldlist(av_o, av_a, Re_ocn2atm, VECTOR=usevector, &
                  ALLTOALL=usealltoall, fldlist=fluxlist)
          end if
          if (present(statelist)) then
             call mct_rearr_rearrange_fldlist(av_o, av_a, Re_ocn2atm, VECTOR=usevector, &
                  ALLTOALL=usealltoall, fldlist=statelist)
          end if
       else
          call mct_rearr_rearrange(av_o, av_a, Re_ocn2atm, VECTOR=usevector, ALLTOALL=usealltoall)
       endif

    else

       ! Normalize input data with fraction of atmosphere on ocean grid
          
       if (present(fractions_o) .and. present(fractions_a)) then

          lsize = mct_aVect_lsize(av_o)
          numats = mct_aVect_nRAttr(av_o)
          ko = mct_aVect_indexRA(fractions_o,"ofrac")
          call mct_aVect_init(temp, av_o, lsize=lsize)
          do j=1,lsize
             do i=1,numats
                temp%rAttr(i,j) = av_o%rAttr(i,j)* fractions_o%rAttr(ko,j)
             end do
          end do
          
          ! Perform remapping
          
          if (present(fluxlist) .or. present(statelist)) then
             if (present (fluxlist)) then 
                lsize = mct_aVect_lsize(av_a)
                call mct_aVect_init (av_a_f, rlist=fluxlist , lsize=lsize)
                call mct_sMat_avMult(temp, sMatp_Fo2a, av_a_f, VECTOR=usevector, rList=fluxlist )
                call mct_aVect_copy (aVin=av_a_f, aVout=av_a, vector=usevector)
                call mct_aVect_clean(av_a_f)
             end if
             if (present(statelist)) then 
                lsize = mct_aVect_lsize(av_a)
                call mct_aVect_init (av_a_s, rlist=statelist, lsize=lsize)
                call mct_sMat_avMult(temp, sMatp_So2a, av_a_s, VECTOR=usevector, rList=statelist)
                call mct_aVect_copy (aVin=av_a_s, aVout=av_a, vector=usevector)
                call mct_aVect_clean(av_a_s)
             end if
          else
             ! --- default is flux mapping
             call mct_sMat_avMult(temp, sMatp_Fo2a, av_a_f, VECTOR=usevector)
          endif

          ! Clean up temporary vector
          call mct_aVect_clean(temp)
          

          ! Denormalize output data with fraction of ocean on atmosphere grid

          lsize = mct_aVect_lsize(av_a)
          numats = mct_aVect_nRAttr(av_a)
          allocate(recip(lsize),stat=ier)
          if(ier/=0) call die(subName,'allocate recip',ier)
          
          ko    = mct_aVect_indexRA(fractions_a, "ofrac")
          do j=1,lsize
             recip(j) = 0.0_R8
             if (fractions_a%rAttr(ko,j) /= 0.0_R8) then
                recip(j)= 1.0_R8 / fractions_a%rAttr(ko,j)
             end if
             do i =1,numats
                   av_a%rAttr(i,j) = av_a%rAttr(i,j) * recip(j)
             end do
          end do

          deallocate(recip,stat=ier)
          if(ier/=0) call die(subName,'deallocate recip',ier)

       else
         
          if (present(fluxlist) .or. present(statelist)) then
             if (present (fluxlist)) then 
                lsize = mct_aVect_lsize(av_a)
                call mct_aVect_init (av_a_f, rlist=fluxlist , lsize=lsize)
                call mct_sMat_avMult(av_o, sMatp_Fo2a, av_a_f, VECTOR=usevector, rList=fluxlist )
                call mct_aVect_copy (aVin=av_a_f, aVout=av_a, vector=usevector)
                call mct_aVect_clean(av_a_f)
             end if
             if (present(statelist)) then 
                lsize = mct_aVect_lsize(av_a)
                call mct_aVect_init (av_a_s, rlist=statelist, lsize=lsize)
                call mct_sMat_avMult(av_o, sMatp_So2a, av_a_s, VECTOR=usevector, rList=statelist)
                call mct_aVect_copy (aVin=av_a_s, aVout=av_a, vector=usevector)
                call mct_aVect_clean(av_a_s)
             end if
          else
             ! --- default is flux mapping
             call mct_sMat_avMult(av_o, sMatp_Fo2a, av_a, VECTOR=usevector)
          endif
          
       end if

   end if


 end subroutine map_ocn2atm_mct

!=======================================================================

 subroutine map_npFixNew4R(buni,buno,fld1,fld2,gsmapi,domi,domo,ni_i,nj_i,mpicom)

   !===============================================================================
   !    Correct the north pole mapping of velocity fields from the atm to ocn
   !    grids.  This assumes the input grid is a regular lat/lon with the north
   !    pole surrounded by the last latitude line in the input array.  The
   !    longitudes in the last latitude must be ordered and equally spaced.
   !
   !    4R is a low memory version of 3R.
   !    This version (New4R) is the same as 3R except it uses a lot less memory
   !    and is a bit faster.  Like 3R, it saves data between calls and so
   !    assumes the input grid remains constant for all calls.  This is bfb
   !    with 3R on bluevista as of 2/15/07.
   !
   !    !REVISION HISTORY:
   !    2007-Feb-12 - T. Craig -- modified New3R to reduce memory
   !    2007-Apr-27 - M. Vertenstein - implemented in sequential system
   !===============================================================================
   
#include <mpif.h>  
   ! 
   ! Arguments
   !
   type(mct_Avect),intent(in)   :: buni    ! input  attribute vec
   type(mct_Avect),intent(out)  :: buno    ! output attribute vec
   character(*)   ,intent(in)   :: fld1    ! name of first input field
   character(*)   ,intent(in)   :: fld2    ! name of second input field
   integer        ,intent(in)   :: ni_i    ! number of longitudes in input grid
   integer        ,intent(in)   :: nj_i    ! number of latitudes in input grid
   type(mct_gsMap),pointer      :: gsmapi  ! input gsmap
   type(mct_gGrid),pointer      :: domi    ! input domain 
   type(mct_gGrid),pointer      :: domo    ! output domain 
   integer        ,intent(in)   :: mpicom  ! mpi communicator group 
   !
   ! Local Variables
   !
   integer(IN)  :: n,m                           ! generic indices
   integer(IN)  :: n1,n2,n3                      ! generic indices
   integer(IN)  :: kui,kvi                       ! field indices
   integer(IN)  :: kuo,kvo                       ! field indices
   integer(IN)  :: kin                           ! index index
   integer(IN)  :: nmin,nmax                     ! indices of highest latitude in input
   integer(IN)  :: npts                          ! local number of points in an aV
   integer(IN)  :: num                           ! number of points at highest latitude
   integer(IN)  :: kloni                         ! longitude index on input domain
   integer(IN)  :: klati                         ! latitude index on input domain
   integer(IN)  :: klono                         ! longitude index on output domain
   integer(IN)  :: klato                         ! latitude index on output domain
   integer(IN)  :: index                         ! index value
   real(R8)     :: rindex                        ! index value
   real(R8)     :: latmax                        ! value of highest latitude
   real(R8)     :: olon,olat                     ! output bundle lon/lat
   real(R8)     :: ilon,ilat                     ! input bundle lon/lat
   real(R8)     :: npu,npv                       ! np velocity fields relative to lon
   real(R8)     :: theta1,theta2                 ! angles for trig functions
   real(R8),allocatable,save :: ilon1(:)         ! lon of input grid at highest latitude
   real(R8),allocatable,save :: ilat1(:)         ! lat of input grid at highest latitude
   real(R8),allocatable,save :: ilon2(:)         ! lon of input grid at highest latitude
   real(R8),allocatable,save :: ilat2(:)         ! lat of input grid at highest latitude
   real(R8),allocatable      :: rarray(:)        ! temporary array
   real(R8),allocatable      :: rarray2(:,:)     ! temporary array
   real(R8)     :: w1,w2,w3,w4                   ! weights
   real(R8)     :: f1,f2,f3,f4                   ! function values
   real(R8)     :: alpha,beta                    ! used to generate weights
   real(R8)     :: rtmp                          ! real temporary
   real(R8),allocatable :: lData(:,:)            ! last lat local input bundle data
                                                 ! also compressed global data
   real(R8)   ,allocatable,save :: alphafound(:) ! list of found alphas
   real(R8)   ,allocatable,save :: betafound(:)  ! list of found betas
   integer(IN),allocatable,save :: mfound(:)     ! list of found ms
   integer(IN),allocatable,save :: nfound(:)     ! list of found ns
   real(R8)   ,allocatable      :: rfound(:)     ! temporary for copy
   integer(IN),allocatable      :: ifound(:)     ! temporary for copy
   integer(IN),save             :: cntfound      ! number found
   integer(IN),save             :: cntf_tot      ! cntfound total
   integer(IN)  :: cnt                           ! loop counter
   logical      :: found                         ! search for new interpolation
   integer(IN)  :: rcode                         ! error code
   integer(IN)  :: np1                           ! n+1 or tmp
   real(R8)     :: ilon1x                        ! tmp
   integer(IN),pointer,save  :: gindex(:)        ! global index 
   logical,save :: first_call = .true.           ! flags 1st invocation of routine
   integer(IN),save   :: tnpf1,tnpf2,tnpf3,tnpf4,tnpf5,tnpf6,tnpf7,tnpf8,tnpf9
   real(R8),parameter :: shr_const_deg2rad = shr_const_pi/180.0_R8  ! deg to rads
   integer(IN) :: mype

   !--- formats ---
   character(*),parameter :: subName = '(map_npFixNew4R) '
   character(*),parameter :: F00 = "('(map_npFixNew4R) ',8a)"
   character(*),parameter :: F01 = "('(map_npFixNew4R) ',a,i12)"
   !-------------------------------------------------------------------------------

   call MPI_COMM_RANK(mpicom,mype,rcode)

   kui   = mct_aVect_indexRA(buni,fld1,perrWith=subName)
   kvi   = mct_aVect_indexRA(buni,fld2,perrWith=subName)
   kuo   = mct_aVect_indexRA(buno,fld1,perrWith=subName)
   kvo   = mct_aVect_indexRA(buno,fld2,perrWith=subName)

!  tcx 3/19/08, don't use GlobGridNum, it's not set properly in models
!   kin   = mct_aVect_indexIA(domi%data,"GlobGridNum",perrWith=subName)
   klati = mct_aVect_indexRA(domi%data,"lat"        ,perrWith=subName)
   kloni = mct_aVect_indexRA(domi%data,"lon"        ,perrWith=subName)
   klato = mct_aVect_indexRA(domo%data,"lat"        ,perrWith=subName)
   klono = mct_aVect_indexRA(domo%data,"lon"        ,perrWith=subName)

   ! ni_i, nj_i should be read in from mapping file

   nmin = (ni_i)*(nj_i-1) + 1
   nmax = ni_i*nj_i
   num  = ni_i

   !---------------------------------------------------------------------------
   ! Initialization only
   !---------------------------------------------------------------------------

   if (first_call) then

     if (loglevel > 0) write(logunit,F00) " compute bilinear weights & indicies for NP region."

     allocate(ilon1(num))
     allocate(ilon2(num))
     allocate(ilat1(num))
     allocate(ilat2(num))

     ilon1 = 0._r8
     ilon2 = 0._r8
     ilat1 = 0._r8
     ilat2 = 0._r8

     npts = mct_aVect_lSize(domi%data)
     call mct_gsMap_orderedPoints(gsMapi, mype, gindex)
     do m=1,npts
       if (gindex(m).ge.nmin) then               ! are on highest latitude
         n = gindex(m) - nmin + 1                ! n is 1->ni_i lon index on highest latitude
         rtmp = domi%data%rAttr(kloni,m)      ! rtmp is longitude value on highest latitude
         ilon1(n) = mod(rtmp+360._R8,360._R8) ! ilon1(n) is longitude val mapped from 0->360
         ilat1(n) = domi%data%rAttr(klati,m)  ! ilat1(n) values should all be the same (i.e. highest lat)
       endif
     enddo

     !--- all gather local data, MPI_SUM is low memory and simple
     !--- but is a performance penalty compared to gatherv and copy
     !--- or a fancy send/recv 

     allocate(rarray(num))
     rarray = ilat1
     call MPI_ALLREDUCE(rarray,ilat1,num,MPI_REAL8,MPI_SUM,mpicom,rcode)
     if (rcode.ne.0) then
       write(logunit,*) trim(subName),' ilat1 rcode error ',rcode
       call shr_sys_abort()
     endif
     rarray = ilon1
     call MPI_ALLREDUCE(rarray,ilon1,num,MPI_REAL8,MPI_SUM,mpicom,rcode)
     if (rcode.ne.0) then
       write(logunit,*) trim(subName),' ilon1 rcode error ',rcode
       call shr_sys_abort()
     endif

     do n = 1,num
       np1 = mod(n,num)+1
       ilat2(n) = ilat1(np1)
       ilon2(n) = ilon1(np1)
       if (ilon2(n) < ilon1(n)) ilon2(n) = ilon2(n) + 360._R8
     enddo

     do n = 1,num
       if (ilat1(n) /= ilat2(n)) then
          write(logunit,*) trim(subname),' ERROR: ilat1 ne ilat2 ',n,ilat1(n),ilat2(n)
          call shr_sys_abort(trim(subname)//' ERROR: ilat1 ne ilat2')
       endif
       if (ilon2(n) < ilon1(n)) then
          write(logunit,*) trim(subname),' ERROR: ilon2 lt ilon1 ',n,ilon1(n),ilon2(n)
          call shr_sys_abort(trim(subname)//' ERROR: ilon2 ilon1 error')
       endif
       ! tcraig check that lon diffs are reasonable 4x average dlon seems like reasonable limit
       if (ilon2(n) - ilon1(n) > (360.0_R8/(num*1.0_R8))*4.0) then
          write(logunit,*) trim(subname),' ERROR: ilon1,ilon2 ',n,ilon1(n),ilon2(n)
          call shr_sys_abort(trim(subname)//' ERROR: ilon2 ilon1 size diff')
       endif
     enddo

     latmax = maxval(ilat1)

     !--- compute weights and save them ---

     npts = mct_aVect_lSize(buno)
     allocate(mfound(npts),nfound(npts),alphafound(npts),betafound(npts))
     cntfound = 0
     do m = 1,npts
       olat = domo%data%rAttr(klato,m)
       if (olat >= latmax) then
         rtmp = domo%data%rAttr(klono,m)
         olon = mod(rtmp,360._R8)
         n = 1
         found = .false.
         do while (n <= num .and. .not.found )
           if (    olon         >= ilon1(n) .and. olon < ilon2(n) .or.   &
                   olon+360._R8 >= ilon1(n) .and. olon < ilon2(n)) then


!tcx ilat2==ilat1 so don't average
!--->        ilat = (ilat1(n) + ilat2(n)) * 0.5_R8
             ilat = ilat1(n)
             if (ilon2(n) == ilon1(n)) then
               alpha = 0.5_R8
             else if (    olon >= ilon1(n) .and. olon < ilon2(n)) then
               alpha = (olon - ilon1(n)) / (ilon2(n) - ilon1(n))
             else if (olon+360._R8>= ilon1(n) .and. olon < ilon2(n)) then
               alpha = (olon+360._R8 - ilon1(n)) / (ilon2(n) - ilon1(n))
             else
               write(logunit,*) subName,' ERROR: olon ',olon,ilon1(n),ilon2(n)
             endif
             if (ilat >= 90._R8) then
               beta  = 1.0_R8
             else
               beta  = (olat - ilat) / (90._R8 - ilat)
             endif

             cntfound = cntfound + 1
             mfound(cntfound) = m
             nfound(cntfound) = n
             alphafound(cntfound) = alpha
             betafound(cntfound) = beta
             found = .true.

           endif
           n = n + 1     ! normal increment
         enddo
         if ( .not.found ) then
           write(logunit,*) subName,' ERROR: found = false ',found,m,olon,olat
         endif
       endif
     end do

     allocate(ifound(npts))
     ifound(1:cntfound) = mfound(1:cntfound)
     deallocate(mfound)
     if (cntfound > 0) then
        allocate(mfound(cntfound))
        mfound(1:cntfound) = ifound(1:cntfound)
     endif

     ifound(1:cntfound) = nfound(1:cntfound)
     deallocate(nfound)
     if (cntfound > 0) then
        allocate(nfound(cntfound))
        nfound(1:cntfound) = ifound(1:cntfound)
     endif
     deallocate(ifound)

     allocate(rfound(npts))
     rfound(1:cntfound) = alphafound(1:cntfound)
     deallocate(alphafound)
     if (cntfound > 0) then
        allocate(alphafound(cntfound))
        alphafound(1:cntfound) = rfound(1:cntfound)
     endif

     rfound(1:cntfound) = betafound(1:cntfound)
     deallocate(betafound)
     if (cntfound > 0) then
        allocate(betafound(cntfound))
        betafound(1:cntfound) = rfound(1:cntfound)
     endif
     deallocate(rfound)

     call MPI_ALLREDUCE(cntfound,cntf_tot,1,MPI_INTEGER,MPI_SUM,mpicom,rcode)
     if (mype == 0) then
        write(logunit,F01) ' total npfix points found = ',cntf_tot
     endif

     first_call = .false.

   endif

   !---------------------------------------------------------------------------
   ! Return if there is nothing to do; must be no points on any pes
   ! If there are any npfix points, all pes must continue to the np u,v calc
   !---------------------------------------------------------------------------

   if (cntf_tot < 1) then
      return
   endif

   !---------------------------------------------------------------------------
   ! Non-initialization, run-time fix
   !---------------------------------------------------------------------------

   !--- barrier not required but interesting for timing. ---
   !  call shr_mpi_barrier(mpicom,subName//" barrier")

   !--- extract index, u, v from buni ---

   allocate(lData(3,num))
   lData = 0._R8
   npts = mct_aVect_lSize(buni)
   do n=1,npts
     if (gindex(n).ge.nmin) then
       m = gindex(n) - nmin + 1
       lData(1,m) = gindex(n)
       lData(2,m) = buni%rAttr(kui,n)
       lData(3,m) = buni%rAttr(kvi,n)
     endif
   enddo

   !--- all gather local data, MPI_SUM is low memory and simple
   !--- but is a performance penalty compared to gatherv and copy
   !--- KLUDGE - this should be looked at when it becomes a performance/memory
   !--- penalty   

   allocate(rarray2(3,num))
   rarray2=lData
   call MPI_ALLREDUCE(rarray2,lData,3*num,MPI_REAL8,MPI_SUM,mpicom,rcode)
   deallocate(rarray2)

   if (rcode.ne.0) then
     write(logunit,*) trim(subName),' rcode error ',rcode
     call shr_sys_abort()
   endif

   do n2=1,num
     if (lData(1,n2).lt.0.1_R8) then
       write(logunit,*) trim(subName),' error allreduce ',n2
     endif
   enddo

   !--- compute npu, npv (pole data) and initialize ilon,ilat arrays ---

   npu = 0._R8
   npv = 0._R8
   do n = 1,num
     theta1 = ilon1(n)*shr_const_deg2rad
     npu = npu + cos(theta1)*lData(2,n) &
               - sin(theta1)*lData(3,n)
     npv = npv + sin(theta1)*lData(2,n) &
               + cos(theta1)*lData(3,n)
   enddo
   npu = npu / real(num,R8)
   npv = npv / real(num,R8)

   !--- compute updated pole vectors ---

!DIR$ CONCURRENT
   do cnt = 1,cntfound
      m     = mfound(cnt)
      n     = nfound(cnt)
      np1   = mod(n,num)+1
      alpha = alphafound(cnt)
      beta  = betafound(cnt)

      w1 = (1.0_R8-alpha)*(1.0_R8-beta)
      w2 = (    alpha)*(1.0_R8-beta)
      w3 = (    alpha)*(    beta)
      w4 = (1.0_R8-alpha)*(    beta)

      theta1 = ilon1(n)*shr_const_deg2rad
      theta2 = ilon2(n)*shr_const_deg2rad

      f1 = lData(2,n)
      f2 = lData(2,np1)
      f3 =  cos(theta1)*npu + sin(theta1)*npv
      f4 =  cos(theta2)*npu + sin(theta2)*npv
      rtmp = w1*f1 + w2*f2 + w3*f3 + w4*f4
      buno%rAttr(kuo,m) = w1*f1 + w2*f2 + w3*f3 + w4*f4

      f1 = lData(3,n)
      f2 = lData(3,np1)
      f3 = -sin(theta1)*npu + cos(theta1)*npv
      f4 = -sin(theta2)*npu + cos(theta2)*npv
      rtmp = w1*f1 + w2*f2 + w3*f3 + w4*f4
      buno%rAttr(kvo,m) = w1*f1 + w2*f2 + w3*f3 + w4*f4
   enddo

   deallocate(lData)
   first_call = .false.

end subroutine map_npFixNew4R

!===============================================================================

end module map_atmocn_mct
