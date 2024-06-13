module seq_rearr_mod

!---------------------------------------------------------------------
!
! Purpose:
!
! Shared routines for computation of rearrangers
!       
! Author: T Craig
!
!---------------------------------------------------------------------

  use shr_sys_mod
  use mct_mod
  use seq_flds_mod
  use seq_cdata_mod
  use seq_comm_mct
  use seq_diag_mct
  use shr_kind_mod, only: R8 => SHR_KIND_R8

  implicit none
  save
  private  ! except

!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------

  public  :: seq_rearr_init
  public  :: seq_rearr_gsmapIdentical
  private :: seq_rearr_gsmapExtend
  private :: seq_rearr_gsmapCreate
  private :: seq_rearr_avExtend
  private :: seq_rearr_avCreate

!--------------------------------------------------------------------------
! Public data
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
! Private data
!--------------------------------------------------------------------------

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

  character(*),parameter :: subName = '(seq_rearr_mct)'
  real(r8),parameter :: c1 = 1.0_r8

#include <mpif.h>

!=======================================================================
   contains
!=======================================================================

  subroutine seq_rearr_init( cdata_old, AV1_old, AV2_old, ID_old, &
                             cdata_new, AV1_new, AV2_new, ID_new, ID_join, &
                             Re_old2new, Re_new2old)

    ! This routine modifies and initializes a bunch of datatypes to
    ! handle rearranging of AVs between different pes for a single
    ! component.  Basically, the "old" data exists on a set of pes
    ! associated with ID_old.  We want to initialize a rearranger
    ! between ID_old and ID_new and datatypes on ID_new that will
    ! supporting moving of data between ID_old and ID_new.  Nothing
    ! is set on ID_new except it's ID, so we need to create a gsmap
    ! on those pes associated with the same global grid (but maybe
    ! with a different decomp).  Then we can initialize a rearranger
    ! between the old and new gsmap as well as initialize the new
    ! AVs and new domain.  In the process, AVs and gsmaps need to be
    ! "extended" to span all pes in the joined ID, instead of just
    ! the pes of their original ID.  This has nothing to do with
    ! changing the gsmap, decomp, or data.  It just involves copying
    ! the gsmap to the new joined pes and allocating AVs on those
    ! extended pes with local size zero and the proper fields.  These
    ! extended versions of gsmaps and AVs are then carried around
    ! in the driver without any problems, but they can be used to
    ! rearrange data across the joined pes.

    implicit none

    !-----------------------------------------------------
    ! 
    ! Arguments
    !
    type(seq_cdata),intent(in)    :: cdata_old
    type(mct_aVect),intent(inout),optional :: AV1_old
    type(mct_aVect),intent(inout),optional :: AV2_old
    integer        ,intent(in)    :: ID_old
    type(seq_cdata),intent(inout) :: cdata_new
    type(mct_aVect),intent(inout),optional :: AV1_new
    type(mct_aVect),intent(inout),optional :: AV2_new
    integer        ,intent(in)    :: ID_new
    integer        ,intent(in)    :: ID_join
    type(mct_rearr),intent(inout) :: Re_old2new
    type(mct_rearr),intent(inout) :: Re_new2old
    !
    ! Local Variables
    !
    character(len=*),parameter :: subname = "(seq_rearr_init) "
    type(mct_gsMap), pointer :: gsmap_new           ! new gsMap
    type(mct_gsMap), pointer :: gsmap_old           ! old gsMap
    type(mct_gGrid), pointer :: dom_new             ! new domain
    type(mct_gGrid), pointer :: dom_old             ! old domain

    integer :: lsize
    integer :: mpicom_old
    integer :: mpicom_new
    integer :: mpicom_join
    logical :: master_join
    integer :: mpigrp_old, mpigrp_join
    integer :: ierr
    integer :: ID

    type(mct_gsMap)          :: gsmap_new_join     ! gsmap_new on joined id
    type(mct_gsMap)          :: gsmap_old_join     ! gsmap_old on joined id
    type(mct_gGrid)          :: dom_join          ! gathered grid
    type(mct_gGrid)          :: dom_empty         ! empty grid
    !-----------------------------------------------------

    ! --- Setup data for use and make sure the old ID is ok

    call seq_cdata_setptrs(cdata_old, gsMap=gsmap_old, dom=dom_old, ID=ID)
    call seq_cdata_setptrs(cdata_new, gsMap=gsmap_new, dom=dom_new)

    if (ID /= ID_old) then
       write(logunit,*) subname,' ERROR: inconsistent ID',ID,ID_old
       call shr_sys_abort()
    endif

    call seq_comm_setptrs(ID_old  ,mpicom=mpicom_old  ,mpigrp=mpigrp_old  )
    call seq_comm_setptrs(ID_new  ,mpicom=mpicom_new                     )
    call seq_comm_setptrs(ID_join,mpicom=mpicom_join,mpigrp=mpigrp_join,iamroot=master_join)

    ! --- Set gsmaps
    ! ---   Extend the old one to now span all pes on ID_join
    ! ---   Create a new gsmap on pes associated with ID_new using info from the old one
    ! ---   Extend the new one to span all pes on ID_join

    call seq_rearr_gsmapExtend(gsmap_old     , mpicom_old  , gsmap_old_join, mpicom_join, ID_join)
    call seq_rearr_gsmapCreate(gsmap_old_join, mpicom_join , gsmap_new     , mpicom_new , ID_new  )
    call seq_rearr_gsmapExtend(gsmap_new     , mpicom_new  , gsmap_new_join, mpicom_join, ID_join)

    ! --- Initialize rearrangers based on join gsmaps

    call mct_rearr_init(gsmap_old_join, gsmap_new_join, mpicom_join, Re_old2new)
    call mct_rearr_init(gsmap_new_join, gsmap_old_join, mpicom_join, Re_new2old)

    ! --- Finish setting up the new information.  extend dom_old and rearrange it to dom_new

    lsize = mct_gsMap_lsize(gsmap_new_join, mpicom_join)
    call mct_gGrid_init( GGrid=dom_new, CoordChars=seq_flds_dom_coord, OtherChars=seq_flds_dom_other, lsize=lsize )
    call mct_avect_zero(dom_new%data)
    call seq_rearr_avExtend(dom_old%data, ID_old, ID_join)
    call mct_rearr_rearrange(dom_old%data, dom_new%data, Re_old2new, VECTOR=usevector, ALLTOALL=usealltoall)

    ! --- Clean temporary gsmaps

    call mct_gsMap_clean(gsmap_old_join)
    call mct_gsMap_clean(gsmap_new_join)

   ! --- Extend old avs and initialize new avs for use in the future

   lsize = 0
   if (seq_comm_iamin(ID_new)) then
      lsize = mct_gsMap_lsize(gsMap_new,mpicom_new)
   endif
   if (present(AV1_old)) then
      call seq_rearr_avExtend(AV1_old, ID_old, ID_join)
      if (present(AV1_new)) then
         call seq_rearr_avCreate(AV1_old, ID_old, AV1_new, ID_join, lsize)
      endif
   endif
   if (present(AV2_old)) then
      call seq_rearr_avExtend(AV2_old, ID_old, ID_join)
      if (present(AV2_new)) then
         call seq_rearr_avCreate(AV2_old, ID_old, AV2_new, ID_join, lsize)
      endif
   endif

  end subroutine seq_rearr_init

!=======================================================================

  subroutine seq_rearr_gsmapExtend(gsmapi, mpicomi, gsmapo, mpicomo, compido)

    !----------------------------------------------------------------
    ! Extend/Convert a gsmap from one mpicom to another mpicom that contains
    ! at least all the pes that gsmap uses, but with different ranks
    !----------------------------------------------------------------
  
    implicit none
    type(mct_gsMap), intent(IN) :: gsmapi
    integer        , intent(IN) :: mpicomi
    type(mct_gsMap), intent(OUT):: gsmapo
    integer        , intent(IN) :: mpicomo
    integer        , intent(IN) :: compido

    character(len=*),parameter :: subname = "(seq_rearr_gsmapExtend) "
    integer :: n
    integer :: ngseg
    integer :: gsize
    integer :: msizei,msizeo
    integer :: mrank,mranko,mrankog   ! sets pe rank of root mpicomi pe in mpicomo
    integer :: mpigrpi,mpigrpo
    integer :: ierr
    integer, pointer :: pei(:),peo(:)
    integer, pointer :: start(:),length(:),peloc(:)

    mranko = -1

    ! --- create the new gsmap on the mpicomi root only

    if (mpicomi /= MPI_COMM_NULL) then
       call mpi_comm_rank(mpicomi,mrank,ierr)
       call shr_mpi_chkerr(ierr,subname//' gsm_cop mpi_comm_rank i')
       if (mrank == 0) then
          call mpi_comm_group(mpicomi,mpigrpi,ierr)
          call shr_mpi_chkerr(ierr,subname//' gsm_cop mpi_comm_group i')
          call mpi_comm_group(mpicomo,mpigrpo,ierr)
          call shr_mpi_chkerr(ierr,subname//' gsm_cop mpi_comm_group o')
          call mpi_comm_size(mpicomi,msizei,ierr)
          call shr_mpi_chkerr(ierr,subname//' gsm_cop mpi_comm_size i')
          call mpi_comm_size(mpicomo,msizeo,ierr)
          call shr_mpi_chkerr(ierr,subname//' gsm_cop mpi_comm_size o')

          ! --- setup the translation of pe numbers from the old gsmap(mpicom)
          ! --- to the new one, pei -> peo

          allocate(pei(0:msizei-1),peo(0:msizei-1))
          do n = 0,msizei-1
             pei(n) = n
          enddo

          peo = -1
          call mpi_group_translate_ranks(mpigrpi,msizei,pei,mpigrpo,peo,ierr)
          call shr_mpi_chkerr(ierr,subname//' gsm_cop mpi_group_translate_ranks')

          do n = 0,msizei-1
             if (peo(n) < 0 .or. peo(n) > msizeo-1) then
                write(logunit,*) subname,' peo out of bounds ',peo(n),msizeo
                call shr_sys_abort()
             endif
          enddo

          mranko = peo(0)

          ! --- compute the new gsmap which has the same start and length values
          ! --- but peloc is now the mapping of pei to peo

          ngseg = gsmapi%ngseg
          gsize = gsmapi%gsize
          allocate(start(ngseg),length(ngseg),peloc(ngseg))
          do n = 1,ngseg
             start(n)  = gsmapi%start(n)
             length(n) = gsmapi%length(n)
             peloc(n)  = peo(gsmapi%pe_loc(n))
          enddo

          ! --- initialize the gsmap on the root pe

          call mct_gsmap_init(gsmapo,compido,ngseg,gsize,start,length,peloc)

          deallocate(pei,peo,start,length,peloc)
       endif
    endif

    ! --- broadcast via allreduce the mpicomi root pe in mpicomo space
    ! --- mranko is -1 except on the root pe where is it peo of that pe

    call mpi_allreduce(mranko,mrankog,1,MPI_INTEGER,MPI_MAX,mpicomo,ierr)
    call shr_mpi_chkerr(ierr,subname//' gsm_cop mpi_allreduce max')

    ! --- broadcast the gsmap to all pes in mpicomo from mrankog

    call mct_gsmap_bcast(gsmapo, mrankog, mpicomo)

! tcx summarize decomp info
#if (1 == 0)
    write(logunit,*) trim(subname),'tcxa ',mpicomi,mpicomo
    call shr_sys_flush(logunit)
    call mpi_barrier(mpicomo,ierr)

    if (mpicomi /= MPI_COMM_NULL) then
       call mpi_comm_rank(mpicomi,mrank,ierr)
       write(logunit,*) 'tcxbi ',mrank
       if (mrank == 0) then
          write(logunit,*) 'tcxci ',gsmapi%ngseg,size(gsmapi%start),gsmapi%gsize,gsmapi%comp_id
          do n = 1,gsmapi%ngseg
             write(logunit,*) 'tcx gsmti ',n,gsmapi%start(n),gsmapi%length(n),gsmapi%pe_loc(n)
          enddo
          call shr_sys_flush(logunit)
      endif
    endif

    if (mpicomo /= MPI_COMM_NULL) then
       call mpi_comm_rank(mpicomo,mrank,ierr)
       write(logunit,*) 'tcxbo ',mrank
       if (mrank == 0) then
          write(logunit,*) 'tcxco ',gsmapo%ngseg,size(gsmapo%start),gsmapo%gsize,gsmapo%comp_id
          do n = 1,gsmapo%ngseg
             write(logunit,*) 'tcx gsmto ',n,gsmapo%start(n),gsmapo%length(n),gsmapo%pe_loc(n)
          enddo
          call shr_sys_flush(logunit)
       endif
    endif

    call shr_sys_flush(logunit)
    call mpi_barrier(mpicomo,ierr)
#endif


  end subroutine seq_rearr_gsmapExtend
!=======================================================================

  subroutine seq_rearr_gsmapCreate(gsmapi, mpicomi, gsmapo, mpicomo, compido)

    !---------------------------------------------------------------------
    ! creates a new gsmap on a subset of pes, requires setting a new decomp
    !---------------------------------------------------------------------
  
    implicit none
    type(mct_gsMap), intent(IN) :: gsmapi
    integer        , intent(IN) :: mpicomi
    type(mct_gsMap), intent(OUT):: gsmapo
    integer        , intent(IN) :: mpicomo
    integer        , intent(IN) :: compido

    character(len=*),parameter :: subname = "(seq_rearr_gsmapCreate) "
    integer :: n,m,k
    integer :: ktot            ! number of active cells in gsmap
    integer :: apesi, apeso    ! number of active pes in gsmap
    integer ::        lsizeo   ! local size for lindex
    integer :: ngsegi,ngsego   ! ngseg of mpicomi, mpicomo
    integer :: gsizei,gsizeo   ! gsize of mpicomi, mpicomo
    integer :: msizei,msizeo   ! size of mpicomi, mpicomo
    integer :: mranki,mranko   ! rank in mpicomi, mpicomo
    integer :: ierr
    integer :: decomp_type
    integer, pointer :: start(:),length(:),peloc(:),perm(:),gindex(:),lindex(:)
    real(r8):: rpeloc

    ! --- create a new gsmap on new pes based on the old gsmap
    ! --- gsmapi must be known on all mpicomo pes, compute the same 
    ! --- thing on all pes in parallel

    if (mpicomo /= MPI_COMM_NULL) then
       call mpi_comm_rank(mpicomi,mranki,ierr)
       call shr_mpi_chkerr(ierr,subname//' mpi_comm_rank i')
       call mpi_comm_size(mpicomi,msizei,ierr)
       call shr_mpi_chkerr(ierr,subname//' mpi_comm_size i')
       call mpi_comm_rank(mpicomo,mranko,ierr)
       call shr_mpi_chkerr(ierr,subname//' mpi_comm_rank o')
       call mpi_comm_size(mpicomo,msizeo,ierr)
       call shr_mpi_chkerr(ierr,subname//' mpi_comm_size o')

       ngsegi = gsmapi%ngseg
       gsizei = gsmapi%gsize
       gsizeo = gsizei
       call mct_gsMap_activepes(gsmapi,apesi)

       decomp_type = 0

       if (msizeo == apesi) then      ! preserve segments and decomp
!          decomp_type = 1     ! better in cpl to have all decomps "same-ish"
          decomp_type = 2
       elseif (ngsegi >= msizeo) then ! preserve segments, new decomp
          decomp_type = 2
       else                           ! new segments
          decomp_type = 3
       endif

!tcx       decomp_type = 3 ! over ride setting above for testing
!       if (mranko == 0) write(logunit,'(2A,4I)') trim(subname),' decomp_type =',decomp_type,ngsegi,msizeo,apesi

       select case (decomp_type)

       case(1)   ! --- preserve segments and decomp ---------------------

          ! -- copy the gsmap and translate the pes
          call mct_gsMap_copy(gsmapi,gsmapo)
          ngsego = ngsegi
          do n = 1,ngsego
             gsmapo%pe_loc(n) = mod(gsmapo%pe_loc(n),msizeo)    ! translate pes 1:1 from old to new
          enddo

       case(2)   ! --- preserve segments, new decomp --------------------

          ! --- preserve segments, sort the start and length, assign a new pe list
          ngsego = ngsegi
          allocate(start(ngsego),length(ngsego),peloc(ngsego),perm(ngsego))
          do n = 1,ngsego
             start(n)  = gsmapi%start(n)
             length(n) = gsmapi%length(n)
          enddo
          ! --- sort gsmap to minimize permute cost in mct
          call mct_indexset(perm)
          call mct_indexsort(ngsego,perm,start)
          call mct_permute(start,perm,ngsego)
          call mct_permute(length,perm,ngsego)
          ! --- give each pe "equal" number of segments, use reals to avoid integer overflow
          do n = 1,ngsego
             rpeloc = (((msizeo*c1)*((n-1)*c1))/(ngsego*c1))      ! give each pe "equal" number of segments, use reals to avoid integer overflow
             peloc(n) = int(rpeloc)
          enddo
          call mct_gsmap_init(gsmapo,ngsego,start,length,peloc,0,mpicomo,compido,gsizeo)
          deallocate(start,length,peloc,perm)

       case(3)   ! --- new segments, new decomp -------------------------

          ! --- new segments, compute gindex, then parse the gridcells out evenly

          k = 0
          do n = 1,ngsegi
          do m = 1,gsmapi%length(n)
             k = k + 1
             if (k > gsizei) then
                write(logunit,*) trim(subname),' ERROR in gindex ',k,gsizei
                call shr_sys_abort()
             endif
          enddo
          enddo
          ktot = k

          allocate(gindex(ktot),perm(ktot))  

          k = 0
          do n = 1,ngsegi
          do m = 1,gsmapi%length(n)
             k = k + 1
             gindex(k) = gsmapi%start(n) + m - 1
          enddo
          enddo
          call mct_indexset(perm)
          call mct_indexsort(ktot,perm,gindex)
          call mct_permute(gindex,perm,ktot)

          k = 0
          do m = 0,msizeo-1
             lsizeo = ktot/msizeo
             if (m < (ktot - lsizeo*msizeo)) lsizeo = lsizeo + 1
             if (mranko == m) then
                allocate(lindex(lsizeo))
                if (k+lsizeo > ktot) then
                   write(logunit,*) trim(subname),' ERROR: decomp out of bounds ',mranko,k,lsizeo,ktot
                   call shr_sys_abort()
                endif
                lindex(1:lsizeo) = gindex(k+1:k+lsizeo)
!                write(logunit,*) trim(subname),' decomp is ',mranko,lsizeo,k+1,k+lsizeo
             endif
             k = k + lsizeo
          enddo
          if (k /= ktot) then
             write(logunit,*) trim(subname),' ERROR: decomp incomplete ',k,ktot
             call shr_sys_abort()
          endif

          call mct_gsmap_init(gsmapo,lindex,mpicomo,compido,size(lindex),gsizeo)
          deallocate(gindex,perm,lindex)

       end select

       if (mranko == 0) then
          write(logunit,102) trim(subname),' created new gsmap decomp_type =',decomp_type
          write(logunit,102) trim(subname),'   ngseg/gsize        = ', &
             mct_gsmap_ngseg(gsmapo),mct_gsmap_gsize(gsmapo)
          call mct_gsmap_activepes(gsmapo,apeso)
          write(logunit,102) trim(subname),'   mpisize/active_pes = ', &
             msizeo,apeso
          write(logunit,102) trim(subname),'   avg seg per pe/ape = ', &
             mct_gsmap_ngseg(gsmapo)/msizeo,mct_gsmap_ngseg(gsmapo)/apeso
          write(logunit,102) trim(subname),'   nlseg/maxnlsegs    = ', &
             mct_gsmap_nlseg(gsmapo,0),mct_gsmap_maxnlseg(gsmapo)
 102      format(2A,2I8)
       endif

!p       if (.not. mct_gsmap_increasing(gsmapo) ) then
!          write(logunit,*) trim(subname),' ERROR: gsmapo not increasing'
!          call shr_sys_abort()
!       endif

    endif

  end subroutine seq_rearr_gsmapCreate
    
!=======================================================================

subroutine seq_rearr_avExtend(AVin,IDin,ID)

  !-----------------------------------------------------------------------
  ! Extend an AV to a larger set of pes or
  ! Initialize an AV on another set of pes
  !-----------------------------------------------------------------------

  implicit none
  type(mct_aVect), intent(INOUT):: AVin
  integer         ,intent(IN)   :: IDin ! ID associated with AVin
  integer        , intent(IN)   :: ID   ! ID to initialize over

  ! Local variables

  character(len=*),parameter :: subname = "(seq_rearr_avExtend) "
  integer :: mpicom
  integer :: rank,rank2
  integer :: lsizei, lsizen
  integer :: srank,srankg
  integer :: ierr
  integer :: nints
  character(len=1024) :: iList,rList


  call seq_comm_setptrs(ID,mpicom=mpicom,iam=rank)

  ! --- lsizen is the size of the newly initialized AV, zero is valid
  ! --- lsizei is -1 on any peszero on any pes where AV is not yet initialized

  lsizei = -1  
  if (seq_comm_iamin(IDin)) lsizei = mct_aVect_lsize(AVin)
  lsizen = 0

  ! --- find a pe that already has AVin allocated, use MPI_MAX to do so
  ! --- set the pe and broadcast it to all other pes using mpi_allreduce

  srank = -1
  srankg = -1
  if (lsizei > 0) srank = rank

  call mpi_allreduce(srank,srankg,1,MPI_INTEGER,MPI_MAX,mpicom,ierr)
  call shr_mpi_chkerr(ierr,subname//' mpi_allreduce max')

  if (srankg < 0) then
    write(logunit,*) subname,' WARNING AVin empty '
    return
  endif

  ! --- set the iList and rList from the broadcast pe (srankg) and 
  ! --- broadcast the lists

  iList = " "
  rList = " "
  if (rank == srankg) then
    if (mct_aVect_nIAttr(AVin) /= 0) iList = mct_aVect_ExportIList2c(AVin)
    if (mct_aVect_nRattr(AVin) /= 0) rList = mct_aVect_ExportRList2c(AVin)
  endif

  call mpi_bcast(iList,len(iList),MPI_CHARACTER,srankg,mpicom,ierr)
  call mpi_bcast(rList,len(rList),MPI_CHARACTER,srankg,mpicom,ierr)

  ! --- now allocate the AV on any pes where the orig size is zero.  those
  ! --- should be pes that either have no data and may have been allocated
  ! --- before (no harm in doing it again) or have never been allocated

  if (lsizei <= 0) then
    if(len_trim(iList) > 0 .and. len_trim(rList) > 0) then
      call mct_aVect_init(AVin,iList=iList,rList=rList,lsize=lsizen)
    elseif (len_trim(iList) > 0 .and. len_trim(rList) == 0) then
      call mct_aVect_init(AVin,iList=iList,lsize=lsizen)
    elseif (len_trim(iList) == 0 .and. len_trim(rList) > 0) then
      call mct_aVect_init(AVin,rList=rList,lsize=lsizen)
    endif
  endif

end subroutine seq_rearr_avExtend

!=======================================================================

subroutine seq_rearr_avCreate(AVin,IDin,AVout,ID,lsize)

  !-----------------------------------------------------------------------
  ! Extend an AV to a larger set of pes or
  ! Initialize an AV on another set of pes
  !-----------------------------------------------------------------------

  implicit none
  type(mct_aVect), intent(INOUT):: AVin
  integer         ,intent(IN)   :: IDin ! ID associated with AVin
  type(mct_aVect), intent(INOUT):: AVout
  integer        , intent(IN)   :: ID   ! ID to initialize over
  integer        , intent(IN)   :: lsize

  ! Local variables

  character(len=*),parameter :: subname = "(seq_rearr_avCreate) "
  integer :: mpicom
  integer :: rank,rank2
  integer :: lsizei, lsizen
  integer :: srank,srankg
  integer :: ierr
  integer :: nints
  character(len=1024) :: iList,rList

  call seq_comm_setptrs(ID,mpicom=mpicom,iam=rank)

  ! --- lsizen is the size of the newly initialized AV, zero is valid

  lsizei = -1  
  if (seq_comm_iamin(IDin)) lsizei = mct_aVect_lsize(AVin)
  lsizen = lsize

  ! --- find a pe that already has AVin allocated, use MPI_MAX to do so
  ! --- set the pe and broadcast it to all other pes

  srank = -1
  srankg = -1
  if (lsizei > 0) srank = rank

  call mpi_allreduce(srank,srankg,1,MPI_INTEGER,MPI_MAX,mpicom,ierr)
  call shr_mpi_chkerr(ierr,subname//' mpi_allreduce max')

  if (srankg < 0) then
    write(logunit,*) subname,' ERROR AVin not initialized '
    call shr_sys_abort()
  endif

  ! --- set the iList and rList from the broadcast pe (srankg) and 
  ! --- broadcast the lists

  iList = " "
  rList = " "
  if (rank == srankg) then
    if (mct_aVect_nIAttr(AVin) /= 0) iList = mct_aVect_ExportIList2c(AVin)
    if (mct_aVect_nRattr(AVin) /= 0) rList = mct_aVect_ExportRList2c(AVin)
  endif

  call mpi_bcast(iList,len(iList),MPI_CHARACTER,srankg,mpicom,ierr)
  call mpi_bcast(rList,len(rList),MPI_CHARACTER,srankg,mpicom,ierr)

  ! --- now allocate the AV on all pes.  the AV should not exist before.
  ! --- If it does, mct should die.

  if(len_trim(iList) > 0 .and. len_trim(rList) > 0) then
    call mct_aVect_init(AVout,iList=iList,rList=rList,lsize=lsizen)
  elseif (len_trim(iList) > 0 .and. len_trim(rList) == 0) then
    call mct_aVect_init(AVout,iList=iList,lsize=lsizen)
  elseif (len_trim(iList) == 0 .and. len_trim(rList) > 0) then
    call mct_aVect_init(AVout,rList=rList,lsize=lsizen)
  endif

end subroutine seq_rearr_avCreate

!=======================================================================
logical function seq_rearr_gsmapIdentical(gsmap1,gsmap2)

  implicit none
  type(mct_gsMap), intent(IN):: gsmap1
  type(mct_gsMap), intent(IN):: gsmap2

  ! Local variables

  character(len=*),parameter :: subname = "(seq_rearr_gsmapIdentical) "
  integer :: n
  logical :: identical

  !-----------------------

  identical = .true.

  ! --- continue compare ---
  if (identical) then
     if (mct_gsMap_gsize(gsmap1) /= mct_gsMap_gsize(gsmap2)) identical = .false.
     if (mct_gsMap_ngseg(gsmap1) /= mct_gsMap_ngseg(gsmap2)) identical = .false.
  endif

  ! --- continue compare ---
  if (identical) then
     do n = 1,mct_gsMap_ngseg(gsmap1)
        if (gsmap1%start(n)  /= gsmap2%start(n) ) identical = .false.
        if (gsmap1%length(n) /= gsmap2%length(n)) identical = .false.
        if (gsmap1%pe_loc(n) /= gsmap2%pe_loc(n)) identical = .false.
     enddo
  endif

  seq_rearr_gsmapIdentical = identical

end function seq_rearr_gsmapIdentical
    
!=======================================================================

end module seq_rearr_mod
