MODULE input_vars

    USE global

    implicit none

! Optimization variables:
    character(len=6) :: rtimePRE_str
    character(len=6) :: rtimeST1_str
    character(len=6) :: rtimeST2_str
    character(len=6) :: rtimeST3_str
    integer(i4b) :: maxiter_PRE
    integer(i4b) :: maxiter_ST1
    integer(i4b) :: maxiter_ST2
    integer(i4b) :: maxiter_ST3
!    integer(i4b) :: strategy_pre  ! move to global
!    integer(i4b) :: strategy
    real(DP) :: RMStgt_pre
    real(DP) :: RMStgt_de
    real(DP) :: RMStgt_sres
    real(DP) :: RMStgt_isres
    real(DP) :: RMSrej_Q23
    real(DP) :: THdRMSE
    real(DP) :: THdObjF
    real(DP) :: isres_alpha

    
#ifdef ERGO
    integer(i4b) :: DE_pop
    integer(i4b) :: DE_rndscale
    integer(i4b) :: DE_strat
    real(DP) :: DE_F_CR
#endif

END MODULE input_vars
