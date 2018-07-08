
! stoch_params.f90                                                                
!
!>  The stochastic models parameters module. 
!
!> @copyright                                                               
!> 2018 Jonathan Demaeyer.
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------
!                                                                           
!>  @remark                                                                 
!>  
!>  
!                                                                           
!---------------------------------------------------------------------------


MODULE stoch_params
  USE params, only: dt
  IMPLICIT NONE

  PUBLIC

  REAL(KIND=8) :: mnuti              !< Multiplicative noise update time interval
  
  REAL(KIND=8) :: t_trans_stoch      !< Transient time period of the stochastic model evolution
  REAL(KIND=8) :: q_ar               !< Atmospheric resolved component noise amplitude
  REAL(KIND=8) :: q_au               !< Atmospheric unresolved component noise amplitude
  REAL(KIND=8) :: q_or               !< Oceanic resolved component noise amplitude
  REAL(KIND=8) :: q_ou               !< Oceanic unresolved component noise amplitude
  REAL(KIND=8) :: dtn                !< Square root of the timestep
  REAL(KIND=8) :: eps_pert           !< Perturbation parameter for the coupling
  REAL(KIND=8) :: tdelta             !< Time separation parameter

  REAL(KIND=8) :: muti                !< Memory update time interval
  REAL(KIND=8) :: meml               !< Time over which the memory kernel is integrated
  REAL(KIND=8) :: t_trans_mem        !< Transient time period to initialize the memory term
  CHARACTER(len=4) :: x_int_mode     !< Integration mode for the resolved component
  REAL(KIND=8) :: dts                !< Intrisic resolved dynamics time step

  INTEGER :: mems                    !< Number of steps in the memory kernel integral
  REAL(KIND=8) :: dtsn               !< Square root of the intrisic resolved dynamics time step

  REAL(KIND=8) :: maxint             !< Upper integration limit of the correlations
  CHARACTER(LEN=4) :: load_mode      !< Loading mode for the correlations
  CHARACTER(LEN=4) :: int_corr_mode  !< Correlation integration mode

  CHARACTER(len=4) :: mode           !< Stochastic mode parameter




CONTAINS

  !> Stochastic parameters initialization routine 
  SUBROUTINE init_stoch_params

    NAMELIST /mtvparams/ mnuti
    NAMELIST /stparams/ q_ar,q_au,q_or,q_ou,eps_pert,tdelta,t_trans_stoch
    NAMELIST /wlparams/ muti,meml,x_int_mode,dts,t_trans_mem
    NAMELIST /corr_init_mode/ load_mode,int_corr_mode,maxint
    NAMELIST /stoch_int_params/ mode


    OPEN(8, file="stoch_params.nml", status='OLD', recl=80, delim='APOSTROPHE')
    READ(8,nml=mtvparams)
    READ(8,nml=wlparams)
    READ(8,nml=stparams)
    READ(8,nml=stoch_int_params)
    READ(8,nml=corr_init_mode)
    CLOSE(8)

    dtn=sqrt(dt)
    dtsn=sqrt(dts)
    mems=ceiling(meml/muti)

    q_au=q_au/tdelta
    q_ou=q_ou/tdelta

  END SUBROUTINE init_stoch_params
END MODULE stoch_params
