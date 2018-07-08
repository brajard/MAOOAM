
! rk2_ss_integrator.f90
!
!>  Module with the stochastic uncoupled resolved nonlinear and tangent linear 
!>  rk2 dynamics integration routines.
!
!> @copyright                                                               
!> 2018 Jonathan Demaeyer.
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------!
!                                                                           
!>  @remark                                                                 
!>  This module actually contains the Heun algorithm routines.
!                                                                           
!---------------------------------------------------------------------------

MODULE rk2_ss_integrator
  USE params, only: ndim,natm
  USE stoch_params, only: q_ar,q_au,q_or,q_ou
  USE tensor, only: sparse_mul3
  USE dec_tensor, only: ss_tensor,ss_tl_tensor
  USE stoch_mod
  IMPLICIT NONE

  PRIVATE

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: dWar,dWor !< Standard gaussian noise buffers

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: anoise !< Additive noise term

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: buf_y1,buf_f0,buf_f1 !< Integration buffers

  PUBLIC :: init_ss_integrator,ss_step,tendencies,ss_tl_step,tl_tendencies

CONTAINS

  !> Subroutine to initialize the uncoupled resolved rk2 integrator
  SUBROUTINE init_ss_integrator
    INTEGER :: AllocStat

    ALLOCATE(buf_y1(0:ndim),buf_f0(0:ndim),buf_f1(0:ndim),STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    ALLOCATE(anoise(0:ndim),STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    ALLOCATE(dWar(0:ndim),dWor(0:ndim),STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    dWar=0.D0
    dWor=0.D0

  END SUBROUTINE init_ss_integrator

  !> Routine computing the tendencies of the uncoupled resolved model
  !> @param t Time at which the tendencies have to be computed. Actually not needed for autonomous systems.
  !> @param y Point at which the tendencies have to be computed.
  !> @param res vector to store the result.
  !> @remark Note that it is NOT safe to pass `y` as a result buffer, 
  !> as this operation does multiple passes.
  SUBROUTINE tendencies(t,y,res)
    REAL(KIND=8), INTENT(IN) :: t
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(IN) :: y
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(OUT) :: res
    CALL sparse_mul3(ss_tensor, y, y, res)
  END SUBROUTINE tendencies

  !> Tendencies for the tangent linear model of the uncoupled resolved dynamics
  !> in point ystar for perturbation deltay.
  !> @param t time
  !> @param y point of the tangent space at which the tendencies have to be computed.
  !> @param ys point in trajectory to which the tangent space belongs.
  !> @param res vector to store the result.
  SUBROUTINE tl_tendencies(t,y,ys,res)
    REAL(KIND=8), INTENT(IN) :: t
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(IN) :: y,ys
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(OUT) :: res
    CALL sparse_mul3(ss_tl_tensor, y, ys, res)
  END SUBROUTINE tl_tendencies

  !> Routine to perform a stochastic integration step of the unresolved
  !> uncoupled dynamics (Heun algorithm).
  !> The incremented time is returned.
  !> @param y Initial point.
  !> @param ys Dummy argument for compatibility.
  !> @param t Actual integration time
  !> @param dt Integration timestep.
  !> @param dtn Stochastic integration timestep (normally square-root of dt).
  !> @param res Final point after the step.
  SUBROUTINE ss_step(y,ys,t,dt,dtn,res)
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(IN) :: y,ys
    REAL(KIND=8), INTENT(INOUT) :: t
    REAL(KIND=8), INTENT(IN) :: dt,dtn
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(OUT) :: res

    CALL stoch_atm_res_vec(dWar)
    CALL stoch_oc_res_vec(dWor)
    anoise=(q_ar*dWar+q_or*dWor)*dtn
    CALL tendencies(t,y,buf_f0)
    buf_y1 = y+dt*buf_f0+anoise
    CALL tendencies(t,buf_y1,buf_f1)
    res=y+0.5*(buf_f0+buf_f1)*dt+anoise
    t=t+dt
  END SUBROUTINE ss_step

  !> Routine to perform a stochastic integration step of the unresolved
  !> uncoupled tangent linear dynamics (Heun algorithm).
  !> The incremented time is returned.
  !> @param y Initial point.
  !> @param ys point in trajectory to which the tangent space belongs.
  !> @param t Actual integration time
  !> @param dt Integration timestep.
  !> @param dtn Stochastic integration timestep (normally square-root of dt).
  !> @param res Final point after the step.
  SUBROUTINE ss_tl_step(y,ys,t,dt,dtn,res)
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(IN) :: y,ys
    REAL(KIND=8), INTENT(INOUT) :: t
    REAL(KIND=8), INTENT(IN) :: dt,dtn
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(OUT) :: res

    CALL stoch_atm_res_vec(dWar)
    CALL stoch_oc_res_vec(dWor)
    anoise=(q_ar*dWar+q_or*dWor)*dtn
    CALL tl_tendencies(t,y,ys,buf_f0)
    buf_y1 = y+dt*buf_f0+anoise
    CALL tl_tendencies(t,buf_y1,ys,buf_f1)
    res=y+0.5*(buf_f0+buf_f1)*dt+anoise
    t=t+dt
  END SUBROUTINE ss_tl_step

END MODULE rk2_ss_integrator
