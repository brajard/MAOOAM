
! rk2_stoch_integrator.f90
!
!>  Module with the stochastic rk2 integration routines.
!
!> @copyright                                                               
!> 2018 Jonathan Demaeyer.
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------!
!                                                                           
!>  @remark                                                                 
!>  This module actually contains the Heun algorithm routines.
!>  There are four modes for this integrator:
!>         - full: use the full dynamics
!>         - ures: use the intrinsic unresolved dynamics
!>         - qfst: use the quadratic terms of the unresolved tendencies
!>         - reso: use the resolved dynamics alone
!                                                                           
!---------------------------------------------------------------------------

MODULE rk2_stoch_integrator
  USE params, only: ndim,natm
  USE stoch_params, only: q_ar,q_au,q_or,q_ou,mode
  USE tensor, only:sparse_mul3,coolist,copy_coo
  USE aotensor_def, only: aotensor
  USE dec_tensor, only: ss_tensor,ff_tensor,Byyy
  USE stoch_mod
  IMPLICIT NONE

  PRIVATE

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: dWar,dWau,dWor,dWou !< Standard gaussian noise buffers

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: buf_y1,buf_f0,buf_f1 !< Integration buffers

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: anoise !< Additive noise term

  TYPE(coolist), DIMENSION(:), ALLOCATABLE :: int_tensor !< Dummy tensor that will hold the tendencies tensor

  PUBLIC :: init_integrator, step

CONTAINS
  
  !> Subroutine to initialize the integrator
  !> @param force Parameter to force the mode of the integrator
  SUBROUTINE init_integrator(force)
    INTEGER :: AllocStat
    CHARACTER*4, INTENT(IN), OPTIONAL :: force
    CHARACTER*4 :: test

    ALLOCATE(buf_y1(0:ndim),buf_f0(0:ndim),buf_f1(0:ndim),STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    ALLOCATE(anoise(0:ndim),STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    ALLOCATE(dWar(0:ndim),dWau(0:ndim),dWor(0:ndim),dWou(0:ndim),STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    ALLOCATE(int_tensor(ndim),STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    dWar=0.D0
    dWor=0.D0
    dWau=0.D0
    dWou=0.D0

    IF (PRESENT(force)) THEN
       test=force
    ELSE
       test=mode
    ENDIF

    SELECT CASE (test)
    CASE('full')
       CALL copy_coo(aotensor,int_tensor)
    CASE('ures')
       CALL copy_coo(ff_tensor,int_tensor)
    CASE('qfst')
       CALL copy_coo(Byyy,int_tensor)
    CASE('reso')
       CALL copy_coo(ss_tensor,int_tensor)
    CASE DEFAULT
       STOP '*** MODE variable not properly defined ***'
    END SELECT

  END SUBROUTINE init_integrator
  
  !> Routine computing the tendencies of the selected model 
  !> @param t Time at which the tendencies have to be computed. Actually not needed for autonomous systems.
  !> @param y Point at which the tendencies have to be computed.
  !> @param res vector to store the result.
  !> @remark Note that it is NOT safe to pass `y` as a result buffer, 
  !> as this operation does multiple passes.
  SUBROUTINE tendencies(t,y,res)
    REAL(KIND=8), INTENT(IN) :: t
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(IN) :: y
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(OUT) :: res
    CALL sparse_mul3(int_tensor, y, y, res)
  END SUBROUTINE tendencies
  
  !> Routine to perform a stochastic step of the selected dynamics (Heun algorithm).
  !> The incremented time is returned.
  !> @param y Initial point.
  !> @param t Actual integration time
  !> @param dt Integration timestep.
  !> @param dtn Stochastic integration timestep (normally square-root of dt).
  !> @param res Final point after the step.
  !> @param tend Partial or full tendencies used to perform the step (used for debugging).
  SUBROUTINE step(y,t,dt,dtn,res,tend)
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(IN) :: y
    REAL(KIND=8), INTENT(INOUT) :: t
    REAL(KIND=8), INTENT(IN) :: dt,dtn
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(OUT) :: res,tend
    
    CALL stoch_atm_res_vec(dWar)
    CALL stoch_atm_unres_vec(dWau)
    CALL stoch_oc_res_vec(dWor)
    CALL stoch_oc_unres_vec(dWou)
    anoise=(q_ar*dWar+q_au*dWau+q_or*dWor+q_ou*dWou)*dtn
    CALL tendencies(t,y,buf_f0)
    CALL sparse_mul3(int_tensor,y,y,tend)
    buf_y1 = y+dt*buf_f0+anoise
    CALL sparse_mul3(int_tensor,buf_y1,buf_y1,buf_f1)
    tend=0.5*(tend+buf_f1)
    CALL tendencies(t,buf_y1,buf_f1)
    res=y+0.5*(buf_f0+buf_f1)*dt+anoise
    t=t+dt
  END SUBROUTINE step

END MODULE rk2_stoch_integrator
