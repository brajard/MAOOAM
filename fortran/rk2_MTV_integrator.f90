
! rk2_MTV_integrator.f90
!
!>  Module with the MTV rk2 integration routines.
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

MODULE rk2_MTV_integrator
  USE params, only: ndim
  USE stoch_params, only: q_ar,q_au,q_or,q_ou,mnuti
  USE tensor, only: sparse_mul2_k,sparse_mul3,sparse_mul4
  USE MTV_int_tensor, only: Ltot,Htot,Mtot,Btot
  USE aotensor_def, only:aotensor
  USE sigma
  USE stoch_mod
  USE rk2_ss_integrator
  IMPLICIT NONE

  PRIVATE

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: buf_y1,buf_f0,buf_f1 !< Integration buffers
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: dW, dWmult           !< Standard gaussian noise buffers
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: dWar,dWau,dWor,dWou  !< Standard gaussian noise buffers
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: anoise,noise         !< Additive noise term
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: noisemult            !< Multiplicative noise term
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: G                    !< G term of the MTV tendencies
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: buf_G                !< Buffer for the G term computation

  LOGICAL :: mult                                                 !< Logical indicating if the sigma1 matrix must be computed for every state change
  LOGICAL :: Q1fill                                               !< Logical indicating if the matrix Q1 is non-zero
  LOGICAL :: compute_mult                                         !< Logical indicating if the Gaussian noise for the multiplicative noise must be computed

  REAL(KIND=8), PARAMETER :: sq2 = sqrt(2.D0)                     !< Hard coded square root of 2
  
  PUBLIC :: init_integrator, step, full_step

CONTAINS

  !> Subroutine to initialize the MTV rk2 integrator
  SUBROUTINE init_integrator
    INTEGER :: AllocStat

    CALL init_ss_integrator ! Initialize the uncoupled resolved dynamics

    ALLOCATE(buf_y1(0:ndim),buf_f0(0:ndim),buf_f1(0:ndim),STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    buf_y1=0.D0
    buf_f1=0.D0
    buf_f0=0.D0
    
    print*, 'Initializing the integrator ...'
    CALL init_sigma(mult,Q1fill)
    CALL init_noise
    CALL init_G
  END SUBROUTINE init_integrator

  !> Routine to initialize the noise vectors and buffers
  SUBROUTINE init_noise
    INTEGER :: AllocStat
    ALLOCATE(dW(0:ndim), dWmult(0:ndim), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    ALLOCATE(dWar(0:ndim),dWau(0:ndim),dWor(0:ndim),dWou(0:ndim),STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    ALLOCATE(anoise(0:ndim), noise(0:ndim), noisemult(0:ndim), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    dW=0.D0
    dWmult=0.D0

    dWar=0.D0
    dWor=0.D0
    dWau=0.D0
    dWou=0.D0

    anoise=0.D0
    noise=0.D0
    noisemult=0.D0

    compute_mult=((Q1fill).OR.(mult))
    
  END SUBROUTINE init_noise
  
  !> Routine to initialize the G term
  SUBROUTINE init_G
    INTEGER :: AllocStat
    ALLOCATE(G(0:ndim), buf_G(0:ndim), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
  END SUBROUTINE init_G

  !> Routine to actualize the G term based on the state y of the MTV system
  !> @param y State of the MTV system 
  SUBROUTINE compG(y)
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(IN) :: y

    G=Htot
    CALL sparse_mul2_k(Ltot,y,buf_G)
    G=G+buf_G
    CALL sparse_mul3(Btot,y,y,buf_G)
    G=G+buf_G
    CALL sparse_mul4(Mtot,y,y,y,buf_G)
    G=G+buf_G
  END SUBROUTINE compG
 
  !> Routine to perform an integration step (Heun algorithm) of the MTV system. The incremented time is returned.
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

    CALL compG(y)
    
    CALL stoch_atm_res_vec(dWar)
    CALL stoch_oc_res_vec(dWor)
    anoise=q_ar*dWar+q_or*dWor
    CALL stoch_vec(dW)
    IF (compute_mult) CALL stoch_vec(dWmult)
    noise(1:ndim)=matmul(sig2,dW(1:ndim))
    IF ((mult).and.(mod(t,mnuti)<dt)) CALL compute_mult_sigma(y)
    IF (compute_mult) noisemult(1:ndim)=matmul(sig1,dWmult(1:ndim))

    CALL tendencies(t,y,buf_f0)
    buf_y1 = y+dt*(buf_f0+G)+(anoise+sq2*(noise+noisemult))*dtn

    buf_f1=G
    CALL compG(buf_y1)
    G=0.5*(G+buf_f1)
    
    IF ((mult).and.(mod(t,mnuti)<dt)) CALL compute_mult_sigma(buf_y1)
    IF (compute_mult) THEN
       buf_f1(1:ndim)=matmul(sig1,dWmult(1:ndim))
       noisemult(1:ndim)=0.5*(noisemult(1:ndim)+buf_f1(1:ndim))
    ENDIF


    CALL tendencies(t,buf_y1,buf_f1)
    buf_f0=0.5*(buf_f0+buf_f1)
    res=y+dt*(buf_f0+G)+(anoise+sq2*(noise+noisemult))*dtn
    ! tend=G+sq2*(noise+noisemult)/dtn
    tend=sq2*noisemult/dtn
    t=t+dt

  END SUBROUTINE step

  !> Routine to perform an integration step (Heun algorithm) of the full stochastic system. The incremented time is returned.
  !> @param y Initial point.
  !> @param t Actual integration time
  !> @param dt Integration timestep.
  !> @param dtn Stochastoc integration timestep (normally square-root of dt).
  !> @param res Final point after the step.
  SUBROUTINE full_step(y,t,dt,dtn,res)
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(IN) :: y
    REAL(KIND=8), INTENT(INOUT) :: t
    REAL(KIND=8), INTENT(IN) :: dt,dtn
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(OUT) :: res
    CALL stoch_atm_res_vec(dWar)
    CALL stoch_atm_unres_vec(dWau)
    CALL stoch_oc_res_vec(dWor)
    CALL stoch_oc_unres_vec(dWou)
    anoise=(q_ar*dWar+q_au*dWau+q_or*dWor+q_ou*dWou)*dtn
    CALL sparse_mul3(aotensor,y,y,buf_f0)
    buf_y1 = y+dt*buf_f0+anoise
    CALL sparse_mul3(aotensor,buf_y1,buf_y1,buf_f1)
    res=y+0.5*(buf_f0+buf_f1)*dt+anoise
    t=t+dt
  END SUBROUTINE full_step

END MODULE rk2_MTV_integrator
