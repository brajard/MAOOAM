
! rk2_WL_integrator.f90
!
!>  Module with the WL rk2 integration routines.
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

MODULE rk2_WL_integrator
  USE params, only: ndim,natm
  USE stoch_params, only:  q_ar,q_au,q_or,q_ou,muti,dts,dtsn,mode
  USE tensor, only: sparse_mul2_k,sparse_mul3
  USE WL_tensor
  USE aotensor_def, only:aotensor
  USE MAR, only: init_MAR,MAR_step,ms
  USE stoch_mod
  USE memory
  USE rk2_ss_integrator
  IMPLICIT NONE

  PRIVATE

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: buf_y1,buf_f0,buf_f1 !< Integration buffers
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: buf_M2,buf_M1,buf_M3,buf_M,buf_M3s !< Dummy buffers holding the terms /f$M_i\f$ of the parameterization
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: anoise          !< Additive noise term
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: dWar,dWau,dWor,dWou  !< Standard gaussian noise buffers

  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: x1 !< Buffer holding the subsequent states of the first MAR
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: x2 !< Buffer holding the subsequent states of the second MAR

  PUBLIC :: init_integrator, step, full_step

CONTAINS
  !> Subroutine that initialize the MARs, the memory unit and the integration buffers
  SUBROUTINE init_integrator
    INTEGER :: AllocStat,i

    CALL init_ss_integrator
    
    print*, 'Initializing the integrator ...'

    IF (mode.ne.'ures') THEN
       PRINT*, '*** Mode set to ',mode,' in stoch_params.nml ***'
       PRINT*, '*** WL configuration only support unresolved mode ***'
       STOP "*** Please change to 'ures' and perform the configuration again ! ***"
    ENDIF
    
    ALLOCATE(buf_y1(0:ndim),buf_f0(0:ndim),buf_f1(0:ndim),STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    ALLOCATE(buf_M1(0:ndim), buf_M2(0:ndim), buf_M3(0:ndim), buf_M(0:ndim), buf_M3s(0:ndim), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    ALLOCATE(dWar(0:ndim),dWau(0:ndim),dWor(0:ndim),dWou(0:ndim),STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    ALLOCATE(anoise(0:ndim), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    buf_y1=0.D0
    buf_f1=0.D0
    buf_f0=0.D0

    dWar=0.D0
    dWor=0.D0
    dWau=0.D0
    dWou=0.D0

    buf_M1=0.D0
    buf_M2=0.D0
    buf_M3=0.D0
    buf_M3s=0.D0
    buf_M=0.D0

    print*, 'Initializing the MARs ...'

    CALL init_MAR

    ALLOCATE(x1(0:ndim,ms), x2(0:ndim,ms), STAT=AllocStat)

    x1=0.D0
    DO i=1,50000
       CALL MAR_step(x1)
    ENDDO
    
    x2=0.D0
    DO i=1,50000
       CALL MAR_step(x2)
    ENDDO

    CALL init_memory

  END SUBROUTINE init_integrator
  
  !> Routine to compute the \f$M_1\f$ term
  !> @param y Present state of the WL system
  SUBROUTINE compute_M1(y)
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(IN) :: y
    buf_M1=0.D0
    IF (M12def) CALL sparse_mul2_k(M12, y, buf_M1)
    buf_M1=buf_M1+M1tot
  END SUBROUTINE compute_M1

  !> Routine to compute the \f$M_2\f$ term
  !> @param y Present state of the WL system
  SUBROUTINE compute_M2(y)
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(IN) :: y
    buf_M=0.D0
    buf_M2=0.D0
    IF (M21def) CALL sparse_mul3(M21, y, x1(0:ndim,1), buf_M)
    IF (M22def) CALL sparse_mul3(M22, x2(0:ndim,1), x2(0:ndim,1), buf_M2)
    buf_M2=buf_M2+buf_M
  END SUBROUTINE compute_M2
  
  !> Routine to perform an integration step (Heun algorithm) of the WL system.
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
    INTEGER :: i

    IF (mod(t,muti)<dt) THEN 
       CALL compute_M3(y,dts,dtsn,.true.,.true.,.true.,muti/2,buf_M3s)
       buf_M3=buf_M3s
       DO i=1,1
          CALL compute_M3(y,dts,dtsn,.false.,.true.,.true.,muti/2,buf_M3s)
          buf_M3=buf_M3+buf_M3s
       ENDDO
       !DO i=1,2
       !   CALL compute_M3(y,dts,dtsn,.false.,.true.,.true.,muti/2,buf_M3s)
       !   buf_M3=buf_M3+buf_M3s
       !ENDDO
       buf_M3=buf_M3/2
    ENDIF


    CALL stoch_atm_res_vec(dWar)
    CALL stoch_oc_res_vec(dWor)
    anoise=(q_ar*dWar+q_or*dWor)*dtn

    CALL tendencies(t,y,buf_f0)
    CALL MAR_step(x1)
    CALL MAR_step(x2)
    CALL compute_M1(y)
    CALL compute_M2(y)
    buf_f0= buf_f0+buf_M1+buf_M2+buf_M3
    buf_y1 = y+dt*buf_f0+anoise

    CALL tendencies(t+dt,buf_y1,buf_f1)
    CALL compute_M1(buf_y1)
    CALL compute_M2(buf_y1)
    !IF (mod(t,muti)<dt) CALL compute_M3(buf_y1,dts,dtsn,.false.,.true.,buf_M3)

    buf_f0=0.5*(buf_f0+buf_f1+buf_M1+buf_M2+buf_M3)
    res=y+dt*buf_f0+anoise

    tend=buf_M3
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

END MODULE rk2_WL_integrator
