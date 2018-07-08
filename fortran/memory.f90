
! memory.f90
!
!> Module that compute the memory term \f$M_3\f$ of the WL parameterization.
!
!> @copyright                                                               
!> 2018 Jonathan Demaeyer.
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------!
!                                                                           
!> @remark                                                                 
!                                                                           
!---------------------------------------------------------------------------!

MODULE memory

  USE tensor
  USE WL_tensor
  USE params, only: ndim
  USE stoch_params, only: muti,mems,x_int_mode
  USE rk2_ss_integrator, only: ss_step,ss_tl_step
  
  IMPLICIT NONE

  PRIVATE

  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: X !< Array storing the previous state of the system
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: Xs !< Array storing the resolved time evolution of the previous state of the system
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: Zs !< Dummy array to replace Xs in case where the evolution is not stored
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: buf_M !< Dummy vector
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: buf_M3 !< Dummy vector to store the \f$M_3\f$ integrand

  INTEGER :: t_index !< Integer storing the time index (current position in the arrays)

  PROCEDURE(ss_step), POINTER :: step !< Procedural pointer pointing on the resolved dynamics step routine

  PUBLIC :: init_memory,compute_M3,test_M3


CONTAINS

  !> Subroutine to initialise the memory
  SUBROUTINE init_memory
    INTEGER :: AllocStat

    t_index=mems

    ALLOCATE(X(0:ndim,mems), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    X=0.D0
    
    IF (B23def) THEN
       ALLOCATE(Xs(0:ndim,mems), Zs(0:ndim,mems), STAT=AllocStat)
       IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

       Xs=0.D0
    ENDIF

    ALLOCATE(buf_M3(0:ndim), buf_M(0:ndim), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
    
    SELECT CASE (x_int_mode)
    CASE('reso')
       step => ss_step
    CASE('tang')
       step => ss_tl_step
    CASE DEFAULT
       STOP '*** X_INT_MODE variable not properly defined in stoch_params.nml ***'
    END SELECT
    
  END SUBROUTINE init_memory

  !> Compute the integrand of \f$M_3\f$ at each time in the past and integrate
  !> to get the memory term
  !> @param y current state
  !> @param dt timestep
  !> @param dtn stochastic timestep
  !> @param savey set if the state is stored in X at the end
  !> @param save_ev set if the result of the resolved time evolution is stored in Xs at the end
  !> @param evolve set if the resolved time evolution is performed
  !> @param inter set over which time interval the resolved time evolution must be computed 
  !> @param h_int result of the integration - give the memory term
  SUBROUTINE compute_M3(y,dt,dtn,savey,save_ev,evolve,inter,h_int)
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(IN) :: y
    REAL(KIND=8), INTENT(IN) :: dt,dtn
    LOGICAL, INTENT(IN) :: savey,save_ev,evolve
    REAL(KIND=8), INTENT(IN) :: inter 
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(OUT) :: h_int
    REAL(KIND=8) :: t
    INTEGER :: i,j
    
    X(:,t_index)=y
    IF (B23def) THEN
       Xs(:,t_index)=y
       Zs(:,t_index)=y
       DO i=1,mems-1
          j=modulo(t_index+i-1,mems)+1
          Zs(:,j)=Xs(:,j)
          IF (evolve) THEN
             IF (dt.lt.inter) THEN
                t=0.D0
                DO WHILE (t+dt<inter)
                   CALL step(Zs(:,j),y,t,dt,dtn,Zs(:,j))
                ENDDO
                CALL step(Zs(:,j),y,t,inter-t,sqrt(inter-t),Zs(:,j))
             ELSE
                CALL step(Zs(:,j),y,t,inter,sqrt(inter),Zs(:,j))
             ENDIF
          ENDIF
          IF (save_ev) Xs(:,j)=Zs(:,j)
       ENDDO
    ENDIF
    

    ! Computing the integral
    h_int=0.D0

    DO i=1,mems
       j=modulo(t_index+i-2,mems)+1
       buf_M3=0.D0
       IF (Ldef) THEN
          CALL sparse_mul3(Ltot(:,i),y,X(:,j),buf_M)
          buf_M3=buf_M3+buf_M
       ENDIF
       IF (B14def) THEN
          CALL sparse_mul3(B14(:,i),X(:,j),X(:,j),buf_M)
          buf_M3=buf_M3+buf_M
       ENDIF
       IF (B23def) THEN
          CALL sparse_mul3(B23(:,i),X(:,j),Zs(:,j),buf_M)
          buf_M3=buf_M3+buf_M
       ENDIF
       IF (Mdef) THEN
          CALL sparse_mul4(Mtot(:,i),X(:,j),X(:,j),Zs(:,j),buf_M)
          buf_M3=buf_M3+buf_M
       ENDIF
       IF ((i.eq.1).or.(i.eq.mems)) THEN
          h_int=h_int+0.5*buf_M3
       ELSE
          h_int=h_int+buf_M3
       ENDIF
    ENDDO

    h_int=muti*h_int
    IF (savey) THEN
       t_index=t_index-1
       IF (t_index.eq.0) t_index=mems
    ENDIF
  END SUBROUTINE compute_M3

  !> Routine to test the #compute_M3 routine
  !> @param y current state
  !> @param dt timestep
  !> @param dtn stochastic timestep
  !> @param h_int result of the integration - give the memory term
  SUBROUTINE test_M3(y,dt,dtn,h_int)
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(IN) :: y
    REAL(KIND=8), INTENT(IN) :: dt,dtn
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(OUT) :: h_int
    INTEGER :: i,j
    
    CALL compute_M3(y,dt,dtn,.true.,.true.,.true.,muti,h_int)
    print*, t_index
    print*, 'X'
    DO i=1,mems
       j=modulo(t_index+i-1,mems)+1
       print*, i,j,X(1,j)
    ENDDO
    
    IF (B23def) THEN
       print*, 'Xs'
       DO i=1,mems
          j=modulo(t_index+i-1,mems)+1
          print*, i,j,Xs(1,j)
       ENDDO
    ENDIF
    print*, 'h_int',h_int
  END SUBROUTINE test_M3

END MODULE memory

