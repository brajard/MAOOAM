
! MAR.f90
!
!> Multidimensional Autoregressive module to generate the correlation for the
!> WL parameterization.
!
!> @copyright                                                               
!> 2018 Jonathan Demaeyer.
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------!
!                                                                           
!> @remark                                                                 
!> Based on the equation
!> \f$y_n = \sum_{i=1}^m y_{n-i} \cdot W_i + Q \cdot \xi_n
!> for an order \f$m\f$ MAR with \xi_n a Gaussian random variable.
!---------------------------------------------------------------------------!

MODULE MAR
  USE params, only: ndim
  USE sf_def, only: n_unres,ind,rind
  USE sqrt_mod, only: init_sqrt,sqrtm
  USE util, only: ireduce
  USE stoch_mod, only: gasdev
  IMPLICIT NONE

  PRIVATE

  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: Q !< Square root of the noise covariance matrix
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: Qred !< Reduce version of Q
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: Rred !< Covariance matrix of the noise
  REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: W !< W_i matrix
  REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: Wred !< Reduce W_i matrix
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: buf_y,dW

  INTEGER :: ms !< order of the MAR

  PUBLIC :: init_MAR,W,Q,ms,MAR_step,MAR_step_red
  PUBLIC :: Wred,Qred,Rred

CONTAINS

  !> Subroutine to initialise the MAR
  SUBROUTINE init_MAR
    INTEGER :: AllocStat,nf,i,info,info2
    INTEGER, DIMENSION(3) :: s

    print*, 'Initializing the MAR integrator...'

    print*, 'Loading the MAR config from files...'

    OPEN(20,file='MAR_R_params.def',status='old')
    READ(20,*) nf,ms
    IF (nf /= n_unres) STOP "*** Dimension in files MAR_R_params.def and sf.nml do not correspond ! ***"
    ALLOCATE(Qred(n_unres,n_unres),Rred(n_unres,n_unres),Wred(ms,n_unres,n_unres), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
    ALLOCATE(Q(ndim,ndim),W(ms,ndim,ndim), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
    ALLOCATE(buf_y(0:ndim), dW(ndim), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
    READ(20,*) Rred
    CLOSE(20)
    
    OPEN(20,file='MAR_W_params.def',status='old')
    READ(20,*) nf,ms
    s=shape(Wred)
    IF (nf /= n_unres) STOP "*** Dimension in files MAR_W_params.def and sf.nml do not correspond ! ***"
    IF (s(1) /= ms) STOP "*** MAR order in files MAR_R_params.def and MAR_W_params.def do not correspond ! ***"
    DO i=1,ms
       READ(20,*) Wred(i,:,:)
    ENDDO
    CLOSE(20)

    CALL init_sqrt
    CALL sqrtm(Rred,Qred,info,info2)
    CALL ireduce(Q,Qred,n_unres,ind,rind)
    
    DO i=1,ms
       CALL ireduce(W(i,:,:),Wred(i,:,:),n_unres,ind,rind)
    ENDDO
    
    ! Kept for internal testing - Uncomment if not needed
    ! DEALLOCATE(Wred,Rred,Qred, STAT=AllocStat)
    ! IF (AllocStat /= 0) STOP "*** Deallocation problem ! ***"

    print*, 'MAR of order',ms,'found!'
    
  END SUBROUTINE init_MAR

  !> Routine to generate one step of the MAR
  !> @param x State vector of the MAR (store the \f$y_i\f$)
  SUBROUTINE MAR_step(x)
    REAL(KIND=8), DIMENSION(0:ndim,ms), INTENT(INOUT) :: x
    INTEGER :: j
    
    CALL stoch_vec(dW)
    buf_y=0.D0
    buf_y(1:ndim)=matmul(Q,dW)
    DO j=1,ms
       buf_y(1:ndim)=buf_y(1:ndim)+matmul(x(1:ndim,j),W(j,:,:))
    ENDDO
    x=eoshift(x,shift=-1,boundary=buf_y,dim=2)
  END SUBROUTINE MAR_step


  !> Routine to generate one step of the reduce MAR
  !> @param xred State vector of the MAR (store the \f$y_i\f$)
  !> @remark For debugging purpose only
  SUBROUTINE MAR_step_red(xred)
    REAL(KIND=8), DIMENSION(0:ndim,ms), INTENT(INOUT) :: xred
    INTEGER :: j
    
    CALL stoch_vec(dW)
    buf_y=0.D0
    buf_y(1:n_unres)=matmul(Qred,dW(1:n_unres))
    DO j=1,ms
       buf_y(1:n_unres)=buf_y(1:n_unres)+matmul(xred(1:n_unres,j),Wred(j,:,:))
    ENDDO
    xred=eoshift(xred,shift=-1,boundary=buf_y,dim=2)
  END SUBROUTINE MAR_step_red

    

  SUBROUTINE stoch_vec(dW)
    REAL(KIND=8), DIMENSION(ndim), INTENT(INOUT) :: dW
    INTEGER :: i
    DO i=1,ndim
       dW(i)=gasdev()
    ENDDO
  END SUBROUTINE stoch_vec



END MODULE MAR
