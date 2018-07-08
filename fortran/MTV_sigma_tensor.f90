
! MTV_sigma_tensor.f90 
!
!> The MTV noise sigma matrices used to integrate the MTV model   
!
!> @copyright                                                               
!> 2018 Jonathan Demaeyer.
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------!
!                                                                           
!> @remark                                                                 
!> See : Franzke, C., Majda, A. J., & Vanden-Eijnden, E. (2005). Low-order  
!>       stochastic mode reduction for a realistic barotropic model climate.
!>       Journal of the atmospheric sciences, 62(6), 1722-1745.             
!                                                                           
!---------------------------------------------------------------------------!

MODULE sigma
  USE params, only:ndim
  USE MTV_int_tensor, only: Q1,Q2,Utot,Vtot
  USE tensor
  USE util, only: printmat,cprintmat,reduce,ireduce
  USE sqrt_mod, only: sqrtm,init_sqrt,chol,sqrtm_svd
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: init_sigma,compute_mult_sigma

  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: sig1 !< \f$\sigma_1(X)\f$ state-dependent noise matrix
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: sig2 !< \f$\sigma_2\f$ state-independent noise matrix
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: sig1r !< Reduced \f$\sigma_1(X)\f$ state-dependent noise matrix
  
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dumb_mat1 !< Dummy matrix
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dumb_mat2 !< Dummy matrix
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dumb_mat3 !< Dummy matrix
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dumb_mat4 !< Dummy matrix
  INTEGER, DIMENSION(:), ALLOCATABLE :: ind1,rind1,ind2,rind2 !< Reduction indices

  INTEGER :: n1,n2

CONTAINS
 
   
  !> Subroutine to initialize the sigma matices
  SUBROUTINE init_sigma(mult,Q1fill)
    LOGICAL, INTENT(OUT) :: mult,Q1fill
    INTEGER :: AllocStat,info1,info2

    CALL init_sqrt

    ALLOCATE(sig1(ndim,ndim), sig2(ndim,ndim), sig1r(ndim,ndim),STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    ALLOCATE(ind1(ndim), rind1(ndim), ind2(ndim), rind2(ndim), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    ALLOCATE(dumb_mat1(ndim,ndim), dumb_mat2(ndim,ndim), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    ALLOCATE(dumb_mat3(ndim,ndim), dumb_mat4(ndim,ndim), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    print*, "Initializing the sigma matrices"

    CALL reduce(Q2,dumb_mat1,n2,ind2,rind2)
    IF (n2 /= 0) THEN
       CALL sqrtm_svd(dumb_mat1(1:n2,1:n2),dumb_mat2(1:n2,1:n2),info1,info2,min(max(n2/2,2),64))
       CALL ireduce(sig2,dumb_mat2,n2,ind2,rind2)
    ELSE
       sig2=0.D0
    ENDIF

    mult=(.not.((tensor_empty(Utot)).and.(tensor4_empty(Vtot))))
    Q1fill=.true.
    CALL reduce(Q1,dumb_mat1,n1,ind1,rind1)
    IF (n1 /= 0) THEN
       
       CALL sqrtm_svd(dumb_mat1(1:n1,1:n1),dumb_mat2(1:n1,1:n1),info1,info2,min(max(n1/2,2),64))
       CALL ireduce(sig1,dumb_mat2,n1,ind1,rind1)
    ELSE
       Q1fill=.false.
       sig1=0.D0
    ENDIF
    sig1r=sig1

  END SUBROUTINE init_sigma

  !> Routine to actualize the matrix \f$\sigma_1\f$ based on the state y of the MTV system
  !> @param y State of the MTV system 
  SUBROUTINE compute_mult_sigma(y)
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(IN) :: y
    INTEGER :: info,info2
    CALL sparse_mul3_mat(Utot,y,dumb_mat1)
    CALL sparse_mul4_mat(Vtot,y,y,dumb_mat2)
    dumb_mat3=dumb_mat1+dumb_mat2+Q1
    CALL reduce(dumb_mat3,dumb_mat1,n1,ind1,rind1)
    IF (n1 /= 0) THEN
       CALL sqrtm_svd(dumb_mat1(1:n1,1:n1),dumb_mat2(1:n1,1:n1),info,info2,min(max(n1/2,2),64))
       ! dumb_mat2=0.D0
       ! CALL chol(0.5*(dumb_mat1(1:n1,1:n1)+transpose(dumb_mat1(1:n1,1:n1))),dumb_mat2(1:n1,1:n1),info)
       IF ((.not.ANY(ISNAN(dumb_mat2))).and.(info.eq.0).and.(.not.ANY(dumb_mat2>HUGE(0.D0)))) THEN
          CALL ireduce(sig1,dumb_mat2,n1,ind1,rind1)
       ELSE
          sig1=sig1r
       ENDIF
    ELSE
       sig1=sig1r
    ENDIF
  END SUBROUTINE compute_mult_sigma


END MODULE sigma
