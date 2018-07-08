
! MTV_int_tensor.f90 
!
!> The MTV tensors used to integrate the MTV model   
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

MODULE MTV_int_tensor

  !-----------------------------------------------------!
  !                                                     !
  ! Preamble and variables declaration                  !
  !                                                     !
  !-----------------------------------------------------!

  USE tensor
  USE dec_tensor
  USE corrmod
  USE int_corr
  USE params, only:ndim
  USE stoch_params, only:mode
  USE util, only: mat_trace, mat_contract

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: init_MTV_int_tensor

  ! Constant vectors
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, PUBLIC :: H1   !< First constant vector
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, PUBLIC :: H2   !< Second constant vector
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, PUBLIC :: H3   !< Third constant vector
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, PUBLIC :: Htot !< Total constant vector

  ! Tensors for the linear terms
  TYPE(coolist), DIMENSION(:), ALLOCATABLE, PUBLIC :: L1   !< First linear tensor
  TYPE(coolist), DIMENSION(:), ALLOCATABLE, PUBLIC :: L2   !< Second linear tensor
  TYPE(coolist), DIMENSION(:), ALLOCATABLE, PUBLIC :: L3   !< Third linear tensor
  TYPE(coolist), DIMENSION(:), ALLOCATABLE, PUBLIC :: Ltot !< Total linear tensor

  ! Tensors for the quadratic terms
  TYPE(coolist), DIMENSION(:), ALLOCATABLE, PUBLIC :: B1 !< First quadratic tensor
  TYPE(coolist), DIMENSION(:), ALLOCATABLE, PUBLIC :: B2 !< Second quadratic tensor
  TYPE(coolist), DIMENSION(:), ALLOCATABLE, PUBLIC :: Btot !< Total quadratic tensor

  TYPE(coolist4), DIMENSION(:), ALLOCATABLE, PUBLIC :: Mtot !< Tensor for the cubic terms

  ! Tensors for the noise covariance matrices
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: Q1 !< Constant terms for the state-dependent noise covariance matrix
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: Q2 !< Constant terms for the state-independent noise covariance matrix
  TYPE(coolist), DIMENSION(:), ALLOCATABLE, PUBLIC :: Utot !< Linear terms for the state-dependent noise covariance matrix
  TYPE(coolist4), DIMENSION(:), ALLOCATABLE, PUBLIC :: Vtot !< Quadratic terms for the state-dependent noise covariance matrix

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: dumb_vec !< Dummy vector
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dumb_mat1 !< Dummy matrix
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dumb_mat2 !< Dummy matrix
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dumb_mat3 !< Dummy matrix
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dumb_mat4 !< Dummy matrix


  !-----------------------------------------------------!
  !                                                     !
  ! End of preamble                                     !
  !                                                     !
  !-----------------------------------------------------!

CONTAINS

  !-----------------------------------------------------!
  !                                                     !
  ! Initialisation routine                              !
  !                                                     !
  !-----------------------------------------------------!

  !> Subroutine to initialise the MTV tensor
  SUBROUTINE init_MTV_int_tensor
    INTEGER :: AllocStat,i,j,k,l

    print*, 'Initializing the decomposition tensors...'
    CALL init_dec_tensor
    print*, "Initializing the correlation matrices and tensors..."
    CALL init_corrint
    print*, "Computing the correlation integrated matrices and tensors..."
    CALL comp_corrint
    
    !H part
    print*, "Computing the H term..."

    ALLOCATE(H1(0:ndim), H2(0:ndim), H3(0:ndim), Htot(0:ndim), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    ALLOCATE(dumb_mat1(ndim,ndim), dumb_mat2(ndim,ndim), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    ALLOCATE(dumb_mat3(ndim,ndim), dumb_mat4(ndim,ndim), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    !H1
    CALL coo_to_mat_ik(Lxy,dumb_mat1)
    dumb_mat2=matmul(dumb_mat1,corrint)
    CALL sparse_mul3_with_mat(Bxxy,dumb_mat2,H1)
    
    ! H2
    H2=0.D0
    IF (mode.ne.'ures') THEN
       CALL coo_to_mat_ik(Lyy,dumb_mat1)
       dumb_mat1=matmul(inv_corr_i_full,dumb_mat1)
    
       DO i=1,ndim
          CALL coo_to_mat_i(i,Bxyy,dumb_mat2)
          CALL sparse_mul4_with_mat_jl(corr2int,dumb_mat2,dumb_mat3)
          CALL sparse_mul4_with_mat_jl(corr2int,transpose(dumb_mat2),dumb_mat4)
          dumb_mat3=dumb_mat3+dumb_mat4
          H2(i)=mat_contract(dumb_mat1,dumb_mat3)
       ENDDO
    ENDIF

    !H3
    H3=0.D0
    CALL sparse_mul3_with_mat(Bxyy,corr_i_full,H3)

    !Htot
    Htot=0.D0
    Htot=H1+H2+H3

    print*, "Computing the L terms..."
    ALLOCATE(L1(ndim), L2(ndim), L3(ndim), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    !L1
    CALL coo_to_mat_ik(Lyx,dumb_mat1)
    CALL coo_to_mat_ik(Lxy,dumb_mat2)
    dumb_mat3=matmul(inv_corr_i_full,corrint)
    dumb_mat4=matmul(dumb_mat2,matmul(transpose(dumb_mat3),dumb_mat1))
    CALL matc_to_coo(dumb_mat4,L1)
    
    !L2
    dumb_mat4=0.D0
    DO i=1,ndim
       DO j=1,ndim
          CALL coo_to_mat_i(i,Bxyy,dumb_mat1)
          CALL sparse_mul4_with_mat_jl(corr2int,dumb_mat1+transpose(dumb_mat1),dumb_mat2)
          
          CALL coo_to_mat_j(j,Byxy,dumb_mat1)
          dumb_mat1=matmul(inv_corr_i_full,dumb_mat1)
          dumb_mat4(i,j)=mat_contract(dumb_mat1,dumb_mat2)
       END DO
    END DO
    CALL matc_to_coo(dumb_mat4,L2)

    !L3
    dumb_mat4=0.D0
    DO i=1,ndim
       DO j=1,ndim
          CALL coo_to_mat_j(j,Bxxy,dumb_mat1)
          CALL coo_to_mat_i(i,Bxxy,dumb_mat2)
          dumb_mat4(i,j)=mat_trace(matmul(dumb_mat1,matmul(corrint,transpose(dumb_mat2))))
       ENDDO
    END DO
    CALL matc_to_coo(dumb_mat4,L3)

    !Ltot

    ALLOCATE(Ltot(ndim), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
    
    CALL add_to_tensor(L1,Ltot)
    CALL add_to_tensor(L2,Ltot)
    CALL add_to_tensor(L3,Ltot)

    print*, "Computing the B terms..."
    ALLOCATE(B1(ndim), B2(ndim), Btot(ndim), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
    ALLOCATE(dumb_vec(ndim), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    ! B1
    CALL coo_to_mat_ik(Lxy,dumb_mat1)
    dumb_mat2=matmul(inv_corr_i_full,corrint)
    
    dumb_mat3=matmul(dumb_mat1,transpose(dumb_mat2))
    DO j=1,ndim
       DO k=1,ndim
          CALL coo_to_vec_jk(j,k,Byxx,dumb_vec)
          dumb_vec=matmul(dumb_mat3,dumb_vec)
          CALL add_vec_jk_to_tensor(j,k,dumb_vec,B1)
       ENDDO
    END DO

    ! B2
    CALL coo_to_mat_ik(Lyx,dumb_mat3)
    dumb_mat2=matmul(inv_corr_i_full,corrint)

    dumb_mat4=matmul(transpose(dumb_mat2),dumb_mat3)
    DO i=1,ndim
       CALL coo_to_mat_i(i,Bxxy,dumb_mat1)
       dumb_mat2=matmul(dumb_mat1,dumb_mat4)
       CALL add_matc_to_tensor(i,dumb_mat2,B2)
    ENDDO

    CALL add_to_tensor(B1,Btot)
    CALL add_to_tensor(B2,Btot)

    !M

    print*, "Computing the M term..."

    ALLOCATE(Mtot(ndim), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    dumb_mat2=matmul(inv_corr_i_full,corrint)

    DO i=1,ndim
       CALL coo_to_mat_i(i,Bxxy,dumb_mat1)
       dumb_mat3=matmul(dumb_mat1,transpose(dumb_mat2))
       DO k=1,ndim
          DO l=1,ndim
             CALL coo_to_vec_jk(k,l,Byxx,dumb_vec)
             dumb_vec=matmul(dumb_mat3,dumb_vec)
             CALL add_vec_ikl_to_tensor4(i,k,l,dumb_vec,Mtot)
          ENDDO
       END DO
    END DO

    !Q

    print*, "Computing the Q terms..."
    ALLOCATE(Q1(ndim,ndim), Q2(ndim,ndim), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    !Q1

    CALL coo_to_mat_ik(Lxy,dumb_mat1)
    Q1=matmul(dumb_mat1,matmul(corrint,transpose(dumb_mat1)))
    
    !Q2

    DO i=1,ndim
       DO j=1,ndim
          CALL coo_to_mat_i(i,Bxyy,dumb_mat1)
          CALL coo_to_mat_i(j,Bxyy,dumb_mat2)
          CALL sparse_mul4_with_mat_jl(corr2int,dumb_mat2,dumb_mat3)
          CALL sparse_mul4_with_mat_jl(corr2int,transpose(dumb_mat2),dumb_mat4)
          dumb_mat2=dumb_mat3+dumb_mat4
          Q2(i,j)=mat_contract(dumb_mat1,dumb_mat2)
       END DO
    END DO
          
    !U
    
    ALLOCATE(Utot(ndim), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    CALL coo_to_mat_ik(Lxy,dumb_mat1)
    DO i=1,ndim
       CALL coo_to_mat_i(i,Bxxy,dumb_mat2)
       dumb_mat3=matmul(dumb_mat1,matmul(corrint,transpose(dumb_mat2)))
       CALL add_matc_to_tensor(i,dumb_mat3,Utot)
    ENDDO

    DO j=1,ndim
       CALL coo_to_mat_i(j,Bxxy,dumb_mat2)
       dumb_mat3=matmul(dumb_mat1,matmul(corrint,transpose(dumb_mat2)))
       DO k=1,ndim
          CALL add_vec_jk_to_tensor(j,k,dumb_mat3(:,k),Utot)
       ENDDO
    ENDDO
    
    !V
    
    ALLOCATE(Vtot(ndim), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    DO i=1,ndim
       DO j=1,ndim
          CALL coo_to_mat_i(i,Bxxy,dumb_mat1)
          CALL coo_to_mat_i(j,Bxxy,dumb_mat2)
          dumb_mat3=matmul(dumb_mat1,matmul(corrint,transpose(dumb_mat2)))
          CALL add_matc_to_tensor4(j,i,dumb_mat3,Vtot)
       ENDDO
    ENDDO

    DEALLOCATE(dumb_mat1, dumb_mat2, STAT=AllocStat)
    IF (AllocStat /= 0)  STOP "*** Problem to deallocate ! ***"

    DEALLOCATE(dumb_mat3, dumb_mat4, STAT=AllocStat)
    IF (AllocStat /= 0)  STOP "*** Problem to deallocate ! ***"

    DEALLOCATE(dumb_vec, STAT=AllocStat)
    IF (AllocStat /= 0)  STOP "*** Problem to deallocate ! ***"

    
  END SUBROUTINE init_MTV_int_tensor

  
END MODULE MTV_int_tensor



