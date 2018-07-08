
! WL_tensor.f90 
!
!>  The WL tensors used to integrate the model
!
!> @copyright                                                               
!> 2018 Jonathan Demaeyer.
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------!
!                                                                           
!> @remark                                                                 
!>   
!> 
!                                                                           
!---------------------------------------------------------------------------!

MODULE WL_tensor

  !-----------------------------------------------------!
  !                                                     !
  ! Preamble and variables declaration                  !
  !                                                     !
  !-----------------------------------------------------!

  USE tensor
  USE dec_tensor
  USE params, only:ndim
  USE stoch_params, only:mems
  USE util, only: mat_trace, mat_contract
  USE corr_tensor
  USE corrmod, only: corr_i_full,mean_full


  IMPLICIT NONE

  PRIVATE

  PUBLIC :: init_WL_tensor

  ! M1 term vectors
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, PUBLIC :: M11    !< First component of the M1 term
  TYPE(coolist), DIMENSION(:), ALLOCATABLE, PUBLIC :: M12   !< Second component of the M1 term
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, PUBLIC :: M13    !< Third component of the M1 term
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, PUBLIC :: M1tot  !< Total \f$M_1\f$ vector 

  ! M2 term tensors
  TYPE(coolist), DIMENSION(:), ALLOCATABLE, PUBLIC :: M21   !< First tensor of the M2 term 
  TYPE(coolist), DIMENSION(:), ALLOCATABLE, PUBLIC :: M22   !< Second tensor of the M2 term

  ! M3 terms tensors 
  ! Tensors for the linear terms (no L3 for WL)
  TYPE(coolist), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: L1    !< First linear tensor  
  TYPE(coolist), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: L2    !< Second linear tensor
  TYPE(coolist), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: L4    !< Fourth linear tensor
  TYPE(coolist), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: L5    !< Fifth linear tensor
  TYPE(coolist), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: Ltot  !< Total linear tensor

  ! Tensors for the quadratic terms
  TYPE(coolist), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: B1    !< First quadratic tensor
  TYPE(coolist), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: B2    !< Second quadratic tensor
  TYPE(coolist), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: B3    !< Third quadratic tensor
  TYPE(coolist), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: B4    !< Fourth quadratic tensor
  TYPE(coolist), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: B14   !< Joint 1st and 4th tensors
  TYPE(coolist), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: B23   !< Joint 2nd and 3rd tensors

  TYPE(coolist4), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: Mtot !< Tensor for the cubic terms

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: dumb_vec !< Dummy vector
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dumb_mat1 !< Dummy matrix
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dumb_mat2 !< Dummy matrix
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dumb_mat3 !< Dummy matrix
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dumb_mat4 !< Dummy matrix

  LOGICAL, PUBLIC :: M12def,M21def,M22def,Ldef,B14def,B23def,Mdef !< Boolean to (de)activate the computation of the terms

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


  !> Subroutine to initialise the WL tensor
  SUBROUTINE init_WL_tensor
    INTEGER :: AllocStat,i,j,k,m
    
    print*, 'Initializing the decompostion tensors...'
    CALL init_dec_tensor
    print*, "Initializing the correlation matrices and tensors..."
    CALL init_corr_tensor
    
    !M1 part
    print*, "Computing the M1 terms..."

    ALLOCATE(M11(0:ndim), M12(ndim), M13(0:ndim), M1tot(0:ndim), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    ALLOCATE(dumb_mat1(ndim,ndim), dumb_mat2(ndim,ndim), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    ALLOCATE(dumb_mat3(ndim,ndim), dumb_mat4(ndim,ndim), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    ALLOCATE(dumb_vec(ndim), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    !M11
    M11=0.D0
    ! CALL coo_to_mat_ik(Lxy,dumb_mat1)
    ! M11(1:ndim)=matmul(dumb_mat1,mean_full(1:ndim))

    !M12
    ! dumb_mat2=0.D0
    ! DO i=1,ndim
    !    CALL coo_to_mat_i(i,Bxxy,dumb_mat1)
    !    dumb_mat2(i,:)=matmul(dumb_mat1,mean_full(1:ndim))
    ! ENDDO
    ! CALL matc_to_coo(dumb_mat2,M12)

    M12def=.not.tensor_empty(M12)

    !M13
    M13=0.D0
    CALL sparse_mul3_with_mat(Bxyy,corr_i_full,M13)

    !M1tot
    M1tot=0.D0
    M1tot=M11+M13

    print*, "Computing the M2 terms..."
    ALLOCATE(M21(ndim), M22(ndim), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    !M21
    CALL copy_coo(Lxy,M21)
    CALL add_to_tensor(Bxxy,M21)

    M21def=.not.tensor_empty(M21)
    
    !M22
    CALL copy_coo(Bxyy,M22)

    M22def=.not.tensor_empty(M22)
    
    !M3 tensor
    print*, "Computing the M3 terms..."
    ! Linear terms
    print*, "Computing the L subterms..."
    ALLOCATE(L1(ndim,mems), L2(ndim,mems), L4(ndim,mems), L5(ndim,mems), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
    
    !L1
    CALL coo_to_mat_ik(Lyx,dumb_mat1)
    CALL coo_to_mat_ik(Lxy,dumb_mat2)
    DO m=1,mems
       CALL coo_to_mat_ik(dY(:,m),dumb_mat3)
       dumb_mat4=matmul(dumb_mat2,matmul(transpose(dumb_mat3),dumb_mat1))
       CALL matc_to_coo(dumb_mat4,L1(:,m))
    ENDDO

    !L2
    DO m=1,mems
       dumb_mat4=0.D0
       DO i=1,ndim
          CALL coo_to_mat_i(i,Bxyy,dumb_mat1)
          CALL sparse_mul4_with_mat_kl(YdYY(:,m),dumb_mat1,dumb_mat2)
          DO j=1,ndim
             CALL coo_to_mat_j(j,Byxy,dumb_mat1)
             dumb_mat4(i,j)=mat_trace(matmul(dumb_mat1,dumb_mat2))
          ENDDO
       END DO
       CALL matc_to_coo(dumb_mat4,L2(:,m))
    ENDDO

    !L4
    ! DO m=1,mems
    !    dumb_mat4=0.D0
    !    DO i=1,ndim
    !       CALL coo_to_mat_i(i,Bxyy,dumb_mat1)
    !       CALL sparse_mul3_with_mat(dYY(:,m),dumb_mat1,dumb_vec) ! Bxyy*dYY
    !       CALL coo_to_mat_ik(Lyx,dumb_mat1)
    !       dumb_mat4(i,:)=matmul(transpose(dumb_mat1),dumb_vec)
    !    ENDDO
    !    CALL matc_to_coo(dumb_mat4,L4(:,m))           
    ! ENDDO

    !L5
    
    ! CALL coo_to_mat_ik(Lxy,dumb_mat1)
    ! DO m=1,mems
    !    dumb_mat4=0.D0
    !    DO i=1,ndim
    !       CALL sparse_mul3_mat(YdY(:,m),dumb_mat1(i,:),dumb_mat2)
    !       DO j=1,ndim
    !          CALL coo_to_mat_j(j,Byxy,dumb_mat3)
    !          dumb_mat4(i,j)=mat_trace(matmul(dumb_mat3,dumb_mat2))
    !       ENDDO
    !    END DO
    !    CALL matc_to_coo(dumb_mat4,L5(:,m))
    ! ENDDO

    !Ltot

    ALLOCATE(Ltot(ndim,mems), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
    
    DO m=1,mems
       CALL add_to_tensor(L1(:,m),Ltot(:,m))
       CALL add_to_tensor(L2(:,m),Ltot(:,m))
       CALL add_to_tensor(L4(:,m),Ltot(:,m))
       CALL add_to_tensor(L5(:,m),Ltot(:,m))
    ENDDO

    Ldef=.not.tensor_empty(Ltot)
       
    print*, "Computing the B terms..."
    ALLOCATE(B1(ndim,mems), B2(ndim,mems), B3(ndim,mems), B4(ndim,mems), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
    
    ! B1
    CALL coo_to_mat_ik(Lxy,dumb_mat1)
    dumb_mat1=transpose(dumb_mat1)
    DO m=1,mems
       CALL coo_to_mat_ik(dY(:,m),dumb_mat2)
       dumb_mat2=matmul(dumb_mat2,dumb_mat1)
       DO j=1,ndim
          DO k=1,ndim
             CALL coo_to_vec_jk(j,k,Byxx,dumb_vec)
             dumb_vec=matmul(dumb_vec,dumb_mat2)
             CALL add_vec_jk_to_tensor(j,k,dumb_vec,B1(:,m))
          ENDDO
       ENDDO
    ENDDO
    
    ! B2
    CALL coo_to_mat_ik(Lyx,dumb_mat3)
    dumb_mat3=transpose(dumb_mat3)
    DO m=1,mems
       DO i=1,ndim
          CALL coo_to_mat_i(i,Bxxy,dumb_mat1)
          CALL coo_to_mat_ik(dY(:,m),dumb_mat2)
          dumb_mat1=matmul(dumb_mat2,transpose(dumb_mat1))
          dumb_mat1=matmul(dumb_mat3,dumb_mat1)
          CALL add_matc_to_tensor(i,dumb_mat1,B2(:,m))
       ENDDO
    ENDDO

    ! B3
    ! DO m=1,mems
    !    DO i=1,ndim
    !       CALL coo_to_mat_i(i,Bxxy,dumb_mat1)
    !       dumb_mat4=0.D0
    !       DO j=1,ndim
    !          CALL coo_to_mat_j(j,YdY(:,m),dumb_mat2)
    !          CALL coo_to_mat_i(j,Byxy,dumb_mat3)
    !          dumb_mat2=matmul(dumb_mat3,dumb_mat2)
    !          dumb_mat4=dumb_mat4+dumb_mat2
    !       ENDDO
    !       dumb_mat4=matmul(dumb_mat4,transpose(dumb_mat1))
    !       CALL add_matc_to_tensor(i,dumb_mat4,B3(:,m))
    !    ENDDO
    ! ENDDO

    ! B4
    ! DO m=1,mems
    !    DO i=1,ndim
    !       CALL coo_to_mat_i(i,Bxyy,dumb_mat1)
    !       CALL sparse_mul3_with_mat(dYY(:,m),dumb_mat1,dumb_vec) ! Bxyy*dYY
    !       DO j=1,ndim
    !          CALL coo_to_mat_j(j,Byxx,dumb_mat1)
    !          dumb_mat4(j,:)=matmul(transpose(dumb_mat1),dumb_vec)
    !       ENDDO
    !       CALL add_matc_to_tensor(i,dumb_mat4,B4(:,m))
    !    ENDDO
    ! ENDDO

    ALLOCATE(B14(ndim,mems), B23(ndim,mems), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    DO m=1,mems
       CALL add_to_tensor(B1(:,m),B14(:,m))
       CALL add_to_tensor(B2(:,m),B23(:,m))
       CALL add_to_tensor(B4(:,m),B14(:,m))
       CALL add_to_tensor(B3(:,m),B23(:,m))
    ENDDO

    B14def=.not.tensor_empty(B14)
    B23def=.not.tensor_empty(B23)

    !M

    print*, "Computing the M term..."

    ALLOCATE(Mtot(ndim,mems), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    DO m=1,mems
       DO i=1,ndim
          CALL coo_to_mat_i(i,Bxxy,dumb_mat1)
          CALL coo_to_mat_ik(dY(:,m),dumb_mat2)
          dumb_mat1=matmul(dumb_mat2,transpose(dumb_mat1))
          DO j=1,ndim
             DO k=1,ndim
                CALL coo_to_vec_jk(j,k,Byxx,dumb_vec)
                dumb_vec=matmul(dumb_vec,dumb_mat1)
                CALL add_vec_ijk_to_tensor4(i,j,k,dumb_vec,Mtot(:,m))
             ENDDO
          END DO
       END DO
    END DO

    Mdef=.not.tensor4_empty(Mtot)
    

    DEALLOCATE(dumb_mat1, dumb_mat2, STAT=AllocStat)
    IF (AllocStat /= 0)  STOP "*** Problem to deallocate ! ***"

    DEALLOCATE(dumb_mat3, dumb_mat4, STAT=AllocStat)
    IF (AllocStat /= 0)  STOP "*** Problem to deallocate ! ***"

    DEALLOCATE(dumb_vec, STAT=AllocStat)
    IF (AllocStat /= 0)  STOP "*** Problem to deallocate ! ***"

    
  END SUBROUTINE init_WL_tensor

  
END MODULE WL_tensor



