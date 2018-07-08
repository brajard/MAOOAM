
! corr_tensor.f90
!
!> Module to compute the correlations and derivatives used to compute
!> the memory term of the WL parameterization
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

MODULE corr_tensor
  USE tensor
  USE params, only: ndim
  USE stoch_params, only:mems,muti
  USE corrmod, only:init_corr,corr_ij,corrcomp,inv_corr_i,mean
  USE util, only:ireduce,vector_outer
  USE sf_def, only: n_unres,ind,rind
  IMPLICIT NONE

  PRIVATE
    
  PUBLIC :: init_corr_tensor

  TYPE(coolist), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: YY !< Coolist holding the  \f$\langle Y \otimes Y^s \rangle\f$ terms
  TYPE(coolist), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: dY !< Coolist holding the  \f$\langle \partial_Y \otimes Y^s \rangle\f$ terms
  TYPE(coolist), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: YdY !< Coolist holding the  \f$\langle Y \otimes \partial_Y \otimes Y^s \rangle\f$ terms
  TYPE(coolist), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: dYY !< Coolist holding the  \f$\langle \partial_Y \otimes Y^s \otimes Y^s \rangle\f$ terms
  TYPE(coolist4), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: YdYY !< Coolist holding the  \f$\langle Y \otimes \partial_Y \otimes Y^s \otimes Y^s \rangle\f$ terms

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: dumb_vec !< Dumb vector to be used in the calculation
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dumb_mat1 !< Dumb matrix to be used in the calculation
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dumb_mat2 !< Dumb matrix to be used in the calculation
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: expm !< Matrix holding the product inv_corr_i*corr_ij at time \f$s\f$
    
CONTAINS
  
  !> Subroutine to initialise the correlations tensors
  SUBROUTINE init_corr_tensor
    INTEGER :: i,j,m,AllocStat

    CALL init_corr

    print*, 'Computing the time correlation tensors...'
    
    ALLOCATE(YY(ndim,mems),dY(ndim,mems), dYY(ndim,mems), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    ALLOCATE(YdY(ndim,mems), YdYY(ndim,mems), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    ALLOCATE(dumb_vec(ndim), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    ALLOCATE(dumb_mat1(ndim,ndim), dumb_mat2(ndim,ndim), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    ALLOCATE(expm(n_unres,n_unres), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    DO m=1,mems
       CALL corrcomp((m-1)*muti)
       
       ! YY
       CALL ireduce(dumb_mat2,corr_ij,n_unres,ind,rind)
       CALL matc_to_coo(dumb_mat2,YY(:,m))

       ! dY
       expm=matmul(inv_corr_i,corr_ij)
       CALL ireduce(dumb_mat2,expm,n_unres,ind,rind)
       CALL matc_to_coo(dumb_mat2,dY(:,m))

       ! YdY
       DO i=1,n_unres
          CALL ireduce(dumb_mat2,mean(i)*expm,n_unres,ind,rind)
          CALL add_matc_to_tensor(ind(i),dumb_mat2,YdY(:,m))
       ENDDO

       ! dYY
       dumb_vec(1:n_unres)=matmul(mean,expm)
       DO i=1,n_unres
          CALL vector_outer(expm(i,:),dumb_vec(1:n_unres),dumb_mat2(1:n_unres,1:n_unres))
          CALL ireduce(dumb_mat1,dumb_mat2+transpose(dumb_mat2),n_unres,ind,rind)
          CALL add_matc_to_tensor(ind(i),dumb_mat1,dYY(:,m))
       ENDDO

       ! YdYY
       DO i=1,n_unres
          DO j=1,n_unres
             CALL vector_outer(corr_ij(i,:),expm(j,:),dumb_mat2(1:n_unres,1:n_unres))
             CALL ireduce(dumb_mat1,dumb_mat2+transpose(dumb_mat2),n_unres,ind,rind)
             CALL add_matc_to_tensor4(ind(i),ind(j),dumb_mat1,YdYY(:,m))
          ENDDO
       ENDDO
    ENDDO

    DEALLOCATE(dumb_mat1, dumb_mat2, STAT=AllocStat)
    IF (AllocStat /= 0)  STOP "*** Problem to deallocate ! ***"

    DEALLOCATE(dumb_vec, STAT=AllocStat)
    IF (AllocStat /= 0)  STOP "*** Problem to deallocate ! ***"


  END SUBROUTINE init_corr_tensor

END MODULE corr_tensor

       
       
      
