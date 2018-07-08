
! test_corr.f90
!
!> Small program to print the correlation and covariance matrices
!     
!> @copyright                                                               
!> 2018 Jonathan Demaeyer.
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------!

PROGRAM test_corr

  USE params, only: ndim
  USE aotensor_def, only: init_aotensor
  USE dec_tensor, only: init_dec_tensor
  USE int_corr
  USE corrmod
  USE tensor, only: print_tensor4
  USE util, only: printmat,reduce

  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dumb

  INTEGER, DIMENSION(:), ALLOCATABLE :: ind,rind
  INTEGER :: n
  
  CALL init_aotensor     ! Initialize the aotensor 
  CALL init_dec_tensor

  CALL init_corrint

  CALL comp_corrint

  ALLOCATE(dumb(ndim,ndim),ind(ndim),rind(ndim))
  
  PRINT*, 'full version'
  PRINT*, 'mean_full'
  PRINT*, mean_full

  PRINT*, 'corr_i_full'
  CALL printmat(corr_i_full)

  PRINT*, 'inv_corr_i_full'
  CALL printmat(inv_corr_i_full)

  PRINT*, 'inversion test'
  CALL printmat(matmul(corr_i_full,inv_corr_i_full))

  PRINT*, 'corrint full'
  CALL printmat(corrint)

  PRINT*, 'inv_corr_i_full*corrint full'
  CALL printmat(matmul(inv_corr_i_full,corrint))

  PRINT*, 'corr2int'
  CALL print_tensor4(corr2int)


  PRINT*, 'reduced part'
  PRINT*, 'mean'
  PRINT*, mean

  PRINT*, 'corr_i'
  CALL printmat(corr_i)

  PRINT*, 'inv_corr_i'
  CALL printmat(inv_corr_i)

  PRINT*, 'inversion test'
  CALL printmat(matmul(corr_i,inv_corr_i))

  PRINT*, 'corrint'
  CALL reduce(corrint,dumb,n,ind,rind)
  CALL printmat(dumb(1:n,1:n))

  PRINT*, 'inv_corr_i*corrint'
  CALL printmat(matmul(inv_corr_i,dumb(1:n,1:n)))

  PRINT*, 'corr2int'
  CALL print_tensor4(corr2int)


END PROGRAM test_corr
