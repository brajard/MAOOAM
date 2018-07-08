
! test_MTV_int_tensor.f90
!
!>  Small program to print the time correlations tensors.
!     
!> @copyright                                                               
!> 2018 Jonathan Demaeyer.
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------!

PROGRAM test_corr_tensor

  USE tensor, only: print_tensor,print_tensor4
  USE aotensor_def, only: init_aotensor
  USE dec_tensor, only: init_dec_tensor
  USE corr_tensor
  USE stoch_params, only: mems

  INTEGER :: m

  CALL init_aotensor     ! Initialize the aotensor 
  CALL init_dec_tensor

  CALL init_corr_tensor

  print*, 'YY'
  DO m=1,mems
     print*, m
     CALL print_tensor(YY(:,m))
  ENDDO
  print*, 'dY'
  DO m=1,mems
     print*, m
     CALL print_tensor(dY(:,m))
  ENDDO
  print*, 'YdY'
  DO m=1,mems
     print*, m
     CALL print_tensor(YdY(:,m))
  ENDDO
  print*, 'dYY'
  DO m=1,mems
     print*, m
     CALL print_tensor(dYY(:,m))
  ENDDO
  print*, 'YdYY'
  DO m=1,mems
     print*, m
     CALL print_tensor4(YdYY(:,m))
  ENDDO


END PROGRAM test_corr_tensor

