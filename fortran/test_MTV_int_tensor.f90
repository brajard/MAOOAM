
! test_MTV_int_tensor.f90
!
!> Small program to print the MTV integrated tensors
!     
!> @copyright                                                               
!> 2018 Jonathan Demaeyer.
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------!

PROGRAM test_MTV_int_tensor

  USE aotensor_def, only: init_aotensor
  USE MTV_int_tensor
  USE tensor, only: print_tensor,print_tensor4
  USE sigma, only: init_sigma
  USE util, only: printmat

  LOGICAL :: mult,Q1fill
  
  CALL init_aotensor     ! Initialize the aotensor 
  CALL init_MTV_int_tensor
  print*, 'H1'
  print*, H1
  print*, 'H2'
  print*, H2
  print*, 'H3'
  print*, H3
  print*, 'Htot'
  print*, Htot
  print*, 'L1'
  CALL print_tensor(L1)
  print*, 'L2'
  CALL print_tensor(L2)
  print*, 'L3'
  CALL print_tensor(L3)
  print*, 'Ltot'
  CALL print_tensor(Ltot)
  print*, 'B1'
  CALL print_tensor(B1)
  print*, 'B2'
  CALL print_tensor(B2)
  print*, 'Btot'
  CALL print_tensor(Btot)
  print*, 'Mtot'
  CALL print_tensor4(Mtot)
  print*, 'Q1'
  CALL printmat(Q1)
  print*, 'Q2'
  CALL printmat(Q2)
  print*, 'Utot'
  CALL print_tensor(Utot)
  print*, 'Vtot'
  CALL print_tensor4(Vtot)

  CALL init_sigma(mult,Q1fill)
  

  print*, 'mult',mult
  print*, 'Q1fill',Q1fill

END PROGRAM test_MTV_int_tensor
