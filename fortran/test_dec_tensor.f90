
! test_dec_tensor.f90
!
!> Small program to print the decomposed tensors.
!     
!> @copyright                                                               
!> 2018 Jonathan Demaeyer.
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------!

PROGRAM test_dec_tensor

  USE params, only: ndim
  USE dec_tensor
  USE tl_ad_tensor, only:tltensor
  USE aotensor_def, only:init_aotensor,aotensor
  USE tensor, only: print_tensor

  IMPLICIT NONE

  ! Program

  CALL init_aotensor     ! Initialize the aotensor 
  CALL init_dec_tensor    ! Compute the tensor

  print*, 'ao_tensor'
  CALL print_tensor(aotensor)

  print*, 'ss_tensor' 
  CALL print_tensor(ss_tensor)

  print*, 'sf_tensor' 
  CALL print_tensor(sf_tensor)

  print*, 'ff_tensor' 
  CALL print_tensor(ff_tensor)

  print*, 'fs_tensor' 
  CALL print_tensor(fs_tensor)

  print*, 'tltensor' 
  CALL print_tensor(tltensor)

  print*, 'ss_tl_tensor' 
  CALL print_tensor(ss_tl_tensor)

  print*, 'Hx'
  CALL print_tensor(Hx)

  print*, 'Lxx'
  CALL print_tensor(Lxx)

  print*, 'Lxy'
  CALL print_tensor(Lxy)

  print*, 'Bxxx'
  CALL print_tensor(Bxxx)

  print*, 'Bxxy'
  CALL print_tensor(Bxxy)

  print*, 'Bxyy'
  CALL print_tensor(Bxyy)

  print*, 'Hy'
  CALL print_tensor(Hy)

  print*, 'Lyx'
  CALL print_tensor(Lyx)

  print*, 'Lyy'
  CALL print_tensor(Lyy)

  print*, 'Byxx'
  CALL print_tensor(Byxx)

  print*, 'Byxy'
  CALL print_tensor(Byxy)

  print*, 'Byyy'
  CALL print_tensor(Byyy)

END PROGRAM test_dec_tensor

