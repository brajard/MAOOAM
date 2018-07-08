
! test_WL_tensor.f90
!
!>  Small program to print the WL tensors.
!     
!> @copyright                                                               
!> 2018 Jonathan Demaeyer.
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------!
PROGRAM test_WL_tensor

  USE aotensor_def, only: init_aotensor
  USE WL_tensor
  USE tensor, only: print_tensor,print_tensor4
  USE stoch_params, only:mems

  INTEGER :: m

  CALL init_aotensor     ! Initialize the aotensor 
  CALL init_WL_tensor    ! Compute the tensor

  print*, 'M12def',M12def
  print*, 'M21def',M21def
  print*, 'M22def',M22def
  print*, 'Ldef',Ldef
  print*, 'B14def',B14def
  print*, 'B23def',B23def
  print*, 'Mdef',Mdef
  

  print*, 'M11'
  print*, M11
  print*, 'M12'
  CALL print_tensor(M12)
  print*, 'M13'
  print*, M13
  print*, 'M1tot'
  print*, M1tot
  print*, 'M21'
  CALL print_tensor(M21)
  print*, 'M22'
  CALL print_tensor(M22)
  print*, "M3 terms..."
  print*, 'L1'
  DO m=1,mems
     print*, m
     CALL print_tensor(L1(:,m))
  ENDDO
  print*, 'L2'
  DO m=1,mems
     print*, m
     CALL print_tensor(L2(:,m))
  ENDDO
  print*, 'L4'
  DO m=1,mems
     print*, m
     CALL print_tensor(L4(:,m))
  ENDDO
  print*, 'L5'
  DO m=1,mems
     print*, m
     CALL print_tensor(L5(:,m))
  ENDDO
  print*, 'Ltot'
  DO m=1,mems
     print*, m
     CALL print_tensor(Ltot(:,m))
  ENDDO
    print*, 'B1'
  DO m=1,mems
     print*, m
     CALL print_tensor(B1(:,m))
  ENDDO
  print*, 'B2'
  DO m=1,mems
     print*, m
     CALL print_tensor(B2(:,m))
  ENDDO
  print*, 'B3'
  DO m=1,mems
     print*, m
     CALL print_tensor(B3(:,m))
  ENDDO
  print*, 'B4'
  DO m=1,mems
     print*, m
     CALL print_tensor(B3(:,m))
  ENDDO
  print*, 'B14'
  DO m=1,mems
     print*, m
     CALL print_tensor(B14(:,m))
  ENDDO
  print*, 'B23'
  DO m=1,mems
     print*, m
     CALL print_tensor(B23(:,m))
  ENDDO
  print*, 'Mtot'
  DO m=1,mems
     print*, m
     CALL print_tensor4(Mtot(:,m))
  ENDDO


END PROGRAM test_WL_tensor

