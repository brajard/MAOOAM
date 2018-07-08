
! test_MTV_sigma_tensor.f90
!
!> Small program to test the MTV noise sigma matrices
!     
!> @copyright                                                               
!> 2018 Jonathan Demaeyer.
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------!

PROGRAM test_sigma

  USE params, only: ndim
  USE sigma
  USE aotensor_def, only: init_aotensor
  USE MTV_int_tensor, only: init_MTV_int_tensor
  USE util, only: printmat,init_random_seed
  USE sqrt_mod, only: init_sqrt,sqrtm

  REAL(KIND=8), DIMENSION(5,5) :: A,B,sqB
  LOGICAL :: mult,Q1fill
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: y
  INTEGER :: info,info2

  ! A=0
  ! A(1,:)=(/5,-4,1,0,0/)
  ! A(2,:)=(/-4,6,-4,1,0/)
  ! A(3,:)=(/ 1,-4,6,-4,1/)
  ! A(4,:)=(/ 0 ,1,-4,6,-4/)
  ! A(5,:)=(/ 0 ,0,1,-4,6/)

  ! B=matmul(A,A)

  CALL init_aotensor     ! Initialize the aotensor 
  CALL init_MTV_int_tensor
  CALL init_sigma(mult,Q1fill)
  CALL init_random_seed
  ! CALL init_sqrt
  
  ! CALL sqrtm(B,sqB,info,info2,2)

  ! print*, info,info2

  ! print*, 'A'
  ! CALL printmat(A)


  ! print*, 'B'
  ! CALL printmat(B)

  ! print*, 'sqB'
  ! CALL printmat(sqB)
  ! print*, 'sqB*sqB'
  ! CALL printmat(matmul(sqB,sqB))

  print*, 'sig1'
  IF (mult) THEN
     ALLOCATE(y(0:ndim))
     DO i=1,100
        CALL random_number(y)
        y=2*y-1.D0
        y(0)=1.D0
        CALL compute_mult_sigma(y)
        PRINT*, i,y
        CALL printmat(sig1)
     ENDDO
  ELSE
     CALL printmat(sig1)
  ENDIF

  print*, 'sig2'
  CALL printmat(sig2)

  print*, 'mult',mult
  print*, 'Q1fill',Q1fill

END PROGRAM test_sigma
