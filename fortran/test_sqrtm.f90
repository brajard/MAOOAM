
! test_sqrtm.f90
!
!> Small program to test the matrix square-root module
!     
!> @copyright                                                               
!> 2018 Jonathan Demaeyer.
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------!

PROGRAM test_sqrtm

  USE params, only: init_params
  USE util, only: printmat
  USE sqrt_mod, only: init_sqrt,sqrtm,sqrtm_svd

  REAL(KIND=8), DIMENSION(3,3) :: A,B,sqB,sqA
  INTEGER :: info,info2

  A=0
  A(1,:)=(/9.1784812096335050E-009,5.7810098811039660E-009,3.8077584903965249E-009/)
  A(2,:)=(/5.7810098811039676E-009,3.6411334819038273E-009,2.3982932421037010E-009/)
  A(3,:)=(/3.8077584903965258E-009,2.3982932421037010E-009,1.5796758080159289E-009/)

  ! B=matmul(A,A)
  CALL init_params
  CALL init_sqrt

  print*, 'Schur decomposition'
  
  CALL sqrtm(A,sqA,info,info2,2)

  print*, info,info2

  print*, 'A'
  CALL printmat(A)

  print*, 'sqA'
  CALL printmat(sqA)

  print*, 'sqA*sqA'
  CALL printmat(matmul(sqA,sqA))

  print*, 'Svd Decomposition'

  CALL sqrtm_svd(A,sqA,info,info2,2)

  print*, 'A'
  CALL printmat(A)

  print*, 'sqA'
  CALL printmat(sqA)

  print*, 'sqA*sqA'
  CALL printmat(matmul(sqA,sqA))



  ! print*, 'B'
  ! CALL printmat(B)

  ! print*, 'sqB'
  ! CALL printmat(sqB)
  ! print*, 'sqB*sqB'
  ! CALL printmat(matmul(sqB,sqB))


END PROGRAM test_sqrtm
