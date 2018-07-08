
! test_MAR.f90
!
!> Small program to test the Multivariate AutoRegressive model
!     
!> @copyright                                                               
!> 2018 Jonathan Demaeyer.
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------!

PROGRAM test_MAR

  USE params, only: ndim,dt,tw
  USE aotensor_def, only: init_aotensor
  USE dec_tensor, only: init_dec_tensor
  USE sf_def, only: n_unres
  USE util, only: printmat
  USE MAR
    
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: x

  REAL(KIND=8) :: t
  INTEGER :: AllocStat,i,n

  CALL init_aotensor     ! Initialize the aotensor 
  CALL init_dec_tensor
  CALL init_MAR

  print*, 'Qred'
  CALL printmat(Qred)
  print*, 'Q'
  CALL printmat(Q)
  DO i=1,ms
     print*, 'Wred',i
     CALL printmat(Wred(i,:,:))
     print*, 'W',i
     CALL printmat(W(i,:,:))
  ENDDO
  
  print*, 'Generating 100 steps of MAR evolution ...'
  ALLOCATE(x(0:ndim,ms), STAT=AllocStat)
  x=0.D0
  DO i=1,100
     CALL MAR_step(x)
  ENDDO
  DO i=1,100
     CALL MAR_step(x)
     print*, i,x(1:ndim,1)
  ENDDO

  ! print*, 'Generating roughly 2500000 days of MAR reduced evolution and saving it to ./tests/test_MAR.dat ...'
  ! OPEN(30,file='./tests/test_MAR.dat')
  ! x=0.D0
  ! t=0.D0
  ! DO WHILE (t<5.D4)
  !    CALL MAR_step_red(x)
  !    t=t+dt
  ! ENDDO
  ! t=0.D0
  ! DO WHILE (t<2.5D7)
  !    CALL MAR_step_red(x)
  !    t=t+dt
  !    IF (mod(t,tw)<dt) WRITE(30,*) t,x(1:n_unres,1)
  ! ENDDO
  ! CLOSE(30)
  print*, 'Test finished'

END PROGRAM test_MAR

