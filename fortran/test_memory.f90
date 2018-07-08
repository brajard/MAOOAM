
! test_memory.f90
!
!> Small program to test the WL memory module
!     
!> @copyright                                                               
!> 2018 Jonathan Demaeyer.
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------!

PROGRAM test_memory
  USE params, only: ndim, dt
  USE stoch_params, only: dtn,muti,dts,dtsn,t_trans_stoch
  USE rk2_ss_integrator, only: init_ss_integrator,ss_step
  USE aotensor_def, only: init_aotensor
  USE WL_tensor, only: init_WL_tensor
  USE memory

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: y,ynew,h

  REAL(KIND=8) :: t=0.D0

  CALL init_aotensor     ! Initialize the aotensor 
  CALL init_WL_tensor    ! Compute the tensor
  CALL init_ss_integrator
  CALL init_memory

  print*, 'Computing 100 timunits of memory term...'

  ALLOCATE(y(0:ndim),ynew(0:ndim),h(0:ndim))
  
  y=0.D0
  y(0)=1.D0

  DO WHILE (t<t_trans_stoch)
     CALL ss_step(y,y,t,dt,dtn,ynew)
     y=ynew
     IF (mod(t,muti)<dt) CALL compute_M3(y,dts,dtsn,.true.,.true.,.true.,muti,h)
  ENDDO
  
  t=0.D0
  
  DO WHILE (t<100.)
     CALL ss_step(y,y,t,dt,dtn,ynew)
     y=ynew
     IF (mod(t,muti)<dt) THEN
        print*,'y',t,y
        CALL test_M3(y,dts,dtsn,h)
     ENDIF
  ENDDO

END PROGRAM test_memory
