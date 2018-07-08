!  maooam_WL.f90
!
!> Fortran 90 implementation of the modular arbitrary-order ocean-atmosphere 
!> model MAOOAM - WL parameterization. 
!
!> @copyright                                                               
!> 2018 Jonathan Demaeyer.
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------!

PROGRAM maooam_wl
  USE params, only: ndim, dt, tw, t_trans, t_run, writeout
  USE stoch_params, only: dtn, t_trans_stoch,muti,dts,dtsn,t_trans_mem
  USE ic_def, only: load_IC, IC
  USE rk2_WL_integrator, only: init_integrator,step,full_step
  USE aotensor_def, only: init_aotensor
  USE WL_tensor, only: init_WL_tensor
  USE memory, only: compute_M3
  USE stat
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: X       !< State variable in the model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: Xnew    !< Updated state variable
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: tend    !< Store the tendencies
  REAL(KIND=8) :: t=0.D0                             !< Time variable
  REAL(KIND=8) :: t_up                     !< Update time for the progress bar  

  PRINT*, 'Model MAOOAM v1.3 WL'
  PRINT*, 'Loading information...'

  CALL init_aotensor     ! Initialize the aotensor 
  CALL init_WL_tensor ! Load the tensors

  CALL load_IC          ! Load the initial condition

  CALL init_integrator  ! Initialize the integrator

  t_up=dt/t_trans*100.D0

  IF (writeout) OPEN(10,file='evol_WL.dat')
  IF (writeout) OPEN(12,file='ptend_WL.dat')
  IF (writeout) OPEN(13,file='mem_WL.dat')

  ALLOCATE(X(0:ndim),Xnew(0:ndim),tend(0:ndim))

  X=IC

  PRINT*, 'Starting an initial transient time full model evolution...'

  DO WHILE (t<t_trans)
     CALL full_step(X,t,dt,dtn,Xnew)
     X=Xnew
     IF (mod(t/t_trans*100.D0,0.1)<t_up) WRITE(*,'(" Progress ",F6.1," %",A,$)') t/t_trans*100.D0,char(13)
  ENDDO

  PRINT*, 'Starting a transient time to initialize memory ...'

  t=0.D0
  t_up=dt/t_trans_mem*100.D0

  DO WHILE (t<t_trans_mem)
     CALL full_step(X,t,dt,dtn,Xnew)
     X=Xnew
     IF (mod(t/t_trans_mem*100.D0,0.1)<t_up) WRITE(*,'(" Progress ",F6.1," %",A,$)') t/t_trans_mem*100.D0,char(13)
     IF (mod(t,muti)<dt) THEN
        CALL compute_M3(X,dts,dtsn,.true.,.true.,.true.,muti,Xnew)
        IF ((mod(t,max(tw,muti))<dt).and.(writeout)) WRITE(13,*) t,Xnew(1:ndim)
     ENDIF
  ENDDO
 
  CLOSE(13)

  PRINT*, 'Starting a transient time WL model evolution...'

  t=0.D0
  t_up=dt/t_trans_stoch*100.D0

  DO WHILE (t<t_trans_stoch)
     CALL step(X,t,dt,dtn,Xnew,tend)
     X=Xnew
     IF (mod(t/t_trans_stoch*100.D0,0.1)<t_up) WRITE(*,'(" Progress ",F6.1," %",A,$)') t/t_trans_stoch*100.D0,char(13)
  ENDDO
 
  PRINT*, 'Starting the WL time evolution...'

  CALL init_stat
  
  t=0.D0
  t_up=dt/t_run*100.D0

  DO WHILE (t<t_run)
     CALL step(X,t,dt,dtn,Xnew,tend)
     X=Xnew
     IF (mod(t,tw)<dt) THEN
        IF (writeout) WRITE(10,*) t,X(1:ndim)
        IF (writeout) WRITE(12,*) t,tend(1:ndim)
        CALL acc(X)
     END IF
     IF (mod(t/t_run*100.D0,0.1)<t_up) WRITE(*,'(" Progress ",F6.1," %",A,$)') t/t_run*100.D0,char(13)
  END DO

  PRINT*, 'Evolution finished.'

  IF (writeout) CLOSE(10)
  IF (writeout) CLOSE(12)

  IF (writeout) OPEN(10,file='mean_field_WL.dat')

  X=mean()
  IF (writeout) WRITE(10,*) X(1:ndim)
  IF (writeout) CLOSE(10)

END PROGRAM maooam_wl
