
!  maooam_stoch.f90
!
!> Fortran 90 implementation of the stochastic modular arbitrary-order 
!> ocean-atmosphere model MAOOAM 
!>
!
!> @copyright                                                               
!> 2018 Jonathan Demaeyer.
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------!
!                                                                           
!>  @remark                                                                 
!>  There are four dynamics modes:
!>         - full: generate the full dynamics
!>         - unres: generate the intrinsic unresolved dynamics
!>         - qfst: generate dynamics given by the quadratic terms
!>                    of the unresolved tendencies
!>         - reso: use the resolved dynamics alone
!                                                                           
!---------------------------------------------------------------------------

PROGRAM maooam_stoch
  USE params, only: ndim, dt, tw, t_trans, t_run, writeout
  USE stoch_params, only: dtn,mode
  USE aotensor_def, only: init_aotensor
  USE dec_tensor, only: init_dec_tensor
  USE IC_def, only: load_IC, IC
  USE rk2_stoch_integrator, only: init_integrator,step
  USE stat
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: X       !< State variable in the model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: Xnew    !< Updated state variable
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: tend    !< Store the tendencies
  REAL(KIND=8) :: t=0.D0                             !< Time variable
  REAL(KIND=8) :: t_up                 !< Update time for the progress bar
  CHARACTER*4 :: force                 !< Selector for the dynamics
  INTEGER :: is,lg                     !< Dummy integers

  PRINT*, 'Model MAOOAM v1.3 stochastic'
  PRINT*, 'Loading information...'

  CALL get_command_argument(1,force,lg,is)

  CALL init_aotensor     ! Initialize the aotensor 
  CALL init_dec_tensor    ! Compute the tensor

  CALL load_IC          ! Load the initial condition

  ! Initialize the integrator
  IF (is.ne.0) THEN
     print*, 'No or bad dynamics specified, using the full one!'
     CALL init_integrator('full')
  ELSE
     IF (force.eq.'para') THEN
        CALL init_integrator
        print*, 'Using the stoch_params.nml specified dynamics: ', mode
     ELSE
        print*, 'Using the ', force, ' dynamics!'
        CALL init_integrator(force)
     ENDIF
  ENDIF
     
  t_up=dt/t_trans*100.D0

  IF (writeout) OPEN(10,file='evol_field.dat')

  ALLOCATE(X(0:ndim),Xnew(0:ndim),tend(0:ndim))

  X=IC

  PRINT*, 'Starting the transient time evolution...'

  DO WHILE (t<t_trans)
     CALL step(X,t,dt,dtn,Xnew,tend)
     X=Xnew
     IF (mod(t/t_trans*100.D0,0.1)<t_up) WRITE(*,'(" Progress ",F6.1," %",A,$)') t/t_trans*100.D0,char(13)
  END DO

  PRINT*, 'Starting the time evolution...'

  CALL init_stat
  
  t=0.D0
  t_up=dt/t_run*100.D0

  IF (writeout) WRITE(10,*) t,X(1:ndim)

  DO WHILE (t<t_run)
     CALL step(X,t,dt,dtn,Xnew,tend)
     X=Xnew
     IF (mod(t,tw)<dt) THEN
        IF (writeout) WRITE(10,*) t,X(1:ndim)
        CALL acc(X)
     END IF
     IF (mod(t/t_run*100.D0,0.1)<t_up) WRITE(*,'(" Progress ",F6.1," %",A,$)') t/t_run*100.D0,char(13)
  END DO

  PRINT*, 'Evolution finished.'

  IF (writeout) CLOSE(10)

  IF (writeout) OPEN(10,file='mean_field.dat')

  X=mean()
  IF (writeout) WRITE(10,*) X(1:ndim)
  IF (writeout) CLOSE(10)

END PROGRAM maooam_stoch
