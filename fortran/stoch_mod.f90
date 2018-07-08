
! stoch_mod.f90
!
!> Utility module containing the stochastic related routines.
!
!> @copyright                                                               
!> 2018 Jonathan Demaeyer.
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------!
!                                                                           
!> @remark                                                                 
!
!---------------------------------------------------------------------------!
MODULE stoch_mod

  USE params, only: ndim,natm
  USE sf_def, only: SF 
  IMPLICIT NONE

  PRIVATE

  INTEGER :: iset=0
  REAL(KIND=8) :: gset
  
  PUBLIC :: gasdev,stoch_vec,stoch_atm_res_vec,stoch_atm_unres_vec
  PUBLIC :: stoch_atm_vec,stoch_oc_vec,stoch_oc_res_vec,stoch_oc_unres_vec
  
CONTAINS

  FUNCTION gasdev()
    REAL(KIND=8) :: gasdev
    REAL(KIND=8) :: fac,rsq,v1,v2,r
    IF (iset.eq.0) THEN
       DO
          CALL random_number(r)
          v1=2.D0*r-1.
          CALL random_number(r)
          v2=2.D0*r-1.
          rsq=v1**2+v2**2
          IF (rsq.lt.1.D0.and.rsq.ne.0.D0) EXIT
       ENDDO
       fac=sqrt(-2.*log(rsq)/rsq)
       gset=v1*fac
       gasdev=v2*fac
       iset=1
    ELSE
       gasdev=gset
       iset=0
    ENDIF
    RETURN
  END FUNCTION gasdev

  !> Routine to fill a vector with standard Gaussian noise process values 
  !> @param dW Vector to fill
  SUBROUTINE stoch_vec(dW)
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(INOUT) :: dW
    INTEGER :: i
    DO i=1,ndim
       dW(i)=gasdev()
    ENDDO
  END SUBROUTINE stoch_vec

  !> routine to fill the atmospheric component of a vector with standard gaussian noise process values
  !> @param dW vector to fill
  subroutine stoch_atm_vec(dW)
    real(kind=8), dimension(0:ndim), intent(inout) :: dW
    integer :: i
    do i=1,2*natm
       dW(i)=gasdev()
    enddo
  end subroutine stoch_atm_vec

  !> routine to fill the resolved atmospheric component of a vector with standard gaussian noise process values
  !> @param dW vector to fill
  subroutine stoch_atm_res_vec(dW)
    real(kind=8), dimension(0:ndim), intent(inout) :: dW
    integer :: i
    dW=0.D0
    do i=1,2*natm
       IF (SF(i)==0) dW(i)=gasdev()
    enddo
  end subroutine stoch_atm_res_vec

  !> routine to fill the unresolved atmospheric component of a vector with standard gaussian noise process values
  !> @param dW vector to fill
  subroutine stoch_atm_unres_vec(dW)
    real(kind=8), dimension(0:ndim), intent(inout) :: dW
    integer :: i
    dW=0.D0
    do i=1,2*natm
       IF (SF(i)==1) dW(i)=gasdev()
    enddo
  end subroutine stoch_atm_unres_vec

  !> routine to fill the oceanic component of a vector with standard gaussian noise process values
  !> @param dW vector to fill
  subroutine stoch_oc_vec(dW)
    real(kind=8), dimension(0:ndim), intent(inout) :: dW
    integer :: i
    do i=2*natm+1,ndim
       dW(i)=gasdev()
    enddo
  end subroutine stoch_oc_vec

  !> routine to fill the resolved oceanic component of a vector with standard gaussian noise process values
  !> @param dW vector to fill
  subroutine stoch_oc_res_vec(dW)
    real(kind=8), dimension(0:ndim), intent(inout) :: dW
    integer :: i
    dW=0.D0
    do i=2*natm+1,ndim
       IF (SF(i)==0) dW(i)=gasdev()
    enddo
  end subroutine stoch_oc_res_vec

  !> routine to fill the unresolved oceanic component of a vector with standard gaussian noise process values
  !> @param dW vector to fill
  subroutine stoch_oc_unres_vec(dW)
    real(kind=8), dimension(0:ndim), intent(inout) :: dW
    integer :: i
    dW=0.D0
    do i=2*natm+1,ndim
       IF (SF(i)==1) dW(i)=gasdev()
    enddo
  end subroutine stoch_oc_unres_vec

END MODULE stoch_mod
