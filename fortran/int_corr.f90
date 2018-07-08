
! int_corr.f90
!
!> Module to compute or load the integrals of the correlation matrices 
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

MODULE int_corr
  USE params, only: ndim
  USE sf_def, only: n_unres,ind
  USE tensor, only: coolist4,load_tensor4_from_file,write_tensor4_to_file
  USE int_comp, only: integrate
  USE corrmod, only: init_corr,corrcomp,corr_ij
  USE stoch_params, only: int_corr_mode

  PRIVATE

  INTEGER :: oi,oj,ok,ol !< Integers that specify the matrices and tensor component considered as a function of time
  REAL(KIND=8), PARAMETER :: real_eps = 2.2204460492503131e-16 !< Small epsilon constant to determine equality with zero

  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: corrint !< Matrix holding the integral of the correlation matrix
  TYPE(coolist4), DIMENSION(:), ALLOCATABLE, PUBLIC :: corr2int !< Tensor holding the integral of the correlation outer product with itself

  PUBLIC :: init_corrint,comp_corrint

CONTAINS

  !> Subroutine to initialise the integrated matrices and tensors
  SUBROUTINE init_corrint
    INTEGER :: AllocStat

    ALLOCATE(corrint(ndim,ndim), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    ALLOCATE(corr2int(ndim), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    CALL init_corr ! Initialize the correlation matrix function

    corrint=0.D0
        
  END SUBROUTINE init_corrint

  !> Function that returns the component oi and oj of the correlation matrix at time s
  !> @param s time at which the function is evaluated
  FUNCTION func_ij(s)
    IMPLICIT NONE
    REAL(KIND=8) :: s,func_ij
    CALL corrcomp(s)
    func_ij=corr_ij(oi,oj)
    RETURN
  END FUNCTION func_ij

  !> Function that returns the component oi,oj,ok and ol of the outer product of the
  !> correlation matrix with itself at time s
  !> @param s time at which the function is evaluated
  FUNCTION func_ijkl(s)
    IMPLICIT NONE
    REAL(KIND=8) :: s,func_ijkl
    CALL corrcomp(s)
    func_ijkl=corr_ij(oi,oj)*corr_ij(ok,ol)
    RETURN
  END FUNCTION func_ijkl

  !> Routine that actually compute or load the integrals
  SUBROUTINE comp_corrint
    IMPLICIT NONE
    INTEGER :: i,j,k,l,n,AllocStat
    REAL(KIND=8) :: ss
    LOGICAL :: ex

    INQUIRE(FILE='corrint.def',EXIST=ex)
    SELECT CASE (int_corr_mode)
    CASE ('file')
       IF (ex) THEN
          OPEN(30,file='corrint.def',status='old')
          READ(30,*) corrint
          CLOSE(30)
       ELSE
          STOP "*** File corrint.def not found ! ***"
       END IF
    CASE ('prog')
       DO i = 1,n_unres
          DO j= 1,n_unres
             oi=i
             oj=j
             CALL integrate(func_ij,ss)
             corrint(ind(i),ind(j))=ss
          END DO
       END DO

       OPEN(30,file='corrint.def')
       WRITE(30,*) corrint
       CLOSE(30)
    END SELECT


    INQUIRE(FILE='corr2int.def',EXIST=ex)
    SELECT CASE (int_corr_mode)
    CASE ('file')
       IF (ex) THEN
          CALL load_tensor4_from_file("corr2int.def",corr2int)
       ELSE
          STOP "*** File corr2int.def not found ! ***"
       END IF
    CASE ('prog')
       DO i = 1,n_unres
          n=0
          DO j= 1,n_unres
             DO k= 1,n_unres
                DO l = 1,n_unres
                   oi=i
                   oj=j
                   ok=k
                   ol=l

                   CALL integrate(func_ijkl,ss)
                   IF (abs(ss)>real_eps) n=n+1
                ENDDO
             ENDDO
          ENDDO
          IF (n/=0) THEN
             ALLOCATE(corr2int(ind(i))%elems(n), STAT=AllocStat)
             IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

             n=0
             DO j= 1,n_unres
                DO k= 1,n_unres
                   DO l = 1,n_unres
                      oi=i
                      oj=j
                      ok=k
                      ol=l

                      CALL integrate(func_ijkl,ss)
                      IF (abs(ss)>real_eps) THEN
                         n=n+1
                         corr2int(ind(i))%elems(n)%j=ind(j)
                         corr2int(ind(i))%elems(n)%k=ind(k)
                         corr2int(ind(i))%elems(n)%l=ind(l)
                         corr2int(ind(i))%elems(n)%v=ss
                      END IF
                   ENDDO
                ENDDO
             ENDDO
             corr2int(ind(i))%nelems=n
          END IF
       ENDDO

       CALL write_tensor4_to_file("corr2int.def",corr2int)
    CASE DEFAULT
       STOP '*** INT_CORR_MODE variable not properly defined in corrmod.nml ***'
    END SELECT

  END SUBROUTINE comp_corrint

END MODULE int_corr
