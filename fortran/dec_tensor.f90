
! dec_tensor.f90
!
!> The resolved-unresolved components decomposition of the tensor
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

MODULE dec_tensor

  !-----------------------------------------------------!
  !                                                     !
  ! Preamble and variables declaration                  !
  !                                                     !
  !-----------------------------------------------------!

  USE tensor, only:coolist
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: init_dec_tensor

  TYPE(coolist), DIMENSION(:), ALLOCATABLE, PUBLIC :: ff_tensor !< Tensor holding the part of the unresolved tensor involving only unresolved variables
  TYPE(coolist), DIMENSION(:), ALLOCATABLE, PUBLIC :: sf_tensor !< Tensor holding the part of the resolved tensor involving unresolved variables
  TYPE(coolist), DIMENSION(:), ALLOCATABLE, PUBLIC :: ss_tensor !< Tensor holding the part of the resolved tensor involving only resolved variables
  TYPE(coolist), DIMENSION(:), ALLOCATABLE, PUBLIC :: fs_tensor !< Tensor holding the part of the unresolved tensor involving resolved variables

  TYPE(coolist), DIMENSION(:), ALLOCATABLE, PUBLIC :: Hx !< Tensor holding the constant part of the resolved tendencies
  TYPE(coolist), DIMENSION(:), ALLOCATABLE, PUBLIC :: Lxx !< Tensor holding the linear part of the resolved tendencies involving the resolved variables
  TYPE(coolist), DIMENSION(:), ALLOCATABLE, PUBLIC :: Lxy !< Tensor holding the linear part of the resolved tendencies involving the unresolved variables
  TYPE(coolist), DIMENSION(:), ALLOCATABLE, PUBLIC :: Bxxx !< Tensor holding the quadratic part of the resolved tendencies involving resolved variables
  TYPE(coolist), DIMENSION(:), ALLOCATABLE, PUBLIC :: Bxxy !< Tensor holding the quadratic part of the resolved tendencies involving both resolved and unresolved variables
  TYPE(coolist), DIMENSION(:), ALLOCATABLE, PUBLIC :: Bxyy !< Tensor holding the quadratic part of the resolved tendencies involving unresolved variables

  TYPE(coolist), DIMENSION(:), ALLOCATABLE, PUBLIC :: Hy !< Tensor holding the constant part of the unresolved tendencies
  TYPE(coolist), DIMENSION(:), ALLOCATABLE, PUBLIC :: Lyx !< Tensor holding the linear part of the unresolved tendencies involving the resolved variables
  TYPE(coolist), DIMENSION(:), ALLOCATABLE, PUBLIC :: Lyy !< Tensor holding the linear part of the unresolved tendencies involving the unresolved variables
  TYPE(coolist), DIMENSION(:), ALLOCATABLE, PUBLIC :: Byxx !< Tensor holding the quadratic part of the unresolved tendencies involving resolved variables
  TYPE(coolist), DIMENSION(:), ALLOCATABLE, PUBLIC :: Byxy !< Tensor holding the quadratic part of the unresolved tendencies involving both resolved and unresolved variables
  TYPE(coolist), DIMENSION(:), ALLOCATABLE, PUBLIC :: Byyy !< Tensor holding the quadratic part of the unresolved tendencies involving unresolved variables

  TYPE(coolist), DIMENSION(:), ALLOCATABLE, PUBLIC :: ss_tl_tensor !< Tensor of the tangent linear model tendencies of the resolved component alone
 
 
  TYPE(coolist), DIMENSION(:), ALLOCATABLE :: dumb !< Dumb coolist to make the computations


  !-----------------------------------------------------!
  !                                                     !
  ! End of preamble                                     !
  !                                                     !
  !-----------------------------------------------------!

CONTAINS

  !-----------------------------------------------------!
  !                                                     !
  ! Function declarations                               !
  !                                                     !
  !-----------------------------------------------------!

  !> Subroutine to suppress from the tensor \f$t_{ijk}\f$ components satisfying SF(j)=v1
  !> and SF(k)=v2
  !> @param t tensor over which the routine acts
  !> @param cst constant which controls if the 0 index is taken as a unresolved or a resolved one
  !> @param v1 first constant of the conditional (0 to suppress resolved, 1 for unresolved)
  !> @param v2 second constant of the conditional (0 to suppress resolved, 1 for unresolved)
  SUBROUTINE suppress_and(t,cst,v1,v2)
    USE params, only:ndim
    USE sf_def, only: SF
    TYPE(coolist), DIMENSION(ndim), INTENT(INOUT) :: t
    INTEGER, INTENT(IN) :: cst,v1,v2
    INTEGER :: i,n,li,liii
    
    SF(0)=cst ! control wether 0 index is considered unresolved or not
    DO i=1,ndim
       n=t(i)%nelems
       DO li=1,n
          ! Clear entries with only resolved variables and shift rest of the items one place down.
          ! Make sure not to skip any entries while shifting!
             
          DO WHILE ((SF(t(i)%elems(li)%j)==v1).and.(SF(t(i)%elems(li)%k)==v2))
             !print*, i,li,t(i)%nelems,n
             DO liii=li+1,n
                t(i)%elems(liii-1)%j=t(i)%elems(liii)%j
                t(i)%elems(liii-1)%k=t(i)%elems(liii)%k
                t(i)%elems(liii-1)%v=t(i)%elems(liii)%v
             ENDDO
             t(i)%nelems=t(i)%nelems-1
             IF (li>t(i)%nelems) exit
          ENDDO
          IF (li>t(i)%nelems) exit
       ENDDO
    ENDDO

  END SUBROUTINE suppress_and

  !> Subroutine to suppress from the tensor \f$t_{ijk}\f$ components satisfying SF(j)=v1
  !> or SF(k)=v2
  !> @param t tensor over which the routine acts
  !> @param cst constant which controls if the 0 index is taken as a unresolved or a resolved one
  !> @param v1 first constant of the conditional (0 to suppress resolved, 1 for unresolved)
  !> @param v2 second constant of the conditional (0 to suppress resolved, 1 for unresolved)
  SUBROUTINE suppress_or(t,cst,v1,v2)
    USE params, only:ndim
    USE sf_def, only: SF
    TYPE(coolist), DIMENSION(ndim), INTENT(INOUT) :: t
    INTEGER, INTENT(IN) :: cst,v1,v2
    INTEGER :: i,n,li,liii
    
    SF(0)=cst ! control wether 0 index is considered unresolved or not
    DO i=1,ndim
       n=t(i)%nelems
       DO li=1,n
          ! Clear entries with only resolved variables and shift rest of the items one place down.
          ! Make sure not to skip any entries while shifting!
             
          DO WHILE ((SF(t(i)%elems(li)%j)==v1).or.(SF(t(i)%elems(li)%k)==v2))
             !print*, i,li,t(i)%nelems,n
             DO liii=li+1,n
                t(i)%elems(liii-1)%j=t(i)%elems(liii)%j
                t(i)%elems(liii-1)%k=t(i)%elems(liii)%k
                t(i)%elems(liii-1)%v=t(i)%elems(liii)%v
             ENDDO
             t(i)%nelems=t(i)%nelems-1
             IF (li>t(i)%nelems) exit
          ENDDO
          IF (li>t(i)%nelems) exit
       ENDDO
    ENDDO

  END SUBROUTINE suppress_or

  !> Subroutine to reorder the tensor \f$t_{ijk}\f$ components : if SF(j)=v then it
  !> return \f$t_{ikj}\f$
  !> @param t tensor over which the routine acts
  !> @param cst constant which controls if the 0 index is taken as a unresolved or a resolved one
  !> @param v constant of the conditional (0 to invert resolved, 1 for unresolved)
  SUBROUTINE reorder(t,cst,v)
    USE params, only:ndim
    USE sf_def, only: SF
    TYPE(coolist), DIMENSION(ndim), INTENT(INOUT) :: t
    INTEGER, INTENT(IN) :: cst,v
    INTEGER :: i,n,li,liii

    SF(0)=cst ! control wether 0 index is considered unresolved or not
    DO i=1,ndim

       n=t(i)%nelems
       DO li=1,n
          IF (SF(t(i)%elems(li)%j)==v) THEN
             liii=t(i)%elems(li)%j
             t(i)%elems(li)%j=t(i)%elems(li)%k
             t(i)%elems(li)%k=liii
          ENDIF
       ENDDO
    ENDDO

  END SUBROUTINE reorder

  !> Subroutine that suppress all the components of a tensor \f$t_{ijk}\f$ where if SF(i)=v 
  !> @param t tensor over which the routine acts
  !> @param cst constant which controls if the 0 index is taken as a unresolved or a resolved one
  !> @param v constant of the conditional (0 to suppress resolved, 1 for unresolved)
  SUBROUTINE init_sub_tensor(t,cst,v)
    USE params, only:ndim
    USE sf_def, only: SF
    TYPE(coolist), DIMENSION(ndim), INTENT(INOUT) :: t
    INTEGER, INTENT(IN) :: cst,v
    INTEGER :: i

    SF(0)=cst ! control wether 0 index is considered unresolved or not
    DO i=1,ndim
       IF (SF(i)==v) t(i)%nelems=0
    ENDDO

  END SUBROUTINE init_sub_tensor

  !-----------------------------------------------------!
  !                                                     !
  ! Initialisation routine                              !
  !                                                     !
  !-----------------------------------------------------!

  !> Subroutine that initialize and compute the decomposed tensors
  SUBROUTINE init_dec_tensor
    USE params, only:ndim
    USE aotensor_def, only:aotensor
    USE sf_def, only: load_SF
    USE tensor, only:copy_coo,add_to_tensor,scal_mul_coo
    USE tl_ad_tensor, only: init_tltensor,tltensor
    USE stoch_params, only: init_stoch_params,mode,tdelta,eps_pert
    INTEGER :: AllocStat 

    CALL init_stoch_params

    CALL init_tltensor  ! and tl tensor
    
    CALL load_SF  ! Load the resolved-unresolved decomposition

    ! Allocating the returned arrays

    ALLOCATE(ff_tensor(ndim),fs_tensor(ndim),sf_tensor(ndim),ss_tensor(ndim), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    ALLOCATE(ss_tl_tensor(ndim), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    ALLOCATE(Hx(ndim),Lxx(ndim),Lxy(ndim),Bxxx(ndim),Bxxy(ndim),Bxyy(ndim), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    ALLOCATE(Hy(ndim),Lyx(ndim),Lyy(ndim),Byxx(ndim),Byxy(ndim),Byyy(ndim), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    ! General decomposition
    ! ff tensor
    ALLOCATE(dumb(ndim), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    IF (mode.ne.'qfst') THEN
        CALL copy_coo(aotensor,dumb) !Copy the tensors
        CALL init_sub_tensor(dumb,0,0)
        CALL suppress_or(dumb,1,0,0) ! Clear entries with resolved variables
        CALL copy_coo(dumb,ff_tensor)
    ELSE
        CALL copy_coo(aotensor,dumb) !Copy the tensors
        CALL init_sub_tensor(dumb,0,0)
        CALL suppress_or(dumb,0,0,0) ! Clear entries with resolved variables and linear and constant terms
        CALL copy_coo(dumb,ff_tensor)
    ENDIF

    AllocStat=0
    DEALLOCATE(dumb, STAT=AllocStat)
    IF (AllocStat /= 0)  STOP "*** Problem to deallocate ! ***"

    ! fs tensor
    ALLOCATE(dumb(ndim), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    IF (mode.ne.'qfst') THEN
        CALL copy_coo(aotensor,dumb) !Copy the tensors
        CALL init_sub_tensor(dumb,0,0)
        CALL suppress_and(dumb,1,1,1) ! Clear entries with only unresolved variables and constant
        CALL copy_coo(dumb,fs_tensor)
    ELSE
        CALL copy_coo(aotensor,dumb) !Copy the tensors
        CALL init_sub_tensor(dumb,0,0)
        CALL suppress_and(dumb,0,1,1) ! Clear entries with only quadratic unresolved variables 
        CALL copy_coo(dumb,fs_tensor)
    ENDIF

    AllocStat=0
    DEALLOCATE(dumb, STAT=AllocStat)
    IF (AllocStat /= 0)  STOP "*** Problem to deallocate ! ***"

    ! sf tensor
    ALLOCATE(dumb(ndim), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"


    CALL copy_coo(aotensor,dumb) !Copy the tensors
    CALL init_sub_tensor(dumb,1,1)
    CALL suppress_and(dumb,0,0,0) ! Clear entries with only unresolved variables and constant
    CALL copy_coo(dumb,sf_tensor)

    AllocStat=0
    DEALLOCATE(dumb, STAT=AllocStat)
    IF (AllocStat /= 0)  STOP "*** Problem to deallocate ! ***"

    ! ss tensor
    ALLOCATE(dumb(ndim), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"


    CALL copy_coo(aotensor,dumb) !Copy the tensors
    CALL init_sub_tensor(dumb,1,1)
    CALL suppress_or(dumb,0,1,1) ! Clear entries with only unresolved variables and constant
    CALL copy_coo(dumb,ss_tensor)

    AllocStat=0
    DEALLOCATE(dumb, STAT=AllocStat)
    IF (AllocStat /= 0)  STOP "*** Problem to deallocate ! ***"

    ! ss tangent linear tensor

    ALLOCATE(dumb(ndim), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    CALL copy_coo(tltensor,dumb) !Copy the tensors
    CALL init_sub_tensor(dumb,1,1)
    CALL suppress_or(dumb,0,1,1) ! Clear entries with only unresolved variables and constant
    CALL copy_coo(dumb,ss_tl_tensor)

    AllocStat=0
    DEALLOCATE(dumb, STAT=AllocStat)
    IF (AllocStat /= 0)  STOP "*** Problem to deallocate ! ***"

    ! Multiply the aotensor part that need to be by the perturbation and time
    ! separation parameter
    
    ALLOCATE(dumb(ndim), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    CALL copy_coo(ss_tensor,dumb)
    CALL scal_mul_coo(1.D0/tdelta**2,ff_tensor)
    CALL scal_mul_coo(eps_pert/tdelta,fs_tensor)
    CALL add_to_tensor(ff_tensor,dumb)
    CALL add_to_tensor(fs_tensor,dumb)
    CALL scal_mul_coo(eps_pert/tdelta,sf_tensor)
    CALL add_to_tensor(sf_tensor,dumb)

    AllocStat=0
    DEALLOCATE(aotensor, STAT=AllocStat)
    IF (AllocStat /= 0)  STOP "*** Problem to deallocate ! ***"

    ALLOCATE(aotensor(ndim), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    CALL copy_coo(dumb,aotensor)

    AllocStat=0
    DEALLOCATE(dumb, STAT=AllocStat)
    IF (AllocStat /= 0)  STOP "*** Problem to deallocate ! ***"

    ! MTV decomposition
    ! Unresolved tensors

    ! Hy tensor
    ALLOCATE(dumb(ndim), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    CALL copy_coo(aotensor,dumb) !Copy the tensors
    CALL init_sub_tensor(dumb,0,0)
    CALL suppress_or(dumb,0,1,1) ! Clear entries with unresolved variables
    CALL suppress_or(dumb,1,0,0) ! Suppress linear and nonlinear resolved terms
    CALL copy_coo(dumb,Hy)

    AllocStat=0
    DEALLOCATE(dumb, STAT=AllocStat)
    IF (AllocStat /= 0)  STOP "*** Problem to deallocate ! ***"

    ! Lyx tensor
    ALLOCATE(dumb(ndim), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    CALL copy_coo(aotensor,dumb) !Copy the tensors
    CALL init_sub_tensor(dumb,0,0)
    CALL suppress_or(dumb,0,1,1) ! Clear entries with unresolved variables
    CALL suppress_and(dumb,1,1,1) ! Clear constant entries
    CALL suppress_and(dumb,1,0,0) ! Clear entries with nonlinear resolved terms
    CALL reorder(dumb,1,0) ! Resolved variables must be the third (k) index
    CALL copy_coo(dumb,Lyx)

    AllocStat=0
    DEALLOCATE(dumb, STAT=AllocStat)
    IF (AllocStat /= 0)  STOP "*** Problem to deallocate ! ***"

    ! Lyy tensor
    ALLOCATE(dumb(ndim), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    CALL copy_coo(aotensor,dumb) !Copy the tensors
    CALL init_sub_tensor(dumb,0,0)
    CALL suppress_or(dumb,1,0,0) ! Clear entries with resolved variables
    CALL suppress_and(dumb,0,1,1) ! Clear entries with nonlinear unresolved terms
    CALL suppress_and(dumb,0,0,0) ! Clear constant entries
    CALL reorder(dumb,0,1) ! Unresolved variables must be the third (k) index
    CALL copy_coo(dumb,Lyy)

    AllocStat=0
    DEALLOCATE(dumb, STAT=AllocStat)
    IF (AllocStat /= 0)  STOP "*** Problem to deallocate ! ***"

    ! Byxy tensor
    ALLOCATE(dumb(ndim), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    CALL copy_coo(aotensor,dumb) !Copy the tensors
    CALL init_sub_tensor(dumb,0,0)
    CALL suppress_and(dumb,1,1,1) ! Clear constant or linear terms and nonlinear unresolved only entries
    CALL suppress_and(dumb,0,0,0) ! Clear entries with only resolved variables and constant
    CALL reorder(dumb,0,1) ! Unresolved variables must be the third (k) index
    CALL copy_coo(dumb,Byxy)

    AllocStat=0
    DEALLOCATE(dumb, STAT=AllocStat)
    IF (AllocStat /= 0)  STOP "*** Problem to deallocate ! ***"

    ! Byyy tensor
    ALLOCATE(dumb(ndim), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    CALL copy_coo(aotensor,dumb) !Copy the tensors
    CALL init_sub_tensor(dumb,0,0)
    CALL suppress_or(dumb,0,0,0) ! Clear entries with resolved variables and linear and constant terms
    CALL copy_coo(dumb,Byyy)

    AllocStat=0
    DEALLOCATE(dumb, STAT=AllocStat)
    IF (AllocStat /= 0)  STOP "*** Problem to deallocate ! ***"

    ! Byxx tensor
    ALLOCATE(dumb(ndim), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    CALL copy_coo(aotensor,dumb) !Copy the tensors
    CALL init_sub_tensor(dumb,0,0)
    CALL suppress_or(dumb,1,1,1) ! Clear entries with unresolved variables and linear and constant terms
    CALL copy_coo(dumb,Byxx)

    AllocStat=0
    DEALLOCATE(dumb, STAT=AllocStat)
    IF (AllocStat /= 0)  STOP "*** Problem to deallocate ! ***"

    ! Resolved tensors

    ! Hx tensor
    ALLOCATE(dumb(ndim), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    CALL copy_coo(aotensor,dumb) !Copy the tensors
    CALL init_sub_tensor(dumb,1,1)
    CALL suppress_or(dumb,1,0,0) ! Clear entries with resolved variables
    CALL suppress_or(dumb,0,1,1) ! Suppress linear and nonlinear unresolved terms
    CALL copy_coo(dumb,Hx)

    AllocStat=0
    DEALLOCATE(dumb, STAT=AllocStat)
    IF (AllocStat /= 0)  STOP "*** Problem to deallocate ! ***"

    ! Lxy tensor
    ALLOCATE(dumb(ndim), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    CALL copy_coo(aotensor,dumb) !Copy the tensors
    CALL init_sub_tensor(dumb,1,1)
    CALL suppress_or(dumb,1,0,0) ! Clear entries with resolved variables
    CALL suppress_and(dumb,0,0,0) ! Clear constant entries
    CALL suppress_and(dumb,0,1,1) ! Clear entries with nonlinear unresolved terms
    CALL reorder(dumb,0,1) ! Resolved variables must be the third (k) index
    CALL copy_coo(dumb,Lxy)

    AllocStat=0
    DEALLOCATE(dumb, STAT=AllocStat)
    IF (AllocStat /= 0)  STOP "*** Problem to deallocate ! ***"

    ! Lxx tensor
    ALLOCATE(dumb(ndim), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    CALL copy_coo(aotensor,dumb) !Copy the tensors
    CALL init_sub_tensor(dumb,1,1)
    CALL suppress_or(dumb,0,1,1) ! Clear entries with unresolved variables
    CALL suppress_and(dumb,1,0,0) ! Clear entries with nonlinear resolved terms
    CALL suppress_and(dumb,1,1,1) ! Clear constant entries
    CALL reorder(dumb,1,0) ! Resolved variables must be the third (k) index
    CALL copy_coo(dumb,Lxx)

    AllocStat=0
    DEALLOCATE(dumb, STAT=AllocStat)
    IF (AllocStat /= 0)  STOP "*** Problem to deallocate ! ***"

    ! Bxxy tensor
    ALLOCATE(dumb(ndim), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    CALL copy_coo(aotensor,dumb) !Copy the tensors
    CALL init_sub_tensor(dumb,1,1)
    CALL suppress_and(dumb,1,1,1) ! Clear constant or linear terms and nonlinear unresolved only entries
    CALL suppress_and(dumb,0,0,0) ! Clear entries with only resolved variables and constant
    CALL reorder(dumb,0,1) ! Unresolved variables must be the third (k) index
    CALL copy_coo(dumb,Bxxy)

    AllocStat=0
    DEALLOCATE(dumb, STAT=AllocStat)
    IF (AllocStat /= 0)  STOP "*** Problem to deallocate ! ***"

    ! Bxxx tensor
    ALLOCATE(dumb(ndim), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    CALL copy_coo(aotensor,dumb) !Copy the tensors
    CALL init_sub_tensor(dumb,1,1)
    CALL suppress_or(dumb,1,1,1) ! Clear entries with unresolved variables and linear and constant terms
    CALL copy_coo(dumb,Bxxx)

    AllocStat=0
    DEALLOCATE(dumb, STAT=AllocStat)
    IF (AllocStat /= 0)  STOP "*** Problem to deallocate ! ***"

    ! Bxyy tensor
    ALLOCATE(dumb(ndim), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    CALL copy_coo(aotensor,dumb) !Copy the tensors
    CALL init_sub_tensor(dumb,1,1)
    CALL suppress_or(dumb,0,0,0) ! Clear entries with resolved variables and linear and constant terms
    CALL copy_coo(dumb,Bxyy)

    AllocStat=0
    DEALLOCATE(dumb, STAT=AllocStat)
    IF (AllocStat /= 0)  STOP "*** Problem to deallocate ! ***"



  END SUBROUTINE init_dec_tensor
END MODULE dec_tensor



