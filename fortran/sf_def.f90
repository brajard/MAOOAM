
! sf_def.f90
!
!>  Module to select the resolved-unresolved components.
!
!> @copyright                                                               
!> 2018 Jonathan Demaeyer
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------!

MODULE sf_def

  USE params, only: natm,noc,ndim
  USE util, only: str,rstr
  USE inprod_analytic, only:awavenum,owavenum
  IMPLICIT NONE

  PRIVATE

  LOGICAL :: exists !< Boolean to test for file existence.
  
  INTEGER, DIMENSION(:), ALLOCATABLE, PUBLIC :: SF              !< Unresolved variable definition vector
  INTEGER, DIMENSION(:), ALLOCATABLE, PUBLIC :: ind,rind        !< Unresolved reduction indices 
  INTEGER, DIMENSION(:), ALLOCATABLE, PUBLIC :: sl_ind,sl_rind  !< Resolved reduction indices 
  INTEGER, PUBLIC :: n_unres !< Number of unresolved variables
  INTEGER, PUBLIC :: n_res !< Number of resolved variables
  INTEGER, DIMENSION(:,:), ALLOCATABLE, PUBLIC :: Bar,Bau,Bor,Bou  !< Filter matrices 

  PUBLIC :: load_SF

CONTAINS

  !> Subroutine to load the unresolved variable defintion vector `SF` from SF.nml if it exists.
  !> If it does not, then write SF.nml with no unresolved variables specified (null vector).
  SUBROUTINE load_SF
    INTEGER :: i,AllocStat,n,ns
    CHARACTER(len=20) :: fm

    NAMELIST /SFlist/ SF

    fm(1:6)='(F3.1)'
   
    IF (ndim == 0) STOP "*** Number of dimensions is 0! ***"
    ALLOCATE(SF(0:ndim), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    INQUIRE(FILE='./SF.nml',EXIST=exists)

    IF (exists) THEN
       OPEN(8, file="SF.nml", status='OLD', recl=80, delim='APOSTROPHE')
       READ(8,nml=SFlist)
       CLOSE(8)
       n_unres=0
       DO i=1,ndim ! Computing the number of unresolved variables
          IF (SF(i)==1) n_unres=n_unres+1
       ENDDO
       IF (n_unres==0) STOP "*** No unresolved variable specified! ***"
       n_res=ndim-n_unres
       ALLOCATE(ind(n_unres), rind(0:ndim), sl_ind(n_res), sl_rind(0:ndim), STAT=AllocStat)
       IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
       ALLOCATE(Bar(0:ndim,0:ndim), Bau(0:ndim,0:ndim), Bor(0:ndim,0:ndim), Bou(0:ndim,0:ndim), STAT=AllocStat)
       IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
       rind=0
       n=1
       ns=1
       DO i=1,ndim
          IF (SF(i)==1) THEN
             ind(n)=i
             rind(i)=n
             n=n+1
          ELSE
             sl_ind(ns)=i
             sl_rind(i)=ns
             ns=ns+1
          ENDIF
       ENDDO
       Bar=0
       Bau=0
       Bor=0
       Bou=0
       DO i=1,2*natm
          IF (SF(i)==1) THEN
             Bau(i,i)=1
          ELSE
             Bar(i,i)=1
          ENDIF
       ENDDO
       DO i=2*natm+1,ndim
          IF (SF(i)==1) THEN
             Bou(i,i)=1
          ELSE
             Bor(i,i)=1
          ENDIF
       ENDDO
    ELSE
       OPEN(8, file="SF.nml", status='NEW')
       WRITE(8,'(a)') "!------------------------------------------------------------------------------!"
       WRITE(8,'(a)') "! Namelist file :                                                              !"
       WRITE(8,'(a)') "! Unresolved variables specification (1 -> unresolved, 0 -> resolved)                          !"
       WRITE(8,'(a)') "!------------------------------------------------------------------------------!"
       WRITE(8,*) ""
       WRITE(8,'(a)') "&SFLIST"
       WRITE(8,*) " ! psi variables"
       DO i=1,natm
          WRITE(8,*) " SF("//TRIM(str(i))//") = 0"//"   ! typ= "&
               &//awavenum(i)%typ//", Nx= "//TRIM(rstr(awavenum(i)&
               &%Nx,fm))//", Ny= "//TRIM(rstr(awavenum(i)%Ny,fm))
       END DO
       WRITE(8,*) " ! theta variables"
       DO i=1,natm
          WRITE(8,*) " SF("//TRIM(str(i+natm))//") = 0"//"   ! typ= "&
               &//awavenum(i)%typ//", Nx= "//TRIM(rstr(awavenum(i)&
               &%Nx,fm))//", Ny= "//TRIM(rstr(awavenum(i)%Ny,fm))
       END DO

       WRITE(8,*) " ! A variables"
       DO i=1,noc
          WRITE(8,*) " SF("//TRIM(str(i+2*natm))//") = 0"//"   ! Nx&
               &= "//TRIM(rstr(owavenum(i)%Nx,fm))//", Ny= "&
               &//TRIM(rstr(owavenum(i)%Ny,fm))
       END DO
       WRITE(8,*) " ! T variables"
       DO i=1,noc
          WRITE(8,*) " SF("//TRIM(str(i+noc+2*natm))//") = 0"//"   &
               &! Nx= "//TRIM(rstr(owavenum(i)%Nx,fm))//", Ny= "&
               &//TRIM(rstr(owavenum(i)%Ny,fm))
       END DO

       WRITE(8,'(a)') "&END"
       WRITE(8,*) ""
       CLOSE(8)
       STOP "*** SF.nml namelist written. Fill in the file and rerun !***"
    ENDIF
  END SUBROUTINE load_SF
END MODULE sf_def
