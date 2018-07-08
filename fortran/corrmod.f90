
! corrmod.f90
!
!> Module to initialize the correlation matrix of the unresolved variables
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

MODULE corrmod
  USE params, only: ndim
  USE stoch_params, only:  q_ar,q_au,q_or,q_ou,maxint,load_mode,int_corr_mode
  USE sf_def, only: n_unres,ind
  USE util, only: invmat

  PRIVATE

  PUBLIC :: init_corr

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, PUBLIC :: mean !< Vector holding the mean of the unresolved dynamics (reduced version)
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, PUBLIC :: mean_full !< Vector holding the mean of the unresolved dynamics (full version)
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: corr_i_full !< Covariance matrix of the unresolved variables (full version)
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: inv_corr_i_full !< Inverse of the covariance matrix of the unresolved variables (full version)
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: corr_i !< Covariance matrix of the unresolved variables (reduced version)
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: inv_corr_i !< Inverse of the covariance matrix of the unresolved variables (reduced version)
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: corr_ij !< Matrix holding the correlation matrix at a given time
  REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: y2 !< Vector holding coefficient of the spline and exponential correlation representation
  REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: ya !< Vector holding coefficient of the spline and exponential correlation representation
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: xa !< Vector holding coefficient of the spline and exponential correlation representation
 
  !> Integers needed by the spline representation of the correlation 
  INTEGER :: nspl,klo,khi

  !> Pointer to the correlation computation routine
  PROCEDURE(corrcomp_from_spline), POINTER, PUBLIC :: corrcomp

CONTAINS
  !> Subroutine to initialise the computation of the correlation
  SUBROUTINE init_corr
    INTEGER :: AllocStat,i,j,k,nf
    REAL(KIND=8), DIMENSION(5) :: dumb
    LOGICAL :: ex

    ! Selection of the loading mode
    SELECT CASE (load_mode)
    CASE ('defi')
       corrcomp => corrcomp_from_def
    CASE ('spli')
       INQUIRE(FILE='corrspline.def',EXIST=ex)
       IF (.not.ex) STOP "*** File corrspline.def not found ! ***"
       OPEN(20,file='corrspline.def',status='old')
       READ(20,*) nf,nspl
       IF (nf /= n_unres) STOP "*** Dimension in files corrspline.def and sf.nml do not correspond ! ***"
       ALLOCATE(xa(nspl), ya(n_unres,n_unres,nspl), y2(n_unres,n_unres,nspl), STAT=AllocStat)
       IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
       READ(20,*) xa
       maxint=xa(nspl)/2
       DO k=1,n_unres*n_unres
          READ(20,*) i,j
          READ(20,*) ya(i,j,:)
          READ(20,*) y2(i,j,:)
       ENDDO
       CLOSE(20)
       corrcomp => corrcomp_from_spline
       klo=1
       khi=nspl
    CASE ('expo')
       INQUIRE(FILE='correxpo.def',EXIST=ex)
       IF (.not.ex) STOP "*** File correxpo.def not found ! ***"
       OPEN(20,file='correxpo.def',status='old')
       READ(20,*) nf,maxint
       IF (nf /= n_unres) STOP "*** Dimension in files correxpo.def and sf.nml do not correspond ! ***"
       ALLOCATE(ya(n_unres,n_unres,5), STAT=AllocStat)
       IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
       DO k=1,n_unres*n_unres
          READ(20,*) i,j,dumb
          ya(i,j,:)=dumb
       ENDDO
       CLOSE(20)
       corrcomp => corrcomp_from_fit
    CASE DEFAULT
       STOP '*** LOAD_MODE variable not properly defined in corrmod.nml ***'
    END SELECT

    ALLOCATE(mean(n_unres),mean_full(0:ndim), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
    
    ALLOCATE(inv_corr_i(n_unres,n_unres), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    ALLOCATE(corr_i(n_unres,n_unres), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    ALLOCATE(corr_ij(n_unres,n_unres), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    ALLOCATE(corr_i_full(ndim,ndim), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    ALLOCATE(inv_corr_i_full(ndim,ndim), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    corr_ij=0.D0

    CALL corrcomp(0.D0)
    corr_i=corr_ij
    inv_corr_i=invmat(corr_i)

    corr_i_full=0.D0
    DO i=1,n_unres
       DO j=1,n_unres
          corr_i_full(ind(i),ind(j))=corr_i(i,j)
       ENDDO
    ENDDO
    
    inv_corr_i_full=0.D0
    DO i=1,n_unres
       DO j=1,n_unres
          inv_corr_i_full(ind(i),ind(j))=inv_corr_i(i,j)
       ENDDO
    ENDDO

    mean=0.D0
    INQUIRE(FILE='mean.def',EXIST=ex)
    IF (ex) THEN
       OPEN(20,file='mean.def',status='old')
       READ(20,*) mean
       CLOSE(20)
    ENDIF

    mean_full=0.D0
    DO i=1,n_unres
       mean_full(ind(i))=mean(i)
    ENDDO

  END SUBROUTINE init_corr
  
  !> Subroutine to compute the correlation of the unresolved variables \f$\langle Y \otimes Y^s \rangle\f$
  !> at time \f$s\f$ from the definition given inside the module
  !> @param s time \f$s\f$ at which the correlation is computed
  SUBROUTINE corrcomp_from_def(s)
    REAL(KIND=8), INTENT(IN) :: s
    REAL(KIND=8) :: y
    INTEGER :: i,j

    ! Definition of the corr_ij matrix as a function of time

    corr_ij(10,10)=((7.66977252307669*Cos(0.039597160512071454*s))/exp(0.06567483898489704*s) + (1.0240906173830213*Cos(&
         &0.07283568782600224*s))/exp(0.017262015588746404*s) - (0.6434985372062336*Sin(0.039597160512071454*s&
         &))/exp(0.06567483898489704*s) + (0.6434985372062335*Sin(0.07283568782600224*s))/exp(0.017262015588746404*s)) 
    corr_ij(10,9)=((0.6434985372062321*Cos(0.039597160512071454*s))/exp(0.06567483898489704*s) - (0.6434985372062324*Co&
         &s(0.07283568782600224*s))/exp(0.017262015588746404*s) + (7.669772523076694*Sin(0.039597160512071454*&
         &s))/exp(0.06567483898489704*s) + (1.024090617383021*Sin(0.07283568782600224*s))/exp(0.017262015588746404*s)) 
    corr_ij(10,8)=0
    corr_ij(10,7)=0
    corr_ij(10,6)=0
    corr_ij(10,5)=((-2.236364132659011*Cos(0.039597160512071454*s))/exp(0.06567483898489704*s) + (6.952804148086198*Cos&
         &(0.07283568782600224*s))/exp(0.017262015588746404*s) - (1.4494534432272481*Sin(0.039597160512071454*&
         &s))/exp(0.06567483898489704*s) - (0.6818177416446283*Sin(0.07283568782600224*s))/exp(0.017262015588746404*s)) 
    corr_ij(10,4)=((1.4494534432272483*Cos(0.039597160512071454*s))/exp(0.06567483898489704*s) + (0.6818177416446293*Co&
         &s(0.07283568782600224*s))/exp(0.017262015588746404*s) - (2.2363641326590127*Sin(0.039597160512071454&
         &*s))/exp(0.06567483898489704*s) + (6.952804148086195*Sin(0.07283568782600224*s))/exp(0.017262015588746404*s)) 
    corr_ij(10,3)=0
    corr_ij(10,2)=0
    corr_ij(10,1)=0
    corr_ij(9,10)=((-0.6434985372062344*Cos(0.039597160512071454*s))/exp(0.06567483898489704*s) + (0.643498537206234*Co&
         &s(0.07283568782600224*s))/exp(0.017262015588746404*s) - (7.669772523076689*Sin(0.039597160512071454*&
         &s))/exp(0.06567483898489704*s) - (1.0240906173830204*Sin(0.07283568782600224*s))/exp(0.017262015588746404*s)) 
    corr_ij(9,9)=((7.66977252307669*Cos(0.039597160512071454*s))/exp(0.06567483898489704*s) + (1.0240906173830204*Cos(&
         &0.07283568782600224*s))/exp(0.017262015588746404*s) - (0.643498537206233*Sin(0.039597160512071454*s)&
         &)/exp(0.06567483898489704*s) + (0.6434985372062327*Sin(0.07283568782600224*s))/exp(0.017262015588746404*s)) 
    corr_ij(9,8)=0
    corr_ij(9,7)=0
    corr_ij(9,6)=0
    corr_ij(9,5)=((-1.4494534432272477*Cos(0.039597160512071454*s))/exp(0.06567483898489704*s) - (0.6818177416446249*C&
         &os(0.07283568782600224*s))/exp(0.017262015588746404*s) + (2.2363641326590105*Sin(0.03959716051207145&
         &4*s))/exp(0.06567483898489704*s) - (6.952804148086195*Sin(0.07283568782600224*s))/exp(0.017262015588746404*s)) 
    corr_ij(9,4)=((-2.2363641326590127*Cos(0.039597160512071454*s))/exp(0.06567483898489704*s) + (6.952804148086194*Co&
         &s(0.07283568782600224*s))/exp(0.017262015588746404*s) - (1.4494534432272486*Sin(0.039597160512071454&
         &*s))/exp(0.06567483898489704*s) - (0.6818177416446249*Sin(0.07283568782600224*s))/exp(0.017262015588746404*s)) 
    corr_ij(9,3)=0
    corr_ij(9,2)=0
    corr_ij(9,1)=0
    corr_ij(8,10)=0

    corr_ij(8,9)=0
    corr_ij(8,8)=(9.135647293470983/exp(0.05076718239027029*s) + 2.2233889637758932/exp(0.016285467000648854*s)) 
    corr_ij(8,7)=0
    corr_ij(8,6)=0
    corr_ij(8,5)=0
    corr_ij(8,4)=0
    corr_ij(8,3)=(-5.938566084855411/exp(0.05076718239027029*s) + 11.97129741027622/exp(0.016285467000648854*s)) 
    corr_ij(8,2)=0
    corr_ij(8,1)=0
    corr_ij(7,10)=0

    corr_ij(7,9)=0
    corr_ij(7,8)=0
    corr_ij(7,7)=((11.518026982819887*Cos(0.029341409726871926*s))/exp(0.04435489221745234*s) + (0.05351107793747755*C&
         &os(0.11425932545092894*s))/exp(0.019700737327669783*s) - (0.14054811601869432*Sin(0.0293414097268719&
         &26*s))/exp(0.04435489221745234*s) + (0.14054811601869702*Sin(0.11425932545092894*s))/exp(0.019700737327669783*s)) 
    corr_ij(7,6)=((0.14054811601869532*Cos(0.029341409726871926*s))/exp(0.04435489221745234*s) - (0.14054811601869702*&
         &Cos(0.11425932545092894*s))/exp(0.019700737327669783*s) + (11.518026982819887*Sin(0.0293414097268719&
         &26*s))/exp(0.04435489221745234*s) + (0.0535110779374777*Sin(0.11425932545092894*s))/exp(0.019700737327669783*s)) 
    corr_ij(7,5)=0
    corr_ij(7,4)=0
    corr_ij(7,3)=0
    corr_ij(7,2)=((-0.732907009016115*Cos(0.029341409726871926*s))/exp(0.04435489221745234*s) + (2.728845031386875*Cos&
         &(0.11425932545092894*s))/exp(0.019700737327669783*s) - (2.4717920234033532*Sin(0.029341409726871926*&
         &s))/exp(0.04435489221745234*s) - (0.24003801347124257*Sin(0.11425932545092894*s))/exp(0.019700737327669783*s)) 
    corr_ij(7,1)=((2.4717920234033532*Cos(0.029341409726871926*s))/exp(0.04435489221745234*s) + (0.2400380134712426*Co&
         &s(0.11425932545092894*s))/exp(0.019700737327669783*s) - (0.7329070090161153*Sin(0.029341409726871926&
         &*s))/exp(0.04435489221745234*s) + (2.728845031386876*Sin(0.11425932545092894*s))/exp(0.019700737327669783*s)) 
    corr_ij(6,10)=0

    corr_ij(6,9)=0
    corr_ij(6,8)=0
    corr_ij(6,7)=((-0.1405481160186977*Cos(0.029341409726871926*s))/exp(0.04435489221745234*s) + (0.14054811601869713*&
         &Cos(0.11425932545092894*s))/exp(0.019700737327669783*s) - (11.518026982819885*Sin(0.0293414097268719&
         &26*s))/exp(0.04435489221745234*s) - (0.05351107793747755*Sin(0.11425932545092894*s))/exp(0.019700737327669783*s)) 
    corr_ij(6,6)=((11.518026982819885*Cos(0.029341409726871926*s))/exp(0.04435489221745234*s) + (0.05351107793747768*C&
         &os(0.11425932545092894*s))/exp(0.019700737327669783*s) - (0.14054811601869832*Sin(0.0293414097268719&
         &26*s))/exp(0.04435489221745234*s) + (0.14054811601869707*Sin(0.11425932545092894*s))/exp(0.019700737327669783*s)) 
    corr_ij(6,5)=0
    corr_ij(6,4)=0
    corr_ij(6,3)=0
    corr_ij(6,2)=((-2.471792023403353*Cos(0.029341409726871926*s))/exp(0.04435489221745234*s) - (0.2400380134712425*Co&
         &s(0.11425932545092894*s))/exp(0.019700737327669783*s) + (0.7329070090161155*Sin(0.029341409726871926&
         &*s))/exp(0.04435489221745234*s) - (2.7288450313868755*Sin(0.11425932545092894*s))/exp(0.019700737327669783*s)) 
    corr_ij(6,1)=((-0.7329070090161154*Cos(0.029341409726871926*s))/exp(0.04435489221745234*s) + (2.728845031386876*Co&
         &s(0.11425932545092894*s))/exp(0.019700737327669783*s) - (2.4717920234033524*Sin(0.029341409726871926&
         &*s))/exp(0.04435489221745234*s) - (0.24003801347124343*Sin(0.11425932545092894*s))/exp(0.019700737327669783*s)) 
    corr_ij(5,10)=((0.5794534449999711*Cos(0.039597160512071454*s))/exp(0.06567483898489704*s) + (4.136986570427212*Cos&
         &(0.07283568782600224*s))/exp(0.017262015588746404*s) - (1.0360597341248128*Sin(0.039597160512071454*&
         &s))/exp(0.06567483898489704*s) + (3.167330918996692*Sin(0.07283568782600224*s))/exp(0.017262015588746404*s)) 
    corr_ij(5,9)=((1.0360597341248134*Cos(0.039597160512071454*s))/exp(0.06567483898489704*s) - (3.1673309189966856*Co&
         &s(0.07283568782600224*s))/exp(0.017262015588746404*s) + (0.5794534449999746*Sin(0.039597160512071454&
         &*s))/exp(0.06567483898489704*s) + (4.1369865704272115*Sin(0.07283568782600224*s))/exp(0.017262015588746404*s)) 
    corr_ij(5,8)=0
    corr_ij(5,7)=0
    corr_ij(5,6)=0
    corr_ij(5,5)=((-0.37825091063447547*Cos(0.039597160512071454*s))/exp(0.06567483898489704*s) + (30.094690926061638*&
         &Cos(0.07283568782600224*s))/exp(0.017262015588746404*s) + (0.16085380971100194*Sin(0.039597160512071&
         &454*s))/exp(0.06567483898489704*s) - (0.1608538097109995*Sin(0.07283568782600224*s))/exp(0.017262015588746404*s)) 
    corr_ij(5,4)=((-0.16085380971100238*Cos(0.039597160512071454*s))/exp(0.06567483898489704*s) + (0.16085380971100127&
         &*Cos(0.07283568782600224*s))/exp(0.017262015588746404*s) - (0.37825091063447586*Sin(0.03959716051207&
         &1454*s))/exp(0.06567483898489704*s) + (30.09469092606163*Sin(0.07283568782600224*s))/exp(0.017262015588746404*s)) 
    corr_ij(5,3)=0
    corr_ij(5,2)=0
    corr_ij(5,1)=0
    corr_ij(4,10)=((-1.0360597341248106*Cos(0.039597160512071454*s))/exp(0.06567483898489704*s) + (3.167330918996689*Co&
         &s(0.07283568782600224*s))/exp(0.017262015588746404*s) - (0.5794534449999716*Sin(0.039597160512071454&
         &*s))/exp(0.06567483898489704*s) - (4.1369865704272115*Sin(0.07283568782600224*s))/exp(0.017262015588746404*s)) 
    corr_ij(4,9)=((0.5794534449999711*Cos(0.039597160512071454*s))/exp(0.06567483898489704*s) + (4.1369865704272115*Co&
         &s(0.07283568782600224*s))/exp(0.017262015588746404*s) - (1.0360597341248114*Sin(0.039597160512071454&
         &*s))/exp(0.06567483898489704*s) + (3.1673309189966843*Sin(0.07283568782600224*s))/exp(0.017262015588746404*s)) 
    corr_ij(4,8)=0
    corr_ij(4,7)=0
    corr_ij(4,6)=0
    corr_ij(4,5)=((0.16085380971100194*Cos(0.039597160512071454*s))/exp(0.06567483898489704*s) - (0.16085380971100371*&
         &Cos(0.07283568782600224*s))/exp(0.017262015588746404*s) + (0.37825091063447497*Sin(0.039597160512071&
         &454*s))/exp(0.06567483898489704*s) - (30.094690926061617*Sin(0.07283568782600224*s))/exp(0.017262015588746404*s)) 
    corr_ij(4,4)=((-0.37825091063447536*Cos(0.039597160512071454*s))/exp(0.06567483898489704*s) + (30.094690926061617*&
         &Cos(0.07283568782600224*s))/exp(0.017262015588746404*s) + (0.16085380971100172*Sin(0.039597160512071&
         &454*s))/exp(0.06567483898489704*s) - (0.16085380971100616*Sin(0.07283568782600224*s))/exp(0.017262015588746404*s)) 
    corr_ij(4,3)=0
    corr_ij(4,2)=0
    corr_ij(4,1)=0
    corr_ij(3,10)=0

    corr_ij(3,9)=0
    corr_ij(3,8)=(0.24013456462471527/exp(0.05076718239027029*s) + 5.792596760796093/exp(0.016285467000648854*s)) 
    corr_ij(3,7)=0
    corr_ij(3,6)=0
    corr_ij(3,5)=0
    corr_ij(3,4)=0
    corr_ij(3,3)=(-0.15609785880208227/exp(0.05076718239027029*s) + 31.18882918422289/exp(0.016285467000648854*s)) 
    corr_ij(3,2)=0
    corr_ij(3,1)=0
    corr_ij(2,10)=0

    corr_ij(2,9)=0
    corr_ij(2,8)=0
    corr_ij(2,7)=((1.6172201305728584*Cos(0.029341409726871926*s))/exp(0.04435489221745234*s) + (0.37871789179790255*C&
         &os(0.11425932545092894*s))/exp(0.019700737327669783*s) + (1.2889451151208258*Sin(0.02934140972687192&
         &6*s))/exp(0.04435489221745234*s) + (1.4228849217537705*Sin(0.11425932545092894*s))/exp(0.019700737327669783*s)) 
    corr_ij(2,6)=((-1.2889451151208255*Cos(0.029341409726871926*s))/exp(0.04435489221745234*s) - (1.4228849217537702*C&
         &os(0.11425932545092894*s))/exp(0.019700737327669783*s) + (1.6172201305728586*Sin(0.02934140972687192&
         &6*s))/exp(0.04435489221745234*s) + (0.3787178917979035*Sin(0.11425932545092894*s))/exp(0.019700737327669783*s)) 
    corr_ij(2,5)=0
    corr_ij(2,4)=0
    corr_ij(2,3)=0
    corr_ij(2,2)=((0.1789135645266575*Cos(0.029341409726871926*s))/exp(0.04435489221745234*s) + (26.817024457844113*Co&
         &s(0.11425932545092894*s))/exp(0.019700737327669783*s) - (0.4268927977731004*Sin(0.029341409726871926&
         &*s))/exp(0.04435489221745234*s) + (0.4268927977730982*Sin(0.11425932545092894*s))/exp(0.019700737327669783*s)) 
    corr_ij(2,1)=((0.4268927977731007*Cos(0.029341409726871926*s))/exp(0.04435489221745234*s) - (0.42689279777309963*C&
         &os(0.11425932545092894*s))/exp(0.019700737327669783*s) + (0.17891356452665746*Sin(0.0293414097268719&
         &26*s))/exp(0.04435489221745234*s) + (26.81702445784412*Sin(0.11425932545092894*s))/exp(0.019700737327669783*s)) 
    corr_ij(1,10)=0

    corr_ij(1,9)=0
    corr_ij(1,8)=0
    corr_ij(1,7)=((1.288945115120824*Cos(0.029341409726871926*s))/exp(0.04435489221745234*s) + (1.4228849217537711*Cos&
         &(0.11425932545092894*s))/exp(0.019700737327669783*s) - (1.617220130572856*Sin(0.029341409726871926*s&
         &))/exp(0.04435489221745234*s) - (0.3787178917979028*Sin(0.11425932545092894*s))/exp(0.019700737327669783*s)) 
    corr_ij(1,6)=((1.6172201305728564*Cos(0.029341409726871926*s))/exp(0.04435489221745234*s) + (0.37871789179790377*C&
         &os(0.11425932545092894*s))/exp(0.019700737327669783*s) + (1.2889451151208242*Sin(0.02934140972687192&
         &6*s))/exp(0.04435489221745234*s) + (1.4228849217537711*Sin(0.11425932545092894*s))/exp(0.019700737327669783*s)) 
    corr_ij(1,5)=0
    corr_ij(1,4)=0
    corr_ij(1,3)=0
    corr_ij(1,2)=((-0.4268927977731002*Cos(0.029341409726871926*s))/exp(0.04435489221745234*s) + (0.4268927977730981*C&
         &os(0.11425932545092894*s))/exp(0.019700737327669783*s) - (0.1789135645266573*Sin(0.02934140972687192&
         &6*s))/exp(0.04435489221745234*s) - (26.81702445784412*Sin(0.11425932545092894*s))/exp(0.019700737327669783*s)) 
    corr_ij(1,1)=((0.1789135645266574*Cos(0.029341409726871926*s))/exp(0.04435489221745234*s) + (26.817024457844113*Co&
         &s(0.11425932545092894*s))/exp(0.019700737327669783*s) - (0.42689279777310024*Sin(0.02934140972687192&
         &6*s))/exp(0.04435489221745234*s) + (0.4268927977730997*Sin(0.11425932545092894*s))/exp(0.019700737327669783*s)) 
    
    corr_ij=q_au**2*corr_ij

  END SUBROUTINE corrcomp_from_def

  !> Subroutine to compute the correlation of the unresolved variables \f$\langle Y \otimes Y^s \rangle\f$
  !> at time \f$s\f$ from the spline representation
  !> @param s time \f$s\f$ at which the correlation is computed
  SUBROUTINE corrcomp_from_spline(s)
    REAL(KIND=8), INTENT(IN) :: s
    REAL(KIND=8) :: y
    INTEGER :: i,j
    corr_ij=0.D0
    DO i=1,n_unres
       DO j=1,n_unres
          CALL splint(xa,ya(i,j,:),y2(i,j,:),nspl,s,y)
          corr_ij(i,j)=y
       END DO
    END DO
  END SUBROUTINE corrcomp_from_spline

  !> Routine to compute the spline representation parameters
  SUBROUTINE splint(xa,ya,y2a,n,x,y)
    INTEGER, INTENT(IN) :: n
    REAL(KIND=8), INTENT(IN), DIMENSION(n) :: xa,y2a,ya
    REAL(KIND=8), INTENT(IN) :: x
    REAL(KIND=8), INTENT(OUT) :: y
    INTEGER :: k
    REAL(KIND=8) :: a,b,h
    if ((khi-klo.gt.1).or.(xa(klo).gt.x).or.(xa(khi).lt.x)) then
       if ((khi-klo.eq.1).and.(xa(klo).lt.x)) then
          khi=klo
          DO WHILE (xa(khi).lt.x)
             khi=khi+1
          END DO
          klo=khi-1
       else
          khi=n
          klo=1
          DO WHILE (khi-klo.gt.1)
             k=(khi+klo)/2
             if(xa(k).gt.x)then
                khi=k
             else
                klo=k
             endif
          END DO
       end if
    !    print*, "search",x,khi-klo,xa(klo),xa(khi)
    ! else
    !    print*, "ok",x,khi-klo,xa(klo),xa(khi)
    endif
    h=xa(khi)-xa(klo)
    if (h.eq.0.) stop 'bad xa input in splint'
    a=(xa(khi)-x)/h
    b=(x-xa(klo))/h
    y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
    return
  END SUBROUTINE splint
  
  !> Exponential fit function
  !> @param s time \f$s\f$ at which the function is evaluated
  !> @param p vector holding the coefficients of the fit function
  FUNCTION fs(s,p)
    REAL(KIND=8), INTENT(IN) :: s
    REAL(KIND=8), DIMENSION(5), INTENT(IN) :: p
    REAL(KIND=8) :: fs
    fs=p(1)*exp(-s/p(2))*cos(p(3)*s+p(4))
    RETURN 
  END FUNCTION fs

  !> Subroutine to compute the correlation of the unresolved variables \f$\langle Y \otimes Y^s \rangle\f$
  !> at time \f$s\f$ from the exponential representation
  !> @param s time \f$s\f$ at which the correlation is computed
  SUBROUTINE corrcomp_from_fit(s)
    REAL(KIND=8), INTENT(IN) :: s
    REAL(KIND=8) :: y
    INTEGER :: i,j

    corr_ij=0.D0
    DO i=1,n_unres
       DO j=1,n_unres
          corr_ij(i,j)=fs(s,ya(i,j,:))
       END DO
    END DO
  END SUBROUTINE corrcomp_from_fit


END MODULE corrmod
