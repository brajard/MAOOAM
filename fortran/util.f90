
! util.f90
!
!>  Utility module
!
!> @copyright                                                               
!> 2018 Lesley De Cruz & Jonathan Demaeyer.
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------!

MODULE util
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: str,rstr,init_random_seed,init_one,isin,piksrt
  PUBLIC :: mat_trace, mat_contract, choldc
  PUBLIC :: printmat,cprintmat,invmat,ireduce,reduce,floordiv,triu,diag,cdiag
  PUBLIC :: vector_outer

CONTAINS
    
  ! SUBROUTINE scalar_allocate(x)
  !   INTEGER :: AllocStat
  !   IF (.NOT. ALLOCATED(x)) THEN 
  !      ALLOCATE(x, STAT=AllocStat)
  
  !> Convert an integer to string.
  CHARACTER(len=20) FUNCTION str(k) 
    INTEGER, INTENT(IN) :: k
    WRITE (str, *) k
    str = ADJUSTL(str)
  END FUNCTION str

  !> Convert a real to string with a given format
  CHARACTER(len=40) FUNCTION rstr(x,fm) 
    REAL(KIND=8), INTENT(IN) :: x
    CHARACTER(len=20), INTENT(IN) :: fm
    WRITE (rstr, TRIM(ADJUSTL(fm))) x
    rstr = ADJUSTL(rstr)
  END FUNCTION rstr

  !> Determine if a character is in a string and where
  !> @remark: return positions in a vector if found and 0 vector if not found
  FUNCTION isin(c,s)
    CHARACTER, INTENT(IN) :: c
    CHARACTER, DIMENSION(:), INTENT(IN) :: s
    INTEGER, DIMENSION(size(s)) :: isin
    INTEGER :: i,j

    isin=0
    j=0
    DO i=size(s),1,-1
       IF (c==s(i)) THEN
          j=j+1
          isin(j)=i
       END IF
    END DO
  END FUNCTION isin

  !> Random generator initialization routine  
  SUBROUTINE init_random_seed()
    USE iso_fortran_env, only: int64
    USE IFPORT !, only: getpid
    IMPLICIT NONE
    INTEGER, ALLOCATABLE :: seed_loc(:)
    INTEGER :: i, n, un, istat, dt(8), pid
    INTEGER(int64) :: t

    CALL random_seed(size = n)
    ALLOCATE(seed_loc(n))
    ! First try IF the OS provides a random number generator
    OPEN(newunit=un, file="/dev/urandom", access="stream", &
         form="unformatted", action="read", status="old", iostat=istat)
    IF (istat == 0) THEN
       READ(un) seed_loc
       CLOSE(un)
    ELSE
       ! Fallback to XOR:ing the current time and pid. The PID is
       ! useful in case one launches multiple instances of the same
       ! program in parallel.
       CALL system_clock(t)
       IF (t == 0) THEN
          CALL date_and_time(values=dt)
          t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
               + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
               + dt(3) * 24_int64 * 60 * 60 * 1000 &
               + dt(5) * 60 * 60 * 1000 &
               + dt(6) * 60 * 1000 + dt(7) * 1000 &
               + dt(8)
       END IF
       pid = getpid()
       t = ieor(t, int(pid, kind(t)))
       DO i = 1, n
          seed_loc(i) = lcg(t)
       END DO
    END IF
    CALL random_seed(put=seed_loc)
  contains
    ! This simple PRNG might not be good enough for real work, but is
    ! sufficient for seeding a better PRNG.
    FUNCTION lcg(s)
      integer :: lcg
      integer(int64) :: s
      IF (s == 0) THEN
         s = 104729
      ELSE
         s = mod(s, 4294967296_int64)
      END IF
      s = mod(s * 279470273_int64, 4294967291_int64)
      lcg = int(mod(s, int(huge(0), int64)), kind(0))
    END FUNCTION lcg
  END SUBROUTINE init_random_seed

  !> Simple card player sorting function
  SUBROUTINE piksrt(k,arr,par)
    INTEGER, INTENT(IN) :: k
    INTEGER, DIMENSION(k), INTENT(INOUT) :: arr
    INTEGER, INTENT(OUT) :: par
    INTEGER :: i,j,a,b

    par=1

    DO j=2,k
       a=arr(j)
       DO i=j-1,1,-1
          if(arr(i).le.a) EXIT
          arr(i+1)=arr(i)
          par=-par
       END DO
       arr(i+1)=a
    ENDDO
    RETURN
  END SUBROUTINE piksrt

  !> Initialize a square matrix A as a unit matrix
  SUBROUTINE init_one(A)
     REAL(KIND=8), DIMENSION(:,:),INTENT(INOUT) :: A
     INTEGER :: i,n
     n=size(A,1)
     A=0.0d0
     DO i=1,n
       A(i,i)=1.0d0
     END DO

  END SUBROUTINE init_one

  FUNCTION mat_trace(A)
    REAL(KIND=8), DIMENSION(:,:) :: A
    REAL(KIND=8) :: mat_trace
    INTEGER :: i,n
    n=size(A,1)
    mat_trace=0.D0
    DO i=1,n
       mat_trace=mat_trace+A(i,i)
    END DO
    RETURN
  END FUNCTION mat_trace

  FUNCTION mat_contract(A,B)
    REAL(KIND=8), DIMENSION(:,:) :: A,B
    REAL(KIND=8) :: mat_contract
    INTEGER :: i,j,n
    n=size(A,1)
    mat_contract=0.D0
    DO i=1,n
       DO j=1,n
          mat_contract=mat_contract+A(i,j)*B(i,j)
       END DO
    ENDDO
    RETURN
  END FUNCTION mat_contract

  SUBROUTINE choldc(a,p)
    REAL(KIND=8), DIMENSION(:,:) :: a
    REAL(KIND=8), DIMENSION(:) :: p
    INTEGER :: n
    INTEGER :: i,j,k
    REAL(KIND=8) :: sum
    n=size(a,1)
    DO i=1,n
       DO j=i,n
          sum=a(i,j)
          DO k=i-1,1,-1
             sum=sum-a(i,k)*a(j,k)
          END DO
          IF (i.eq.j) THEN
             IF (sum.le.0.) stop 'choldc failed'
             p(i)=sqrt(sum)
          ELSE
             a(j,i)=sum/p(i)
          ENDIF
       END DO
    END DO
    RETURN
  END SUBROUTINE choldc

  SUBROUTINE printmat(A) ! to be moved to util
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: A
    INTEGER :: i
    
    DO i=1,SIZE(A,1)
       print*, A(i,:)
    END DO
  END SUBROUTINE printmat

  SUBROUTINE cprintmat(A) ! to be moved to util
    COMPLEX(KIND=16), DIMENSION(:,:), INTENT(IN) :: A
    INTEGER :: i
    
    DO i=1,SIZE(A,1)
       print*, A(i,:)
    END DO
  END SUBROUTINE cprintmat

  FUNCTION invmat(A) RESULT(Ainv)
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: A
    REAL(KIND=8), DIMENSION(SIZE(A,1),SIZE(A,2)) :: Ainv

    REAL(KIND=8), DIMENSION(SIZE(A,1)) :: work  ! work array for LAPACK
    INTEGER, DIMENSION(SIZE(A,1)) :: ipiv   ! pivot indices
    INTEGER :: n, info

    ! Store A in Ainv to prevent it from being overwritten by LAPACK
    Ainv = A
    n = size(A,1)

    ! DGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    CALL DGETRF(n, n, Ainv, n, ipiv, info)

    IF (info /= 0) THEN
       STOP 'Matrix is numerically singular!'
    ENDIF

    ! DGETRI computes the inverse of a matrix using the LU factorization
    ! computed by DGETRF.
    CALL DGETRI(n, Ainv, n, ipiv, work, n, info)

    IF (info /= 0) THEN
       STOP 'Matrix inversion failed!'
    ENDIF
    END FUNCTION invmat 

    SUBROUTINE triu(A,T)
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: A
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: T
    INTEGER i,j
    T=0.D0
    DO i=1,SIZE(A,1)
       DO j=i,SIZE(A,1)
          T(i,j)=A(i,j)
       END DO
    END DO
  END SUBROUTINE triu
    
  SUBROUTINE diag(A,d)
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: A
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: d
    INTEGER :: i
    
    DO i=1,SIZE(A,1)
       d(i)=A(i,i)
    END DO
  END SUBROUTINE diag

  SUBROUTINE cdiag(A,d)
    COMPLEX(KIND=16), DIMENSION(:,:), INTENT(IN) :: A
    COMPLEX(KIND=16), DIMENSION(:), INTENT(OUT) :: d
    INTEGER :: i
    
    DO i=1,SIZE(A,1)
       d(i)=A(i,i)
    END DO
  END SUBROUTINE cdiag

       
  FUNCTION floordiv(i,j)
    INTEGER :: i,j,floordiv
    floordiv=int(floor(real(i)/real(j)))
    RETURN
  END FUNCTION floordiv

  SUBROUTINE reduce(A,Ared,n,ind,rind)
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: A
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: Ared
    INTEGER, INTENT(OUT) :: n
    INTEGER, DIMENSION(:), INTENT(OUT) :: ind,rind
    LOGICAL, DIMENSION(SIZE(A,1)) :: sel
    INTEGER :: i,j

    ind=0
    rind=0
    sel=.FALSE.
    n=0
    DO i=1,SIZE(A,1)
       IF (ANY(A(i,:)/=0)) THEN
          n=n+1
          sel(i)=.TRUE.
          ind(n)=i
          rind(i)=n
       ENDIF
    END DO
    Ared=0.D0
    DO i=1,SIZE(A,1)
       DO j=1,SIZE(A,1)
          IF (sel(i).and.sel(j)) Ared(rind(i),rind(j))=A(i,j)
       ENDDO
    ENDDO
  END SUBROUTINE reduce

  SUBROUTINE ireduce(A,Ared,n,ind,rind)
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: A
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: Ared
    INTEGER, INTENT(IN) :: n
    INTEGER, DIMENSION(:), INTENT(IN) :: ind,rind
    INTEGER :: i,j
    A=0.D0
    DO i=1,n
       DO j=1,n
          A(ind(i),ind(j))=Ared(i,j)
       END DO
    END DO
  END SUBROUTINE ireduce

  SUBROUTINE vector_outer(u,v,A)
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u,v
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: A
    INTEGER :: i,j

    A=0.D0
    DO i=1,SIZE(u)
       DO j=1,SIZE(v)
          A(i,j)=u(i)*v(j)
       ENDDO
    ENDDO
  END SUBROUTINE vector_outer

END MODULE util
