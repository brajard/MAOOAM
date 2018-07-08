
! sqrt_mod.f90
!
!> Utility module with various routine to compute matrix square root.
!
!> @copyright                                                               
!> 2018 Jonathan Demaeyer.
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------!
!                                                                           
!> @remark                                                                 
!> Mainly based on the numerical recipes and from:
!> Edvin Deadman, Nicholas J. Higham, Rui Ralha (2013)
!> "Blocked Schur Algorithms for Computing the Matrix Square Root",
!> Lecture Notes in Computer Science, 7782. pp. 171-182.
!
!---------------------------------------------------------------------------!

MODULE sqrt_mod
  USE params, only:ndim
  USE util, only: triu,diag,cdiag,floordiv
  IMPLICIT NONE
  
  PRIVATE

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: work
  ! COMPLEX(KIND=16), DIMENSION(:), ALLOCATABLE :: zwork
  
  INTEGER :: lwork

  REAL(KIND=8), PARAMETER :: real_eps = 2.2204460492503131e-16

  PUBLIC :: sqrtm,init_sqrt,chol,sqrtm_svd

CONTAINS

  SUBROUTINE init_sqrt
    INTEGER :: AllocStat
    lwork=10
    lwork=ndim*lwork

    ! print*, lwork
    
    IF (ALLOCATED(work)) THEN
       DEALLOCATE(work, STAT=AllocStat)
       IF (AllocStat /= 0) STOP "*** Deallocation problem ! ***"
    ENDIF
    ALLOCATE(work(lwork), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    ! ALLOCATE(zwork(lwork), STAT=AllocStat)
    ! IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
  END SUBROUTINE init_sqrt
  
  !> Routine to compute a real square-root of a matrix
  !> @param A Matrix whose square root to evaluate.
  !> @param sqA Square root of `A`.
  !> @param info Information code returned by the Lapack routines.
  !> @param info_triu Information code returned by the triangular matrix Lapack routines.
  !> @param bs Optional blocksize specification variable.
  SUBROUTINE sqrtm(A,sqA,info,info_triu,bs)
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: A
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: sqA
    INTEGER, INTENT(IN), OPTIONAL :: bs
    INTEGER, INTENT(OUT) :: info,info_triu
    REAL(KIND=8), DIMENSION(SIZE(A,1),SIZE(A,1)) :: T,Z,R
    COMPLEX(KIND=16), DIMENSION(SIZE(A,1),SIZE(A,1)) :: Tz,Zz,Rz
    REAL(KIND=8), DIMENSION(SIZE(A,1)) :: wr,wi
    LOGICAL, DIMENSION(SIZE(A,1)) :: bwork
    LOGICAL :: selectev
    INTEGER :: n
    INTEGER :: sdim=0
    n=SIZE(A,1)
    T=A
    ! print*, n,size(work,1)
    CALL DGEES('v','n',selectev,n,T,n,sdim,wr,wi,Z,n,work,lwork,bwork,info)
    ! print*, 'Z'
    ! CALL printmat(Z)
    ! print*, 'T'
    ! CALL printmat(T)
    ! CALL DGEES('V','N',SIZE(T,1),T,SIZE(T,1),0,wr,wi,Z,SIZE(Z,1),work,lwork,info)
    ! print*, info
    CALL triu(T,R)
    IF (ANY(T /= R)) THEN
       ! print*, 'T'
       ! CALL printmat(T)
       ! print*, 'Z'
       ! CALL printmat(Z)
       CALL rsf2csf(T,Z,Tz,Zz)
       ! print*, 'Tz'
       ! CALL printmat(dble(Tz))
       ! print*, 'iTz'
       ! CALL printmat(dble(aimag(Tz)))
       ! print*, 'Zz'
       ! CALL printmat(dble(Zz))
       ! print*, 'iZz'
       ! CALL printmat(dble(aimag(Zz)))
       IF (PRESENT(bs)) THEN
          CALL csqrtm_triu(Tz,Rz,info_triu,bs)
       ELSE
          CALL csqrtm_triu(Tz,Rz,info_triu)
       END IF
       Rz=matmul(Zz,matmul(Rz,conjg(transpose(Zz))))
       ! print*, 'sqAz'
       ! CALL printmat(dble(Rz))
       ! print*, 'isqAz'
       ! CALL printmat(dble(aimag(Rz)))
       sqA=dble(Rz)
    ELSE
       IF (PRESENT(bs)) THEN
          CALL sqrtm_triu(T,R,info_triu,bs)
       ELSE
          CALL sqrtm_triu(T,R,info_triu)
       END IF
       sqA=matmul(Z,matmul(R,transpose(Z)))
    ENDIF

  END SUBROUTINE sqrtm
  
  FUNCTION selectev(a,b)
    REAL(KIND=8) :: a,b
    LOGICAL selectev
    selectev=.false.
    ! IF (a>b) selectev=.true.
    RETURN
  END FUNCTION selectev
    

  SUBROUTINE sqrtm_triu(A,sqA,info,bs)
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: A
    INTEGER, INTENT(IN), OPTIONAL :: bs
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: sqA
    INTEGER, INTENT(OUT) :: info
    REAL(KIND=8), DIMENSION(SIZE(A,1)) :: A_diag
    REAL(KIND=8), DIMENSION(SIZE(A,1),SIZE(A,1)) :: R,Sm,Rii,Rjj
    INTEGER, DIMENSION(2*SIZE(A,1),2) :: start_stop_pairs
    REAL(KIND=8) :: s,denom,scale
    INTEGER :: i,j,k,start,n,sstop,m
    INTEGER :: istart,istop,jstart,jstop
    INTEGER :: nblocks,blocksize
    INTEGER :: bsmall,blarge,nlarge,nsmall

    blocksize=64
    IF (PRESENT(bs)) blocksize=bs
    n=SIZE(A,1)
    ! print*,  blocksize

    CALL diag(A,A_diag)
    R=0.D0
    DO i=1,n
       R(i,i)=sqrt(A_diag(i))
    ENDDO
    

    nblocks=max(floordiv(n,blocksize),1)
    bsmall=floordiv(n,nblocks)
    nlarge=mod(n,nblocks)
    blarge=bsmall+1
    nsmall=nblocks-nlarge
    IF (nsmall*bsmall + nlarge*blarge /= n) STOP 'Sqrtm: Internal inconsistency'

    ! print*, nblocks,bsmall,nsmall,blarge,nlarge

    start=1
    DO i=1,nsmall
       start_stop_pairs(i,1)=start
       start_stop_pairs(i,2)=start+bsmall-1
       start=start+bsmall
    ENDDO
    DO i=nsmall+1,nsmall+nlarge
       start_stop_pairs(i,1)=start
       start_stop_pairs(i,2)=start+blarge-1
       start=start+blarge
    ENDDO
    
    ! DO i=1,SIZE(start_stop_pairs,1)
    !    print*, i
    !    print*, start_stop_pairs(i,1),start_stop_pairs(i,2)
    ! END DO
    
    DO k=1,nsmall+nlarge
       start=start_stop_pairs(k,1)
       sstop=start_stop_pairs(k,2)
       DO j=start,sstop
          DO i=j-1,start,-1
             s=0.D0
             IF (j-i>1) s= dot_product(R(i,i+1:j-1),R(i+1:j-1,j))
             denom= R(i,i)+R(j,j)
             IF (denom==0.D0) STOP 'Sqrtm: Failed to find the matrix square root'
             R(i,j)=(A(i,j)-s)/denom
          END DO
       END DO
    END DO

    ! print*, 'R'
    ! CALL printmat(R)
    
    DO j=1,nblocks
       jstart=start_stop_pairs(j,1)
       jstop=start_stop_pairs(j,2)
       DO i=j-1,1,-1
          istart=start_stop_pairs(i,1)
          istop=start_stop_pairs(i,2)
          Sm=0.D0
          Sm(istart:istop,jstart:jstop)=A(istart:istop,jstart:jstop)
          IF (j-i>1) Sm(istart:istop,jstart:jstop) = Sm(istart:istop&
               &,jstart:jstop) - matmul(R(istart:istop,istop:jstart)&
               &,R(istop:jstart,jstart:jstop))
          Rii=0.D0
          Rii = R(istart:istop, istart:istop)
          Rjj=0.D0
          Rjj = R(jstart:jstop, jstart:jstop)
          m=istop-istart+1
          n=jstop-jstart+1
          k=1
          ! print*, m,n
          ! print*, istart,istop
          ! print*, jstart,jstop

          ! print*, 'Rii',Rii(istart:istop, istart:istop)
          ! print*, 'Rjj',Rjj(jstart:jstop,jstart:jstop)
          ! print*, 'Sm',Sm(istart:istop,jstart:jstop)

          CALL dtrsyl('N','N',k,m,n,Rii(istart:istop, istart:istop),m&
               &,Rjj(jstart:jstop,jstart:jstop),n,Sm(istart:istop&
               &,jstart:jstop),m,scale,info)
          R(istart:istop,jstart:jstop)=Sm(istart:istop,jstart:jstop)*scale
       ENDDO
    ENDDO
    sqA=R
  END SUBROUTINE sqrtm_triu

  SUBROUTINE csqrtm_triu(A,sqA,info,bs)
    COMPLEX(KIND=16), DIMENSION(:,:), INTENT(IN) :: A
    INTEGER, INTENT(IN), OPTIONAL :: bs
    COMPLEX(KIND=16), DIMENSION(:,:), INTENT(OUT) :: sqA
    INTEGER, INTENT(OUT) :: info
    COMPLEX(KIND=16), DIMENSION(SIZE(A,1)) :: A_diag
    COMPLEX(KIND=16), DIMENSION(SIZE(A,1),SIZE(A,1)) :: R,Sm,Rii,Rjj
    INTEGER, DIMENSION(2*SIZE(A,1),2) :: start_stop_pairs
    COMPLEX(KIND=16) :: s,denom,scale
    INTEGER :: i,j,k,start,n,sstop,m
    INTEGER :: istart,istop,jstart,jstop
    INTEGER :: nblocks,blocksize
    INTEGER :: bsmall,blarge,nlarge,nsmall

    blocksize=64
    IF (PRESENT(bs)) blocksize=bs
    n=SIZE(A,1)
    ! print*,  blocksize

    CALL cdiag(A,A_diag)
    R=0.D0
    DO i=1,n
       R(i,i)=sqrt(A_diag(i))
    ENDDO
    

    nblocks=max(floordiv(n,blocksize),1)
    bsmall=floordiv(n,nblocks)
    nlarge=mod(n,nblocks)
    blarge=bsmall+1
    nsmall=nblocks-nlarge
    IF (nsmall*bsmall + nlarge*blarge /= n) STOP 'Sqrtm: Internal inconsistency'

    ! print*, nblocks,bsmall,nsmall,blarge,nlarge

    start=1
    DO i=1,nsmall
       start_stop_pairs(i,1)=start
       start_stop_pairs(i,2)=start+bsmall-1
       start=start+bsmall
    ENDDO
    DO i=nsmall+1,nsmall+nlarge
       start_stop_pairs(i,1)=start
       start_stop_pairs(i,2)=start+blarge-1
       start=start+blarge
    ENDDO
    
    ! DO i=1,SIZE(start_stop_pairs,1)
    !    print*, i
    !    print*, start_stop_pairs(i,1),start_stop_pairs(i,2)
    ! END DO
    
    DO k=1,nsmall+nlarge
       start=start_stop_pairs(k,1)
       sstop=start_stop_pairs(k,2)
       DO j=start,sstop
          DO i=j-1,start,-1
             s=0.D0
             IF (j-i>1) s= dot_product(R(i,i+1:j-1),R(i+1:j-1,j))
             denom= R(i,i)+R(j,j)
             IF (denom==0.D0) STOP 'Sqrtm: Failed to find the matrix square root'
             R(i,j)=(A(i,j)-s)/denom
          END DO
       END DO
    END DO

    ! print*, 'R'
    ! CALL printmat(R)
    
    DO j=1,nblocks
       jstart=start_stop_pairs(j,1)
       jstop=start_stop_pairs(j,2)
       DO i=j-1,1,-1
          istart=start_stop_pairs(i,1)
          istop=start_stop_pairs(i,2)
          Sm=0.D0
          Sm(istart:istop,jstart:jstop)=A(istart:istop,jstart:jstop)
          IF (j-i>1) Sm(istart:istop,jstart:jstop) = Sm(istart:istop&
               &,jstart:jstop) - matmul(R(istart:istop,istop:jstart)&
               &,R(istop:jstart,jstart:jstop))
          Rii=0.D0
          Rii = R(istart:istop, istart:istop)
          Rjj=0.D0
          Rjj = R(jstart:jstop, jstart:jstop)
          m=istop-istart+1
          n=jstop-jstart+1
          k=1
          ! print*, m,n
          ! print*, istart,istop
          ! print*, jstart,jstop

          ! print*, 'Rii',Rii(istart:istop, istart:istop)
          ! print*, 'Rjj',Rjj(jstart:jstop,jstart:jstop)
          ! print*, 'Sm',Sm(istart:istop,jstart:jstop)

          CALL ztrsyl('N','N',k,m,n,Rii(istart:istop, istart:istop),m&
               &,Rjj(jstart:jstop,jstart:jstop),n,Sm(istart:istop&
               &,jstart:jstop),m,scale,info)
          R(istart:istop,jstart:jstop)=Sm(istart:istop,jstart:jstop)*scale
       ENDDO
    ENDDO
    sqA=R
  END SUBROUTINE csqrtm_triu

  SUBROUTINE rsf2csf(T,Z,Tz,Zz)
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: T,Z
    COMPLEX(KIND=16), DIMENSION(:,:), INTENT(OUT) :: Tz,Zz
    INTEGER, PARAMETER :: nb=2
    COMPLEX(KIND=16), DIMENSION(nb) :: w
    !COMPLEX(KIND=16), DIMENSION(nb,nb) :: vl,vr
    COMPLEX(KIND=16) :: r,c,s,mu
    COMPLEX(KIND=16), DIMENSION(nb,nb) :: G,Gc
    INTEGER :: N,m!,info
    !REAL(KIND=8), DIMENSION(2*nb) :: rwork
    !REAL(KIND=8), DIMENSION(2*nb) :: ztwork

    ! print*, lwork
    Tz=cmplx(T,KIND=16)
    Zz=cmplx(Z,KIND=16)
    N=SIZE(T,1)
    DO m=N,2,-1
       IF (abs(Tz(m,m-1)) > real_eps*(abs(Tz(m-1,m-1)) + abs(Tz(m,m)))) THEN
          G=Tz(m-1:m,m-1:m)
          ! CALL printmat(dble(G))
          ! CALL zgeev('N','N',nb,G,nb,w,vl,nb,vr,nb,ztwork,2*nb,rwork,info)
          ! CALL cprintmat(G)
          ! print*, m,w,info
          s=G(1,1)+G(2,2)
          c=G(1,1)*G(2,2)-G(1,2)*G(2,1)
          w(1)=s/2+sqrt(s**2/4-c)
          mu=w(1)-Tz(m,m)
          r=sqrt(mu*conjg(mu)+Tz(m,m-1)*conjg(Tz(m,m-1)))
          c=mu/r
          s=Tz(m,m-1)/r
          G(1,1)=conjg(c)
          G(1,2)=s
          G(2,1)=-s
          G(2,2)=c
          Gc=conjg(transpose(G))
          Tz(m-1:m,m-1:N)=matmul(G,Tz(m-1:m,m-1:N))
          Tz(1:m,m-1:m)=matmul(Tz(1:m,m-1:m),Gc)
          Zz(:,m-1:m)=matmul(Zz(:,m-1:m),Gc)
       END IF
       Tz(m,m-1)=cmplx(0.D0,KIND=16)
    END DO
  END SUBROUTINE rsf2csf
  
  !> Routine to perform a Cholesky decomposition
  !> @param A Matrix whose decomposition is evaluated.
  !> @param sqA Cholesky decomposition of `A`.
  !> @param info Information code returned by the Lapack routines.
  SUBROUTINE chol(A,sqA,info)
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: A
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: sqA
    INTEGER, INTENT(OUT) :: info

    sqA=A
    CALL DPOTRF('L',SIZE(sqA,1),sqA,SIZE(sqA,1),info)
  END SUBROUTINE chol

  !> Routine to compute a real square-root of a matrix via a SVD decomposition
  !> @param A Matrix whose square root to evaluate.
  !> @param sqA Square root of `A`.
  !> @param info Information code returned by the Lapack routines.
  !> @param info_triu Not used (present for compatibility).
  !> @param bs Not used (present for compatibility).
  SUBROUTINE sqrtm_svd(A,sqA,info,info_triu,bs)
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: A
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: sqA
    INTEGER, INTENT(IN), OPTIONAL :: bs
    INTEGER, INTENT(OUT) :: info,info_triu
    REAL(KIND=8), DIMENSION(SIZE(A,1)) :: S
    REAL(KIND=8), DIMENSION(SIZE(A,1),SIZE(A,1)) :: Sq,U,VT
    INTEGER :: i,n
    
    sqA=A
    n=SIZE(sqA,1)
    CALL DGESVD('A','A',n,n,sqA,n,S,U,n,VT,n,work,lwork,info)
    Sq=0.D0
    DO i=1,n
       Sq(i,i)=sqrt(S(i))
    ENDDO
    sqA=matmul(U,matmul(Sq,VT))
  END SUBROUTINE sqrtm_svd

END MODULE sqrt_mod

  
    
  
  
