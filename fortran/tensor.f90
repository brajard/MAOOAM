
! tensor.f90
!
!>  Tensor utility module
!
!> @copyright                                                               
!> 2015-2018 Lesley De Cruz & Jonathan Demaeyer.
!> See LICENSE.txt for license information.
!
!---------------------------------------------------------------------------!
!
!> @remark
!> coolist is a type and also means "coordinate list"
!
!---------------------------------------------------------------------------!


MODULE tensor
  USE params, only: ndim
  IMPLICIT NONE

  PRIVATE

  !> Coordinate list element type. Elementary elements of the sparse tensors.
  TYPE :: coolist_elem
     INTEGER :: j !< Index \f$j\f$ of the element
     INTEGER :: k !< Index \f$k\f$ of the element
     REAL(KIND=8) :: v !< Value of the element
  END TYPE coolist_elem

  !> 4d coordinate list element type. Elementary elements of the 4d sparse tensors.
  TYPE :: coolist_elem4
     INTEGER :: j,k,l
     REAL(KIND=8) :: v
  END TYPE coolist_elem4

  !> Coordinate list. Type used to represent the sparse tensor.
  TYPE, PUBLIC :: coolist
     TYPE(coolist_elem), DIMENSION(:), ALLOCATABLE :: elems !< Lists of elements tensor::coolist_elem
     INTEGER :: nelems = 0 !< Number of elements in the list.
  END TYPE coolist

  !> 4d coordinate list. Type used to represent the rank-4 sparse tensor.
  TYPE, PUBLIC :: coolist4
     TYPE(coolist_elem4), DIMENSION(:), ALLOCATABLE :: elems
     INTEGER :: nelems = 0
  END TYPE coolist4
  
  !> Parameter to test the equality with zero.
  REAL(KIND=8), PARAMETER :: real_eps = 2.2204460492503131e-16

  PUBLIC :: sparse_mul3,sparse_mul2,copy_coo,mat_to_coo,jsparse_mul,jsparse_mul_mat,simplify
  PUBLIC :: add_elem,add_to_tensor,print_tensor,load_tensor_from_file,write_tensor_to_file,add_check

  PUBLIC :: sparse_mul4,sparse_mul2_j,sparse_mul2_k
  PUBLIC :: tensor_to_coo,tensor4_to_coo4
  PUBLIC :: print_tensor4
  PUBLIC :: coo_to_mat_ik,coo_to_mat_ij,coo_to_mat_i,sparse_mul4_with_mat_jl,matc_to_coo,coo_to_mat_j
  PUBLIC :: add_matc_to_tensor,coo_to_vec_jk,add_vec_jk_to_tensor,add_vec_ikl_to_tensor4_perm!,tensor_ij_perm
  PUBLIC :: add_matc_to_tensor4,load_tensor4_from_file,write_tensor4_to_file
  PUBLIC :: tensor_empty,tensor4_empty
  PUBLIC :: sparse_mul4_mat,sparse_mul3_mat,sparse_mul3_with_mat,add_vec_ikl_to_tensor4
  PUBLIC :: sparse_mul4_with_mat_kl,add_vec_ijk_to_tensor4,scal_mul_coo

CONTAINS
    
  !> Routine to copy a coolist.
  !> @param src Source coolist
  !> @param dst Destination coolist
  !> @remark The destination tensor have to be an empty tensor, i.e. with unallocated list of elements and nelems set to 0.
  SUBROUTINE copy_coo(src,dst)
    TYPE(coolist), DIMENSION(ndim), INTENT(IN) :: src
    TYPE(coolist), DIMENSION(ndim), INTENT(OUT) :: dst
    INTEGER :: i,j,AllocStat
    
    DO i=1,ndim
       IF (dst(i)%nelems/=0) STOP "*** copy_coo : Destination coolist not empty ! ***"
       ALLOCATE(dst(i)%elems(src(i)%nelems), STAT=AllocStat)
       IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
       DO j=1,src(i)%nelems
          dst(i)%elems(j)%j=src(i)%elems(j)%j
          dst(i)%elems(j)%k=src(i)%elems(j)%k
          dst(i)%elems(j)%v=src(i)%elems(j)%v
       ENDDO
       dst(i)%nelems=src(i)%nelems
    ENDDO
  END SUBROUTINE copy_coo

  !> Routine to convert a matrix to a tensor.
  !> @param src Source matrix
  !> @param dst Destination tensor
  !> @remark The destination tensor have to be an empty tensor, i.e. with unallocated list of elements and nelems set to 0.
  SUBROUTINE mat_to_coo(src,dst)
    REAL(KIND=8), DIMENSION(0:ndim,0:ndim), INTENT(IN) :: src
    TYPE(coolist), DIMENSION(ndim), INTENT(OUT) :: dst
    INTEGER :: i,j,n,AllocStat
    DO i=1,ndim
       n=0
       DO j=1,ndim
          IF (ABS(src(i,j))>real_eps) n=n+1
       ENDDO
       IF (n/=0) THEN
          IF (dst(i)%nelems/=0) STOP "*** mat_to_coo : Destination coolist not empty ! ***"
          ALLOCATE(dst(i)%elems(n), STAT=AllocStat)
          IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
          n=0
          DO j=1,ndim
             IF (ABS(src(i,j))>real_eps) THEN
                n=n+1
                dst(i)%elems(n)%j=j
                dst(i)%elems(n)%k=0
                dst(i)%elems(n)%v=src(i,j)
             ENDIF
          ENDDO
       ENDIF
       dst(i)%nelems=n
    ENDDO
  END SUBROUTINE mat_to_coo
  
  !> Sparse multiplication of a tensor with two vectors:  \f${\displaystyle \sum_{j,k=0}^{ndim}} \mathcal{T}_{i,j,k} \, a_j \,b_k\f$.
  !> @param coolist_ijk a coordinate list (sparse tensor) of which index
  !> 2 and 3 will be contracted.
  !> @param arr_j the vector to be contracted with index 2 of coolist_ijk
  !> @param arr_k the vector to be contracted with index 3 of coolist_ijk
  !> @param res vector (buffer) to store the result of the contraction
  !> @remark Note that it is NOT safe to pass `arr_j`/`arr_k` as a result buffer, 
  !> as this operation does multiple passes.
  SUBROUTINE sparse_mul3(coolist_ijk, arr_j, arr_k, res)
    TYPE(coolist), DIMENSION(ndim), INTENT(IN):: coolist_ijk
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(IN)  :: arr_j, arr_k
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(OUT) :: res
    INTEGER :: i,j,k,n
    res=0.D0
    DO i=1,ndim
       DO n=1,coolist_ijk(i)%nelems
         j=coolist_ijk(i)%elems(n)%j
         k=coolist_ijk(i)%elems(n)%k
         res(i) = res(i) + coolist_ijk(i)%elems(n)%v * arr_j(j)*arr_k(k)
      END DO
   END DO
  END SUBROUTINE sparse_mul3

  !> Sparse multiplication of two tensors to determine the Jacobian:
  !> \f[J_{i,j} = {\displaystyle \sum_{k=0}^{ndim}} \left( \mathcal{T}_{i,j,k} + \mathcal{T}_{i,k,j} \right) \, a_k.\f]
  !> It's implemented slightly differently: for every \f$\mathcal{T}_{i,j,k}\f$, we add to \f$J_{i,j}\f$ as follows:
  !> \f[J_{i,j} = J_{i,j} + \mathcal{T}_{i,j,k} \, a_k \\ J_{i,k} = J_{i,k} + \mathcal{T}_{i,j,k} \, a_j\f]
  !> This version return a coolist (sparse tensor).
  !> @param coolist_ijk a coordinate list (sparse tensor) of which index 
  !> 2 or 3 will be contracted.
  !> @param arr_j the vector to be contracted with index 2 and then index 3 of ffi_coo_ijk
  !> @param jcoo_ij a coolist (sparse tensor) to store the result of the contraction
  SUBROUTINE jsparse_mul(coolist_ijk, arr_j, jcoo_ij)
    TYPE(coolist), DIMENSION(ndim), INTENT(IN):: coolist_ijk
    TYPE(coolist), DIMENSION(ndim), INTENT(OUT):: jcoo_ij
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(IN)  :: arr_j
    REAL(KIND=8) :: v
    INTEGER :: i,j,k,n,nj,AllocStat
    DO i=1,ndim
       IF (jcoo_ij(i)%nelems/=0) STOP "*** jsparse_mul : Destination coolist not empty ! ***"
       nj=2*coolist_ijk(i)%nelems
       ALLOCATE(jcoo_ij(i)%elems(nj), STAT=AllocStat)
       IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
       nj=0
       DO n=1,coolist_ijk(i)%nelems
          j=coolist_ijk(i)%elems(n)%j
          k=coolist_ijk(i)%elems(n)%k
          v=coolist_ijk(i)%elems(n)%v
          IF (j /=0) THEN
             nj=nj+1
             jcoo_ij(i)%elems(nj)%j=j
             jcoo_ij(i)%elems(nj)%k=0
             jcoo_ij(i)%elems(nj)%v=v*arr_j(k)
          END IF

          IF (k /=0) THEN
             nj=nj+1
             jcoo_ij(i)%elems(nj)%j=k
             jcoo_ij(i)%elems(nj)%k=0
             jcoo_ij(i)%elems(nj)%v=v*arr_j(j)
          END IF
       END DO
       jcoo_ij(i)%nelems=nj
    END DO
  END SUBROUTINE jsparse_mul

  !> Sparse multiplication of two tensors to determine the Jacobian:
  !> \f[J_{i,j} = {\displaystyle \sum_{k=0}^{ndim}} \left( \mathcal{T}_{i,j,k} + \mathcal{T}_{i,k,j} \right) \, a_k.\f]
  !> It's implemented slightly differently: for every \f$\mathcal{T}_{i,j,k}\f$, we add to \f$J_{i,j}\f$ as follows:
  !> \f[J_{i,j} = J_{i,j} + \mathcal{T}_{i,j,k} \, a_k \\ J_{i,k} = J_{i,k} + \mathcal{T}_{i,j,k} \, a_j\f]
  !> This version return a matrix.
  !> @param coolist_ijk a coordinate list (sparse tensor) of which index 
  !> 2 or 3 will be contracted.
  !> @param arr_j the vector to be contracted with index 2 and then index 3 of ffi_coo_ijk
  !> @param jcoo_ij a matrix to store the result of the contraction
  SUBROUTINE jsparse_mul_mat(coolist_ijk, arr_j, jcoo_ij)
    TYPE(coolist), DIMENSION(ndim), INTENT(IN):: coolist_ijk
    REAL(KIND=8), DIMENSION(ndim,ndim), INTENT(OUT):: jcoo_ij
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(IN)  :: arr_j
    REAL(KIND=8) :: v
    INTEGER :: i,j,k,n
    jcoo_ij=0.D0
    DO i=1,ndim
       DO n=1,coolist_ijk(i)%nelems
          j=coolist_ijk(i)%elems(n)%j
          k=coolist_ijk(i)%elems(n)%k
          v=coolist_ijk(i)%elems(n)%v
          IF (j /=0) jcoo_ij(i,j)=jcoo_ij(i,j)+v*arr_j(k)
          IF (k /=0) jcoo_ij(i,k)=jcoo_ij(i,k)+v*arr_j(j)
       END DO
    END DO
  END SUBROUTINE jsparse_mul_mat

  !> Sparse multiplication of a 2d sparse tensor with a vector:  \f${\displaystyle \sum_{j=0}^{ndim}} \mathcal{T}_{i,j,k} \, a_j \f$.
  !> @param coolist_ij a coordinate list (sparse tensor) of which index
  !> 2 will be contracted.
  !> @param arr_j the vector to be contracted with index 2 of coolist_ijk
  !> @param res vector (buffer) to store the result of the contraction
  !> @remark Note that it is NOT safe to pass `arr_j` as a result buffer, 
  !> as this operation does multiple passes.
  SUBROUTINE sparse_mul2(coolist_ij, arr_j, res)
    TYPE(coolist), DIMENSION(ndim), INTENT(IN):: coolist_ij
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(IN)  :: arr_j
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(OUT) :: res
    INTEGER :: i,j,n
    res=0.D0
    DO i=1,ndim
       DO n=1,coolist_ij(i)%nelems
         j=coolist_ij(i)%elems(n)%j
         res(i) = res(i) + coolist_ij(i)%elems(n)%v * arr_j(j)
      END DO
   END DO
 END SUBROUTINE sparse_mul2
 
 !> Routine to simplify a coolist (sparse tensor). For each index \f$i\f$, it upper triangularize the matrix
 !> \f[\mathcal{T}_{i,j,k} \qquad 0 \leq j,k \leq ndim.\f]
 !> @param tensor a coordinate list (sparse tensor) which will be simplified.
 SUBROUTINE simplify(tensor)
   TYPE(coolist), DIMENSION(ndim), INTENT(INOUT):: tensor
   INTEGER :: i,j,k
   INTEGER :: li,lii,liii,n
   DO i= 1,ndim
      n=tensor(i)%nelems
      DO li=n,2,-1
         j=tensor(i)%elems(li)%j
         k=tensor(i)%elems(li)%k
         DO lii=li-1,1,-1
            IF (((j==tensor(i)%elems(lii)%j).AND.(k==tensor(i)&
                 &%elems(lii)%k)).OR.((j==tensor(i)%elems(lii)%k).AND.(k==tensor(i)%elems(lii)%j))) THEN
               ! Found another entry with the same i,j,k: merge both into
               ! the one listed first (of those two). 
               tensor(i)%elems(lii)%v=tensor(i)%elems(lii)%v+tensor(i)%elems(li)%v
               IF (j>k) THEN
                  tensor(i)%elems(lii)%j=tensor(i)%elems(li)%k
                  tensor(i)%elems(lii)%k=tensor(i)%elems(li)%j
               ENDIF
               
               ! Shift the rest of the items one place down.
               DO liii=li+1,n
                  tensor(i)%elems(liii-1)%j=tensor(i)%elems(liii)%j
                  tensor(i)%elems(liii-1)%k=tensor(i)%elems(liii)%k
                  tensor(i)%elems(liii-1)%v=tensor(i)%elems(liii)%v
               END DO
               tensor(i)%nelems=tensor(i)%nelems-1
               ! Here we should stop because the li no longer points to the
               ! original i,j,k element
               EXIT
            ENDIF
         ENDDO
      ENDDO
      n=tensor(i)%nelems
      DO li=1,n
         ! Clear new "almost" zero entries and shift rest of the items one place down.
         ! Make sure not to skip any entries while shifting!
         DO WHILE (ABS(tensor(i)%elems(li)%v) < real_eps)
            DO liii=li+1,n
               tensor(i)%elems(liii-1)%j=tensor(i)%elems(liii)%j
               tensor(i)%elems(liii-1)%k=tensor(i)%elems(liii)%k
               tensor(i)%elems(liii-1)%v=tensor(i)%elems(liii)%v
            ENDDO
            tensor(i)%nelems=tensor(i)%nelems-1
            if (li > tensor(i)%nelems) THEN
               EXIT
            ENDIF
         ENDDO
      ENDDO

      n=tensor(i)%nelems
      DO li=1,n
         ! Upper triangularize
         j=tensor(i)%elems(li)%j
         k=tensor(i)%elems(li)%k
         IF (j>k) THEN
            tensor(i)%elems(li)%j=k
            tensor(i)%elems(li)%k=j
         ENDIF
      ENDDO


   ENDDO
 END SUBROUTINE simplify

    
  !> Subroutine to add element to a coolist.
  !> @param t destination tensor
  !> @param i tensor \f$i\f$ index
  !> @param j tensor \f$j\f$ index
  !> @param k tensor \f$k\f$ index
  !> @param v value to add
  SUBROUTINE add_elem(t,i,j,k,v)
    TYPE(coolist), DIMENSION(ndim), INTENT(INOUT) :: t
    INTEGER, INTENT(IN) :: i,j,k
    REAL(KIND=8), INTENT(IN) :: v
    INTEGER :: n
    IF (ABS(v) .ge. real_eps) THEN
       n=(t(i)%nelems)+1
       t(i)%elems(n)%j=j
       t(i)%elems(n)%k=k
       t(i)%elems(n)%v=v
       t(i)%nelems=n
    END IF
  END SUBROUTINE add_elem

  !> Subroutine to add element to a coolist and check for overflow.
  !> Once the t buffer tensor is full, add it to the destination buffer.
  !> @param t temporary buffer tensor for the destination tensor
  !> @param i tensor \f$i\f$ index
  !> @param j tensor \f$j\f$ index
  !> @param k tensor \f$k\f$ index
  !> @param v value to add
  !> @param dst destination tensor
  SUBROUTINE add_check(t,i,j,k,v,dst)
    TYPE(coolist), DIMENSION(ndim), INTENT(INOUT) :: t
    TYPE(coolist), DIMENSION(ndim), INTENT(INOUT) :: dst
    INTEGER, INTENT(IN) :: i,j,k
    REAL(KIND=8), INTENT(IN) :: v
    INTEGER :: n
    CALL add_elem(t,i,j,k,v)
    IF (t(i)%nelems==size(t(i)%elems)) THEN
       CALL add_to_tensor(t,dst)
       DO n=1,ndim
           t(n)%nelems=0
       ENDDO
    ENDIF
  END SUBROUTINE add_check

  !> Routine to add a rank-3 tensor to another one.
  !> @param src Tensor to add
  !> @param dst Destination tensor
  SUBROUTINE add_to_tensor(src,dst)
    TYPE(coolist), DIMENSION(ndim), INTENT(IN) :: src
    TYPE(coolist), DIMENSION(ndim), INTENT(INOUT) :: dst
    TYPE(coolist_elem), DIMENSION(:), ALLOCATABLE :: celems
    INTEGER :: i,j,n,AllocStat

    DO i=1,ndim
       IF (src(i)%nelems/=0) THEN
          IF (dst(i)%nelems==0) THEN
             IF (ALLOCATED(dst(i)%elems)) THEN
                DEALLOCATE(dst(i)%elems, STAT=AllocStat)
                IF (AllocStat /= 0) STOP "*** Deallocation problem ! ***"
             ENDIF
             ALLOCATE(dst(i)%elems(src(i)%nelems), STAT=AllocStat)
             IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
             n=0
          ELSE
             n=dst(i)%nelems
             ALLOCATE(celems(n), STAT=AllocStat)
             DO j=1,n
                celems(j)%j=dst(i)%elems(j)%j
                celems(j)%k=dst(i)%elems(j)%k
                celems(j)%v=dst(i)%elems(j)%v
             ENDDO
             DEALLOCATE(dst(i)%elems, STAT=AllocStat)
             IF (AllocStat /= 0) STOP "*** Deallocation problem ! ***"
             ALLOCATE(dst(i)%elems(src(i)%nelems+n), STAT=AllocStat)
             IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
             DO j=1,n
                dst(i)%elems(j)%j=celems(j)%j
                dst(i)%elems(j)%k=celems(j)%k
                dst(i)%elems(j)%v=celems(j)%v
             ENDDO
             DEALLOCATE(celems, STAT=AllocStat)
             IF (AllocStat /= 0) STOP "*** Deallocation problem ! ***"
          ENDIF
          DO j=1,src(i)%nelems
             dst(i)%elems(n+j)%j=src(i)%elems(j)%j
             dst(i)%elems(n+j)%k=src(i)%elems(j)%k
             dst(i)%elems(n+j)%v=src(i)%elems(j)%v
          ENDDO
          dst(i)%nelems=src(i)%nelems+n
       ENDIF
    ENDDO

  END SUBROUTINE add_to_tensor

  !> Routine to print a rank 3 tensor coolist.
  !> @param t coolist to print
  SUBROUTINE print_tensor(t,s)
    USE util, only: str
    TYPE(coolist), DIMENSION(ndim), INTENT(IN) :: t
    CHARACTER, INTENT(IN), OPTIONAL :: s
    CHARACTER :: r
    INTEGER :: i,n,j,k
    IF (PRESENT(s)) THEN
       r=s
    ELSE
       r="t"
    END IF
    DO i=1,ndim
       DO n=1,t(i)%nelems
          j=t(i)%elems(n)%j
          k=t(i)%elems(n)%k
          IF( ABS(t(i)%elems(n)%v) .GE. real_eps) THEN
             write(*,"(A,ES12.5)") r//"["//TRIM(str(i))//"]["//TRIM(str(j)) &
                  &//"]["//TRIM(str(k))//"] = ",t(i)%elems(n)%v
          END IF
       END DO
    END DO
  END SUBROUTINE print_tensor

  !> Load a rank-4 tensor coolist from a file definition
  !> @param s Destination filename
  !> @param t The coolist to write
  SUBROUTINE write_tensor_to_file(s,t)
    CHARACTER (LEN=*), INTENT(IN) :: s
    TYPE(coolist), DIMENSION(ndim), INTENT(IN) :: t
    INTEGER :: i,j,k,n
    OPEN(30,file=s)
    DO i=1,ndim
       WRITE(30,*) i,t(i)%nelems
       DO n=1,t(i)%nelems
          j=t(i)%elems(n)%j
          k=t(i)%elems(n)%k
          WRITE(30,*) i,j,k,t(i)%elems(n)%v
       END DO
    END DO
    CLOSE(30)
  END SUBROUTINE write_tensor_to_file

  !> Load a rank-4 tensor coolist from a file definition
  !> @param s Filename of the tensor definition file
  !> @param t The loaded coolist
  !> @remark The destination tensor have to be an empty tensor, i.e. with unallocated list of elements and nelems set to 0.
  SUBROUTINE load_tensor_from_file(s,t)
    CHARACTER (LEN=*), INTENT(IN) :: s
    TYPE(coolist), DIMENSION(ndim), INTENT(OUT) :: t
    INTEGER :: i,ir,j,k,n,AllocStat
    REAL(KIND=8) :: v
    OPEN(30,file=s,status='old')
    DO i=1,ndim
       READ(30,*) ir,n
       IF (n /= 0) THEN
          ALLOCATE(t(i)%elems(n), STAT=AllocStat)
          IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
          t(i)%nelems=n
       ENDIF
       DO n=1,t(i)%nelems
          READ(30,*) ir,j,k,v
          t(i)%elems(n)%j=j
          t(i)%elems(n)%k=k
          t(i)%elems(n)%v=v
       ENDDO
    END DO
    CLOSE(30)
  END SUBROUTINE load_tensor_from_file



  !> Routine to add a matrix to a rank-3 tensor.
  !> @param i   Add to tensor component i
  !> @param src Matrix to add
  !> @param dst Destination tensor
  SUBROUTINE add_matc_to_tensor(i,src,dst)
    INTEGER, INTENT(IN) :: i
    REAL(KIND=8), DIMENSION(ndim,ndim), INTENT(IN) :: src
    TYPE(coolist), DIMENSION(ndim), INTENT(INOUT) :: dst
    TYPE(coolist_elem), DIMENSION(:), ALLOCATABLE :: celems
    INTEGER :: j,k,r,n,nsrc,AllocStat

    nsrc=0
    DO j=1,ndim
       DO k=1,ndim
          IF (ABS(src(j,k))>real_eps) nsrc=nsrc+1
       END DO
    END DO

    IF (dst(i)%nelems==0) THEN
       IF (ALLOCATED(dst(i)%elems)) THEN
          DEALLOCATE(dst(i)%elems, STAT=AllocStat)
          IF (AllocStat /= 0) STOP "*** Deallocation problem ! ***"
       ENDIF
       ALLOCATE(dst(i)%elems(nsrc), STAT=AllocStat)
       IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
       n=0
    ELSE
       n=dst(i)%nelems
       ALLOCATE(celems(n), STAT=AllocStat)
       DO j=1,n
          celems(j)%j=dst(i)%elems(j)%j
          celems(j)%k=dst(i)%elems(j)%k
          celems(j)%v=dst(i)%elems(j)%v
       ENDDO
       DEALLOCATE(dst(i)%elems, STAT=AllocStat)
       IF (AllocStat /= 0) STOP "*** Deallocation problem ! ***"
       ALLOCATE(dst(i)%elems(nsrc+n), STAT=AllocStat)
       IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
       DO j=1,n
          dst(i)%elems(j)%j=celems(j)%j
          dst(i)%elems(j)%k=celems(j)%k
          dst(i)%elems(j)%v=celems(j)%v
       ENDDO
       DEALLOCATE(celems, STAT=AllocStat)
       IF (AllocStat /= 0) STOP "*** Deallocation problem ! ***"
    ENDIF
    r=0
    DO j=1,ndim
       DO k=1,ndim
          IF (ABS(src(j,k))>real_eps) THEN
             r=r+1
             dst(i)%elems(n+r)%j=j
             dst(i)%elems(n+r)%k=k
             dst(i)%elems(n+r)%v=src(j,k)
          ENDIF
       ENDDO
    END DO
    dst(i)%nelems=nsrc+n
    
  END SUBROUTINE add_matc_to_tensor


  !> Routine to add a matrix to a rank-4 tensor.
  !> @param i   Add to tensor component i,j
  !> @param j   Add to tensor component i,j
  !> @param src Matrix to add
  !> @param dst Destination tensor
  SUBROUTINE add_matc_to_tensor4(i,j,src,dst)
    INTEGER, INTENT(IN) :: i,j
    REAL(KIND=8), DIMENSION(ndim,ndim), INTENT(IN) :: src
    TYPE(coolist4), DIMENSION(ndim), INTENT(INOUT) :: dst
    TYPE(coolist_elem4), DIMENSION(:), ALLOCATABLE :: celems
    INTEGER :: k,l,r,n,nsrc,AllocStat

    nsrc=0
    DO k=1,ndim
       DO l=1,ndim
          IF (ABS(src(k,l))>real_eps) nsrc=nsrc+1
       END DO
    END DO

    IF (dst(i)%nelems==0) THEN
       IF (ALLOCATED(dst(i)%elems)) THEN
          DEALLOCATE(dst(i)%elems, STAT=AllocStat)
          IF (AllocStat /= 0) STOP "*** Deallocation problem ! ***"
       ENDIF
       ALLOCATE(dst(i)%elems(nsrc), STAT=AllocStat)
       IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
       n=0
    ELSE
       n=dst(i)%nelems
       ALLOCATE(celems(n), STAT=AllocStat)
       DO k=1,n
          celems(k)%j=dst(i)%elems(k)%j
          celems(k)%k=dst(i)%elems(k)%k
          celems(k)%l=dst(i)%elems(k)%l
          celems(k)%v=dst(i)%elems(k)%v
       ENDDO
       DEALLOCATE(dst(i)%elems, STAT=AllocStat)
       IF (AllocStat /= 0) STOP "*** Deallocation problem ! ***"
       ALLOCATE(dst(i)%elems(nsrc+n), STAT=AllocStat)
       IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
       DO k=1,n
          dst(i)%elems(k)%j=celems(k)%j
          dst(i)%elems(k)%k=celems(k)%k
          dst(i)%elems(k)%l=celems(k)%l
          dst(i)%elems(k)%v=celems(k)%v
       ENDDO
       DEALLOCATE(celems, STAT=AllocStat)
       IF (AllocStat /= 0) STOP "*** Deallocation problem ! ***"
    ENDIF
    r=0
    DO k=1,ndim
       DO l=1,ndim
          IF (ABS(src(k,l))>real_eps) THEN
             r=r+1
             dst(i)%elems(n+r)%j=j
             dst(i)%elems(n+r)%k=k
             dst(i)%elems(n+r)%l=l
             dst(i)%elems(n+r)%v=src(k,l)
          ENDIF
       ENDDO
    END DO
    dst(i)%nelems=nsrc+n
    
  END SUBROUTINE add_matc_to_tensor4


  !> Routine to add a vector to a rank-3 tensor.
  !> @param j,k   Add to tensor component j and k
  !> @param src Vector to add
  !> @param dst Destination tensor
  SUBROUTINE add_vec_jk_to_tensor(j,k,src,dst)
    INTEGER, INTENT(IN) :: j,k
    REAL(KIND=8), DIMENSION(ndim), INTENT(IN) :: src
    TYPE(coolist), DIMENSION(ndim), INTENT(INOUT) :: dst
    TYPE(coolist_elem), DIMENSION(:), ALLOCATABLE :: celems
    INTEGER :: i,l,r,n,nsrc,AllocStat

    DO i=1,ndim
       nsrc=0
       IF (ABS(src(i))>real_eps) nsrc=1
       IF (dst(i)%nelems==0) THEN
          IF (ALLOCATED(dst(i)%elems)) THEN
             DEALLOCATE(dst(i)%elems, STAT=AllocStat)
             IF (AllocStat /= 0) STOP "*** Deallocation problem ! ***"
          ENDIF
          ALLOCATE(dst(i)%elems(nsrc), STAT=AllocStat)
          IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
          n=0
       ELSE
          n=dst(i)%nelems
          ALLOCATE(celems(n), STAT=AllocStat)
          DO l=1,n
             celems(l)%j=dst(i)%elems(l)%j
             celems(l)%k=dst(i)%elems(l)%k
             celems(l)%v=dst(i)%elems(l)%v
          ENDDO
          DEALLOCATE(dst(i)%elems, STAT=AllocStat)
          IF (AllocStat /= 0) STOP "*** Deallocation problem ! ***"
          ALLOCATE(dst(i)%elems(nsrc+n), STAT=AllocStat)
          IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
          DO l=1,n
             dst(i)%elems(l)%j=celems(l)%j
             dst(i)%elems(l)%k=celems(l)%k
             dst(i)%elems(l)%v=celems(l)%v
          ENDDO
          DEALLOCATE(celems, STAT=AllocStat)
          IF (AllocStat /= 0) STOP "*** Deallocation problem ! ***"
       ENDIF
       r=0
       IF (ABS(src(i))>real_eps) THEN
          r=r+1
          dst(i)%elems(n+r)%j=j
          dst(i)%elems(n+r)%k=k
          dst(i)%elems(n+r)%v=src(i)
       ENDIF
       dst(i)%nelems=nsrc+n
    END DO
    
    
  END SUBROUTINE add_vec_jk_to_tensor

  !> Routine to add a vector to a rank-4 tensor plus permutation.
  !> @param i,k,l   Add to tensor component i,k and l
  !> @param src Vector to add
  !> @param dst Destination tensor
  SUBROUTINE add_vec_ikl_to_tensor4_perm(i,k,l,src,dst)
    INTEGER, INTENT(IN) :: i,k,l
    REAL(KIND=8), DIMENSION(ndim), INTENT(IN) :: src
    TYPE(coolist4), DIMENSION(ndim), INTENT(INOUT) :: dst
    TYPE(coolist_elem4), DIMENSION(:), ALLOCATABLE :: celems
    INTEGER :: j,ne,r,n,nsrc,AllocStat

    nsrc=0
    DO j=1,ndim
       IF (ABS(src(j))>real_eps) nsrc=nsrc+1
    ENDDO
    nsrc=nsrc*3
    IF (dst(i)%nelems==0) THEN
       IF (ALLOCATED(dst(i)%elems)) THEN
          DEALLOCATE(dst(i)%elems, STAT=AllocStat)
          IF (AllocStat /= 0) STOP "*** Deallocation problem ! ***"
       ENDIF
       ALLOCATE(dst(i)%elems(nsrc), STAT=AllocStat)
       IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
       n=0
    ELSE
       n=dst(i)%nelems
       ALLOCATE(celems(n), STAT=AllocStat)
       DO ne=1,n
          celems(ne)%j=dst(i)%elems(ne)%j
          celems(ne)%k=dst(i)%elems(ne)%k
          celems(ne)%l=dst(i)%elems(ne)%l
          celems(ne)%v=dst(i)%elems(ne)%v
       ENDDO
       DEALLOCATE(dst(i)%elems, STAT=AllocStat)
       IF (AllocStat /= 0) STOP "*** Deallocation problem ! ***"
       ALLOCATE(dst(i)%elems(nsrc+n), STAT=AllocStat)
       IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
       DO ne=1,n
          dst(i)%elems(ne)%j=celems(ne)%j
          dst(i)%elems(ne)%k=celems(ne)%k
          dst(i)%elems(ne)%l=celems(ne)%l
          dst(i)%elems(ne)%v=celems(ne)%v
       ENDDO
       DEALLOCATE(celems, STAT=AllocStat)
       IF (AllocStat /= 0) STOP "*** Deallocation problem ! ***"
    ENDIF
    r=0
    DO j=1,ndim
       IF (ABS(src(j))>real_eps) THEN
          r=r+1
          dst(i)%elems(n+r)%j=j
          dst(i)%elems(n+r)%k=k
          dst(i)%elems(n+r)%l=l
          dst(i)%elems(n+r)%v=src(j)
          r=r+1
          dst(i)%elems(n+r)%j=k
          dst(i)%elems(n+r)%k=l
          dst(i)%elems(n+r)%l=j
          dst(i)%elems(n+r)%v=src(j)
          r=r+1
          dst(i)%elems(n+r)%j=l
          dst(i)%elems(n+r)%k=j
          dst(i)%elems(n+r)%l=k
          dst(i)%elems(n+r)%v=src(j)
       ENDIF
    ENDDO
    dst(i)%nelems=nsrc+n
  END SUBROUTINE add_vec_ikl_to_tensor4_perm

  !> Routine to add a vector to a rank-4 tensor.
  !> @param i,k,l   Add to tensor component i,k and l
  !> @param src Vector to add
  !> @param dst Destination tensor
  SUBROUTINE add_vec_ikl_to_tensor4(i,k,l,src,dst)
    INTEGER, INTENT(IN) :: i,k,l
    REAL(KIND=8), DIMENSION(ndim), INTENT(IN) :: src
    TYPE(coolist4), DIMENSION(ndim), INTENT(INOUT) :: dst
    TYPE(coolist_elem4), DIMENSION(:), ALLOCATABLE :: celems
    INTEGER :: j,ne,r,n,nsrc,AllocStat

    nsrc=0
    DO j=1,ndim
       IF (ABS(src(j))>real_eps) nsrc=nsrc+1
    ENDDO

    IF (dst(i)%nelems==0) THEN
       IF (ALLOCATED(dst(i)%elems)) THEN
          DEALLOCATE(dst(i)%elems, STAT=AllocStat)
          IF (AllocStat /= 0) STOP "*** Deallocation problem ! ***"
       ENDIF
       ALLOCATE(dst(i)%elems(nsrc), STAT=AllocStat)
       IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
       n=0
    ELSE
       n=dst(i)%nelems
       ALLOCATE(celems(n), STAT=AllocStat)
       DO ne=1,n
          celems(ne)%j=dst(i)%elems(ne)%j
          celems(ne)%k=dst(i)%elems(ne)%k
          celems(ne)%l=dst(i)%elems(ne)%l
          celems(ne)%v=dst(i)%elems(ne)%v
       ENDDO
       DEALLOCATE(dst(i)%elems, STAT=AllocStat)
       IF (AllocStat /= 0) STOP "*** Deallocation problem ! ***"
       ALLOCATE(dst(i)%elems(nsrc+n), STAT=AllocStat)
       IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
       DO ne=1,n
          dst(i)%elems(ne)%j=celems(ne)%j
          dst(i)%elems(ne)%k=celems(ne)%k
          dst(i)%elems(ne)%l=celems(ne)%l
          dst(i)%elems(ne)%v=celems(ne)%v
       ENDDO
       DEALLOCATE(celems, STAT=AllocStat)
       IF (AllocStat /= 0) STOP "*** Deallocation problem ! ***"
    ENDIF
    r=0
    DO j=1,ndim
       IF (ABS(src(j))>real_eps) THEN
          r=r+1
          dst(i)%elems(n+r)%j=j
          dst(i)%elems(n+r)%k=k
          dst(i)%elems(n+r)%l=l
          dst(i)%elems(n+r)%v=src(j)
       ENDIF
    ENDDO
    dst(i)%nelems=nsrc+n
  END SUBROUTINE add_vec_ikl_to_tensor4

  !> Routine to add a vector to a rank-4 tensor.
  !> @param i,j,k   Add to tensor component i,j and k
  !> @param src Vector to add
  !> @param dst Destination tensor
  SUBROUTINE add_vec_ijk_to_tensor4(i,j,k,src,dst)
    INTEGER, INTENT(IN) :: i,j,k
    REAL(KIND=8), DIMENSION(ndim), INTENT(IN) :: src
    TYPE(coolist4), DIMENSION(ndim), INTENT(INOUT) :: dst
    TYPE(coolist_elem4), DIMENSION(:), ALLOCATABLE :: celems
    INTEGER :: l,ne,r,n,nsrc,AllocStat

    nsrc=0
    DO l=1,ndim
       IF (ABS(src(l))>real_eps) nsrc=nsrc+1
    ENDDO

    IF (dst(i)%nelems==0) THEN
       IF (ALLOCATED(dst(i)%elems)) THEN
          DEALLOCATE(dst(i)%elems, STAT=AllocStat)
          IF (AllocStat /= 0) STOP "*** Deallocation problem ! ***"
       ENDIF
       ALLOCATE(dst(i)%elems(nsrc), STAT=AllocStat)
       IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
       n=0
    ELSE
       n=dst(i)%nelems
       ALLOCATE(celems(n), STAT=AllocStat)
       DO ne=1,n
          celems(ne)%j=dst(i)%elems(ne)%j
          celems(ne)%k=dst(i)%elems(ne)%k
          celems(ne)%l=dst(i)%elems(ne)%l
          celems(ne)%v=dst(i)%elems(ne)%v
       ENDDO
       DEALLOCATE(dst(i)%elems, STAT=AllocStat)
       IF (AllocStat /= 0) STOP "*** Deallocation problem ! ***"
       ALLOCATE(dst(i)%elems(nsrc+n), STAT=AllocStat)
       IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
       DO ne=1,n
          dst(i)%elems(ne)%j=celems(ne)%j
          dst(i)%elems(ne)%k=celems(ne)%k
          dst(i)%elems(ne)%l=celems(ne)%l
          dst(i)%elems(ne)%v=celems(ne)%v
       ENDDO
       DEALLOCATE(celems, STAT=AllocStat)
       IF (AllocStat /= 0) STOP "*** Deallocation problem ! ***"
    ENDIF
    r=0
    DO l=1,ndim
       IF (ABS(src(l))>real_eps) THEN
          r=r+1
          dst(i)%elems(n+r)%j=j
          dst(i)%elems(n+r)%k=k
          dst(i)%elems(n+r)%l=l
          dst(i)%elems(n+r)%v=src(l)
       ENDIF
    ENDDO
    dst(i)%nelems=nsrc+n
  END SUBROUTINE add_vec_ijk_to_tensor4


 

  !> Routine to convert a rank-3 tensor from matrix to coolist representation.
  !> @param src Source matrix
  !> @param dst Destination coolist
  !> @remark The destination coolist have to be an empty one, i.e. with unallocated list of elements and nelems set to 0.
  SUBROUTINE tensor_to_coo(src,dst)
    REAL(KIND=8), DIMENSION(ndim,0:ndim,0:ndim), INTENT(IN) :: src
    TYPE(coolist), DIMENSION(ndim), INTENT(OUT) :: dst
    INTEGER :: i,j,k,n,AllocStat
    
    DO i=1,ndim
       n=0
       DO j=0,ndim
          DO k=0,ndim
             IF (ABS(src(i,j,k))>real_eps) n=n+1
          ENDDO
       ENDDO
       IF (n/=0) THEN
          IF (dst(i)%nelems/=0) STOP "*** tensor_to_coo : Destination coolist not empty ! ***"
          ALLOCATE(dst(i)%elems(n), STAT=AllocStat)
          IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
          n=0
          DO j=0,ndim
             DO k=0,ndim
                IF (ABS(src(i,j,k))>real_eps) THEN
                   n=n+1
                   dst(i)%elems(n)%j=j
                   dst(i)%elems(n)%k=k
                   dst(i)%elems(n)%v=src(i,j,k)
                ENDIF
             ENDDO
          ENDDO
       ENDIF
       dst(i)%nelems=n
    ENDDO
  END SUBROUTINE tensor_to_coo

  !> Routine to convert a rank-4 tensor from matrix to coolist representation.
  !> @param src Source matrix
  !> @param dst Destination coolist
  !> @remark The destination coolist have to be an empty one, i.e. with unallocated list of elements and nelems set to 0.
  SUBROUTINE tensor4_to_coo4(src,dst)
    REAL(KIND=8), DIMENSION(ndim,0:ndim,0:ndim,0:ndim), INTENT(IN) :: src
    TYPE(coolist4), DIMENSION(ndim), INTENT(OUT) :: dst
    INTEGER :: i,j,k,l,n,AllocStat
    
    DO i=1,ndim
       n=0
       DO j=0,ndim
          DO k=0,ndim
             DO l=0,ndim
                IF (ABS(src(i,j,k,l))>real_eps) n=n+1
             ENDDO
          ENDDO
       ENDDO
       IF (n/=0) THEN
          IF (dst(i)%nelems/=0) STOP "*** tensor_to_coo : Destination coolist not empty ! ***"
          ALLOCATE(dst(i)%elems(n), STAT=AllocStat)
          IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
          n=0
          DO j=0,ndim
             DO k=0,ndim
                DO l=0,ndim
                   IF (ABS(src(i,j,k,l))>real_eps) THEN
                      n=n+1
                      dst(i)%elems(n)%j=j
                      dst(i)%elems(n)%k=k
                      dst(i)%elems(n)%l=l
                      dst(i)%elems(n)%v=src(i,j,k,l)
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
       ENDIF
       dst(i)%nelems=n
    ENDDO
  END SUBROUTINE tensor4_to_coo4

  !> Routine to print a rank-4 tensor coolist.
  !> @param t coolist to print
  SUBROUTINE print_tensor4(t)
    USE util, only: str
    TYPE(coolist4), DIMENSION(ndim), INTENT(IN) :: t
    INTEGER :: i,n,j,k,l
    DO i=1,ndim
       DO n=1,t(i)%nelems
          j=t(i)%elems(n)%j
          k=t(i)%elems(n)%k
          l=t(i)%elems(n)%l
          IF( ABS(t(i)%elems(n)%v) .GE. real_eps) THEN
             write(*,"(A,ES12.5)") "tensor["//TRIM(str(i))//"]["//TRIM(str(j)) &
                  &//"]["//TRIM(str(k))//"]["//TRIM(str(l))//"] = ",t(i)%elems(n)%v
          END IF
       END DO
    END DO
  END SUBROUTINE print_tensor4


  !> Sparse multiplication of a rank-3 tensor coolist with a vector: 
  !> \f${\displaystyle \sum_{k=0}^{ndim}} \mathcal{T}_{i,j,k} \, b_k\f$.
  !> Its output is a matrix.
  !> @param coolist_ijk a coolist (sparse tensor) of which index k will be contracted.
  !> @param arr_k the vector to be contracted with index k of coolist_ijk
  !> @param res matrix (buffer) to store the result of the contraction
  !> @remark Note that it is NOT safe to pass `arr_k` as a result buffer, 
  !> as this operation does multiple passes.
  SUBROUTINE sparse_mul3_mat(coolist_ijk, arr_k, res)
    TYPE(coolist), DIMENSION(ndim), INTENT(IN):: coolist_ijk
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(IN)  :: arr_k
    REAL(KIND=8), DIMENSION(ndim,ndim), INTENT(OUT) :: res
    INTEGER :: i,j,k,n
    res=0.D0
    DO i=1,ndim
       DO n=1,coolist_ijk(i)%nelems
         j=coolist_ijk(i)%elems(n)%j
         IF (j /= 0) THEN
            k=coolist_ijk(i)%elems(n)%k
            res(i,j) = res(i,j) + coolist_ijk(i)%elems(n)%v * arr_k(k)
         ENDIF
      END DO
   END DO
 END SUBROUTINE sparse_mul3_mat


  !> Sparse multiplication of a rank-4 tensor coolist with three vectors:  \f${\displaystyle \sum_{j,k,l=0}^{ndim}} \mathcal{T}_{i,j,k,l} \, a_j \,b_k \, c_l \f$.
  !> @param coolist_ijkl a coolist (sparse tensor) of which index j, k and l will be contracted.
  !> @param arr_j the vector to be contracted with index j of coolist_ijkl
  !> @param arr_k the vector to be contracted with index k of coolist_ijkl
  !> @param arr_l the vector to be contracted with index l of coolist_ijkl
  !> @param res vector (buffer) to store the result of the contraction
  !> @remark Note that it is NOT safe to pass `arr_j`/`arr_k`/`arr_l` as a result buffer, 
  !> as this operation does multiple passes.
  SUBROUTINE sparse_mul4(coolist_ijkl, arr_j, arr_k, arr_l, res)
    TYPE(coolist4), DIMENSION(ndim), INTENT(IN):: coolist_ijkl
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(IN)  :: arr_j, arr_k, arr_l
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(OUT) :: res
    INTEGER :: i,j,k,n,l
    res=0.D0
    DO i=1,ndim
       DO n=1,coolist_ijkl(i)%nelems
         j=coolist_ijkl(i)%elems(n)%j
         k=coolist_ijkl(i)%elems(n)%k
         l=coolist_ijkl(i)%elems(n)%l
         res(i) = res(i) + coolist_ijkl(i)%elems(n)%v * arr_j(j)*arr_k(k)*arr_l(l)
      END DO
   END DO
 END SUBROUTINE sparse_mul4

  !> Sparse multiplication of a tensor with two vectors:  \f${\displaystyle \sum_{k,l=0}^{ndim}} \mathcal{T}_{i,j,k,l}  \,b_k \, c_l \f$.
  !> @param coolist_ijkl a coordinate list (sparse tensor) of which index
  !>  3 and 4 will be contracted.
  !> @param arr_k the vector to be contracted with index 3 of coolist_ijkl
  !> @param arr_l the vector to be contracted with index 4 of coolist_ijkl
  !> @param res matrix (buffer) to store the result of the contraction
  !> @remark Note that it is NOT safe to pass `arr_k`/`arr_l` as a result buffer, 
  !> as this operation does multiple passes.
  SUBROUTINE sparse_mul4_mat(coolist_ijkl, arr_k, arr_l, res)
    TYPE(coolist4), DIMENSION(ndim), INTENT(IN):: coolist_ijkl
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(IN)  :: arr_k, arr_l
    REAL(KIND=8), DIMENSION(ndim,ndim), INTENT(OUT) :: res
    INTEGER :: i,j,k,n,l
    res=0.D0
    DO i=1,ndim
       DO n=1,coolist_ijkl(i)%nelems
         j=coolist_ijkl(i)%elems(n)%j
         IF (j /= 0) THEN
            k=coolist_ijkl(i)%elems(n)%k
            l=coolist_ijkl(i)%elems(n)%l
            res(i,j) = res(i,j) + coolist_ijkl(i)%elems(n)%v * arr_k(k) * arr_l(l)
         ENDIF
      END DO
   END DO
 END SUBROUTINE sparse_mul4_mat


  !> Sparse multiplication of a 3d sparse tensor with a vectors:  \f${\displaystyle \sum_{j=0}^{ndim}} \mathcal{T}_{i,j,k} \, a_j \f$.
  !> @param coolist_ijk a coordinate list (sparse tensor) of which index
  !> 2 will be contracted.
  !> @param arr_j the vector to be contracted with index 2 of coolist_ijk
  !> @param res vector (buffer) to store the result of the contraction
  !> @remark Note that it is NOT safe to pass `arr_j` as a result buffer, 
  !> as this operation does multiple passes.
  SUBROUTINE sparse_mul2_j(coolist_ijk, arr_j, res)
    TYPE(coolist), DIMENSION(ndim), INTENT(IN):: coolist_ijk
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(IN)  :: arr_j
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(OUT) :: res
    INTEGER :: i,j,n
    res=0.D0
    DO i=1,ndim
       DO n=1,coolist_ijk(i)%nelems
         j=coolist_ijk(i)%elems(n)%j
         res(i) = res(i) + coolist_ijk(i)%elems(n)%v * arr_j(j)
      END DO
   END DO
 END SUBROUTINE sparse_mul2_j

  !> Sparse multiplication of a rank-3 sparse tensor coolist with a vector:  \f${\displaystyle \sum_{k=0}^{ndim}} \mathcal{T}_{i,j,k} \, a_k \f$.
  !> @param coolist_ijk a coordinate list (sparse tensor) of which index
  !> k will be contracted.
  !> @param arr_k the vector to be contracted with index k of coolist_ijk
  !> @param res vector (buffer) to store the result of the contraction
  !> @remark Note that it is NOT safe to pass `arr_k` as a result buffer, 
  !> as this operation does multiple passes.
  SUBROUTINE sparse_mul2_k(coolist_ijk, arr_k, res)
    TYPE(coolist), DIMENSION(ndim), INTENT(IN):: coolist_ijk
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(IN)  :: arr_k
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(OUT) :: res
    INTEGER :: i,k,n
    res=0.D0
    DO i=1,ndim
       DO n=1,coolist_ijk(i)%nelems
         k=coolist_ijk(i)%elems(n)%k
         res(i) = res(i) + coolist_ijk(i)%elems(n)%v * arr_k(k)
      END DO
   END DO
 END SUBROUTINE sparse_mul2_k

 
 !> Routine to convert a rank-3 tensor coolist component into a matrix with i and k indices.
 !> @param src Source tensor
 !> @param dst Destination matrix
 SUBROUTINE coo_to_mat_ik(src,dst)
   TYPE(coolist), DIMENSION(ndim), INTENT(IN) :: src
   REAL(KIND=8), DIMENSION(ndim,ndim), INTENT(OUT) :: dst
   INTEGER :: i,n
   
   dst=0.D0
   DO i=1,ndim
      DO n=1,src(i)%nelems
         dst(i,src(i)%elems(n)%k)=src(i)%elems(n)%v
      ENDDO
   ENDDO
 END SUBROUTINE coo_to_mat_ik

 !> Routine to convert a rank-3 tensor coolist component into a matrix with i and j indices.
 !> @param src Source tensor
 !> @param dst Destination matrix
 SUBROUTINE coo_to_mat_ij(src,dst)
   TYPE(coolist), DIMENSION(ndim), INTENT(IN) :: src
   REAL(KIND=8), DIMENSION(ndim,ndim), INTENT(OUT) :: dst
   INTEGER :: i,n
   
   dst=0.D0
   DO i=1,ndim
      DO n=1,src(i)%nelems
         dst(i,src(i)%elems(n)%j)=src(i)%elems(n)%v
      ENDDO
   ENDDO
 END SUBROUTINE coo_to_mat_ij

 ! SUBROUTINE tensor_perm_ij(t)
 !   TYPE(coolist), DIMENSION(ndim), INTENT(INOUT) :: t
 !   INTEGER :: i,j,k,n
   
 !   DO i=1,ndim
 !      DO n=1,t(i)%nelems
 !         j=t(i)%elems(n)%j
 !         k=t(i)%elems(n)%k
         
 !         t(i)%elems(n)%v
 !      ENDDO
 !   ENDDO
 ! END SUBROUTINE tensor_perm_ij

!!! not so cool

 !> Routine to convert a rank-3 tensor coolist component into a matrix.
 !> @param i   Component to convert
 !> @param src Source tensor
 !> @param dst Destination matrix
 SUBROUTINE coo_to_mat_i(i,src,dst)
   INTEGER, INTENT(IN) :: i
   TYPE(coolist), DIMENSION(ndim), INTENT(IN) :: src
   REAL(KIND=8), DIMENSION(ndim,ndim), INTENT(OUT) :: dst
   INTEGER :: n
   
   dst=0.D0
   DO n=1,src(i)%nelems
      dst(src(i)%elems(n)%j,src(i)%elems(n)%k)=src(i)%elems(n)%v
   ENDDO
 END SUBROUTINE coo_to_mat_i

 !> Routine to convert a rank-3 tensor coolist component into a vector.
 !> @param j   Component j,k to convert
 !> @param k   Component j,k to convert
 !> @param src Source tensor
 !> @param dst Destination vector
 SUBROUTINE coo_to_vec_jk(j,k,src,dst)
   INTEGER, INTENT(IN) :: j,k
   TYPE(coolist), DIMENSION(ndim), INTENT(IN) :: src
   REAL(KIND=8), DIMENSION(ndim), INTENT(OUT) :: dst
   INTEGER :: i,n
   
   dst=0.D0
   DO i=1,ndim
      DO n=1,src(i)%nelems
         IF ((src(i)%elems(n)%j==j).and.(src(i)%elems(n)%k==k)) dst(i)=src(i)%elems(n)%v
      END DO
   ENDDO
 END SUBROUTINE coo_to_vec_jk


 !> Routine to convert a rank-3 tensor coolist component into a matrix.
 !> @param j   Component to convert
 !> @param src Source tensor
 !> @param dst Destination matrix
 SUBROUTINE coo_to_mat_j(j,src,dst)
   INTEGER, INTENT(IN) :: j
   TYPE(coolist), DIMENSION(ndim), INTENT(IN) :: src
   REAL(KIND=8), DIMENSION(ndim,ndim), INTENT(OUT) :: dst
   INTEGER :: i,n
   
   dst=0.D0
   DO i=1,ndim
      DO n=1,src(i)%nelems
         IF (src(i)%elems(n)%j==j) dst(i,src(i)%elems(n)%k)=src(i)%elems(n)%v
      ENDDO
   END DO
 END SUBROUTINE coo_to_mat_j


  !> Sparse multiplication of a rank-4 tensor coolist with a matrix :  \f${\displaystyle \sum_{j,l=0}^{ndim}} \mathcal{T}_{i,j,k,l} \, m_{j,l} \f$.
  !> @param coolist_ijkl a coolist (sparse tensor) of which index j and l will be contracted.
  !> @param mat_jl the matrix to be contracted with indices j and l of coolist_ijkl
  !> @param res matrix (buffer) to store the result of the contraction
  !> @remark Note that it is NOT safe to pass `mat_jl` as a result buffer, 
  !> as this operation does multiple passes.
 SUBROUTINE sparse_mul4_with_mat_jl(coolist_ijkl,mat_jl,res)
   TYPE(coolist4), DIMENSION(ndim), INTENT(IN):: coolist_ijkl
   REAL(KIND=8), DIMENSION(ndim,ndim), INTENT(IN)  :: mat_jl
   REAL(KIND=8), DIMENSION(ndim,ndim), INTENT(OUT) :: res
   INTEGER i,j,k,l,n

   res=0.D0
   DO i=1,ndim
      DO n=1,coolist_ijkl(i)%nelems
         j=coolist_ijkl(i)%elems(n)%j
         k=coolist_ijkl(i)%elems(n)%k
         l=coolist_ijkl(i)%elems(n)%l

         res(i,k) = res(i,k) + coolist_ijkl(i)%elems(n)%v * mat_jl(j,l)
      ENDDO
   END DO

 END SUBROUTINE sparse_mul4_with_mat_jl

 !> Sparse multiplication of a rank-4 tensor coolist with a matrix :  \f${\displaystyle \sum_{j,l=0}^{ndim}} \mathcal{T}_{i,j,k,l} \, m_{k,l} \f$.
 !> @param coolist_ijkl a coolist (sparse tensor) of which index k and l will be contracted.
 !> @param mat_kl the matrix to be contracted with indices k and l of coolist_ijkl
 !> @param res matrix (buffer) to store the result of the contraction
 !> @remark Note that it is NOT safe to pass `mat_kl` as a result buffer, 
 !> as this operation does multiple passes.
 SUBROUTINE sparse_mul4_with_mat_kl(coolist_ijkl,mat_kl,res)
   TYPE(coolist4), DIMENSION(ndim), INTENT(IN):: coolist_ijkl
   REAL(KIND=8), DIMENSION(ndim,ndim), INTENT(IN)  :: mat_kl
   REAL(KIND=8), DIMENSION(ndim,ndim), INTENT(OUT) :: res
   INTEGER i,j,k,l,n

   res=0.D0
   DO i=1,ndim
      DO n=1,coolist_ijkl(i)%nelems
         j=coolist_ijkl(i)%elems(n)%j
         k=coolist_ijkl(i)%elems(n)%k
         l=coolist_ijkl(i)%elems(n)%l

         res(i,j) = res(i,j) + coolist_ijkl(i)%elems(n)%v * mat_kl(k,l)
      ENDDO
   END DO

 END SUBROUTINE sparse_mul4_with_mat_kl

  !> Sparse multiplication of a rank-3 tensor coolist with a matrix:  \f${\displaystyle \sum_{j,k=0}^{ndim}} \mathcal{T}_{i,j,k} \, m_{j,k}\f$.
  !> @param coolist_ijk a coolist (sparse tensor) of which index
  !> j and k will be contracted.
  !> @param mat_jk the matrix to be contracted with index j and k of coolist_ijk
  !> @param res vector (buffer) to store the result of the contraction
  !> @remark Note that it is NOT safe to pass `mat_jk` as a result buffer, 
  !> as this operation does multiple passes.
 SUBROUTINE sparse_mul3_with_mat(coolist_ijk,mat_jk,res)
   TYPE(coolist), DIMENSION(ndim), INTENT(IN):: coolist_ijk
   REAL(KIND=8), DIMENSION(ndim,ndim), INTENT(IN)  :: mat_jk
   REAL(KIND=8), DIMENSION(0:ndim), INTENT(OUT) :: res
   INTEGER i,j,k,n

   res=0.D0
   DO i=1,ndim
      DO n=1,coolist_ijk(i)%nelems
         j=coolist_ijk(i)%elems(n)%j
         k=coolist_ijk(i)%elems(n)%k

         res(i) = res(i) + coolist_ijk(i)%elems(n)%v * mat_jk(j,k)
      ENDDO
   END DO

 END SUBROUTINE sparse_mul3_with_mat


 !> Routine to convert a matrix to a rank-3 tensor.
 !> @param src Source matrix
 !> @param dst Destination tensor
 !> @remark The destination tensor have to be an empty tensor, i.e. with unallocated list of elements and nelems set to 0.
 !> @remark The j component will be set to 0.
 SUBROUTINE matc_to_coo(src,dst)
   REAL(KIND=8), DIMENSION(ndim,ndim), INTENT(IN) :: src
   TYPE(coolist), DIMENSION(ndim), INTENT(OUT) :: dst
   INTEGER :: i,j,n,AllocStat
   DO i=1,ndim
      n=0
      DO j=1,ndim
         IF (ABS(src(i,j))>real_eps) n=n+1
      ENDDO
      IF (n/=0) THEN
         IF (dst(i)%nelems/=0) STOP "*** mat_to_coo : Destination coolist not empty ! ***"
         ALLOCATE(dst(i)%elems(n), STAT=AllocStat)
         IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
         n=0
         DO j=1,ndim
            IF (ABS(src(i,j))>real_eps) THEN
               n=n+1
               dst(i)%elems(n)%j=0
               dst(i)%elems(n)%k=j
               dst(i)%elems(n)%v=src(i,j)
            ENDIF
         ENDDO
      ENDIF
      dst(i)%nelems=n
   ENDDO
 END SUBROUTINE matc_to_coo
  
 !> Routine to multiply a rank-3 tensor by a scalar
 !> @param s The scalar
 !> @param t The tensor
 SUBROUTINE scal_mul_coo(s,t)
   REAL(KIND=8), INTENT(IN) :: s 
   TYPE(coolist), DIMENSION(ndim), INTENT(INOUT) :: t
   INTEGER :: i,li,n
   DO i=1,ndim
      n=t(i)%nelems
      DO li=1,n
        t(i)%elems(li)%v=s*t(i)%elems(li)%v
      ENDDO
   ENDDO
 END SUBROUTINE scal_mul_coo

  !> Test if a rank-3 tensor coolist is empty
  !> @param t rank-3 tensor coolist to be tested 
  FUNCTION tensor_empty(t)
    TYPE(coolist), DIMENSION(ndim), INTENT(IN) :: t
    LOGICAL :: tensor_empty
    INTEGER :: i
    tensor_empty=.true.
    DO i=1,ndim
       IF (t(i)%nelems /= 0) THEN
          tensor_empty=.false.
          RETURN
       ENDIF
    END DO
    RETURN
  END FUNCTION tensor_empty

  !> Test if a rank-4 tensor coolist is empty
  !> @param t rank-4 tensor coolist to be tested 
  FUNCTION tensor4_empty(t)
    TYPE(coolist4), DIMENSION(ndim), INTENT(IN) :: t
    LOGICAL :: tensor4_empty
    INTEGER :: i
    tensor4_empty=.true.
    DO i=1,ndim
       IF (t(i)%nelems /= 0) THEN
          tensor4_empty=.false.
          RETURN
       ENDIF
    END DO
    RETURN
  END FUNCTION tensor4_empty

  !> Load a rank-4 tensor coolist from a file definition
  !> @param s Filename of the tensor definition file
  !> @param t The loaded coolist
  !> @remark The destination tensor have to be an empty tensor, i.e. with unallocated list of elements and nelems set to 0.
  SUBROUTINE load_tensor4_from_file(s,t)
    CHARACTER (LEN=*), INTENT(IN) :: s
    TYPE(coolist4), DIMENSION(ndim), INTENT(OUT) :: t
    INTEGER :: i,ir,j,k,l,n,AllocStat
    REAL(KIND=8) :: v
    OPEN(30,file=s,status='old')
    DO i=1,ndim
       READ(30,*) ir,n
       IF (n /= 0) THEN
          ALLOCATE(t(i)%elems(n), STAT=AllocStat)
          IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
          t(i)%nelems=n
       ENDIF
       DO n=1,t(i)%nelems
          READ(30,*) ir,j,k,l,v
          t(i)%elems(n)%j=j
          t(i)%elems(n)%k=k
          t(i)%elems(n)%l=l
          t(i)%elems(n)%v=v
       ENDDO
    END DO
    CLOSE(30)
  END SUBROUTINE load_tensor4_from_file

  !> Load a rank-4 tensor coolist from a file definition
  !> @param s Destination filename
  !> @param t The coolist to write
  SUBROUTINE write_tensor4_to_file(s,t)
    CHARACTER (LEN=*), INTENT(IN) :: s
    TYPE(coolist4), DIMENSION(ndim), INTENT(IN) :: t
    INTEGER :: i,j,k,l,n
    OPEN(30,file=s)
    DO i=1,ndim
       WRITE(30,*) i,t(i)%nelems
       DO n=1,t(i)%nelems
          j=t(i)%elems(n)%j
          k=t(i)%elems(n)%k
          l=t(i)%elems(n)%l
          WRITE(30,*) i,j,k,l,t(i)%elems(n)%v
       END DO
    END DO
    CLOSE(30)
  END SUBROUTINE write_tensor4_to_file
    
END MODULE tensor

