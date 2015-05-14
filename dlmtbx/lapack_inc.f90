!!!
!!! lapack and blas interfaces

!DEC$ IF .NOT. DEFINED(_lapack_include)
!DEC$ DEFINE _lapack_include

!!! BLAS
! dcopy
INTERFACE 
   SUBROUTINE DCOPY(N,DX,INCX,DY,INCY)
     INTEGER(KIND=4) :: N
     REAL(KIND=8) :: DX(*)
     INTEGER(KIND=4) :: INCX
     REAL(KIND=8) :: DY(*)
     INTEGER(KIND=4) :: INCY
   END SUBROUTINE DCOPY
END INTERFACE

INTERFACE 
   SUBROUTINE DSCAL(N,A,X,INCX)
     INTEGER(KIND=4) :: N
     REAL(KIND=8) :: A
     REAL(KIND=8) :: X(*)
     INTEGER(KIND=4) :: INCX
   END SUBROUTINE DSCAL
END INTERFACE


INTERFACE 
   SUBROUTINE DSYMV(UPLO,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
     INTEGER(KIND=4) :: LDA
     CHARACTER(LEN=1) :: UPLO
     INTEGER(KIND=4) :: N
     REAL(KIND=8) :: ALPHA
     REAL(KIND=8) :: A(LDA,*)
     REAL(KIND=8) :: X(*)
     INTEGER(KIND=4) :: INCX
     REAL(KIND=8) :: BETA
     REAL(KIND=8) :: Y(*)
     INTEGER(KIND=4) :: INCY
   END SUBROUTINE DSYMV
END INTERFACE

INTERFACE 
   SUBROUTINE DTRMV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
     INTEGER(KIND=4) :: LDA
     CHARACTER(LEN=1) :: UPLO
     CHARACTER(LEN=1) :: TRANS
     CHARACTER(LEN=1) :: DIAG
     INTEGER(KIND=4) :: N
     REAL(KIND=8) :: A(LDA,*)
     REAL(KIND=8) :: X(*)
     INTEGER(KIND=4) :: INCX
   END SUBROUTINE DTRMV
END INTERFACE

INTERFACE 
   SUBROUTINE DGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
     INTEGER(KIND=4) :: LDA
     CHARACTER(LEN=1) :: TRANS
     INTEGER(KIND=4) :: M
     INTEGER(KIND=4) :: N
     REAL(KIND=8) :: ALPHA
     REAL(KIND=8) :: A(LDA,*)
     REAL(KIND=8) :: X(*)
     INTEGER(KIND=4) :: INCX
     REAL(KIND=8) :: BETA
     REAL(KIND=8) :: Y(*)
     INTEGER(KIND=4) :: INCY
   END SUBROUTINE DGEMV
END INTERFACE
 
interface
     subroutine dgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
       implicit none
       character(len=1), intent(in) :: transa, transb
       integer(kind=4), intent(in):: lda, ldb,ldc
       integer(kind=4), intent(in):: m,n,k
       real(kind=8), intent(in) :: alpha, beta
       real(kind=8), intent(in) :: a(lda,*), b(ldb,*)
       real(kind=8), intent(inout) :: c(ldc,*)
     end subroutine dgemm
end interface

interface
     subroutine dsymm(side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc)
       implicit none
       character(len=1), intent(in) :: side, uplo
       integer(kind=4), intent(in):: lda, ldb,ldc
       integer(kind=4), intent(in):: m,n
       real(kind=8), intent(in) :: alpha, beta
       real(kind=8), intent(in) :: a(lda,*), b(ldb,*)
       real(kind=8), intent(inout) :: c(ldc,*)
     end subroutine dsymm
end interface


 ! rank one update of symmetric matrix
INTERFACE 
   SUBROUTINE DSYR(UPLO,N,ALPHA,X,INCX,A,LDA)
     INTEGER(KIND=4) :: LDA
     CHARACTER(LEN=1) :: UPLO
     INTEGER(KIND=4) :: N
     REAL(KIND=8) :: ALPHA
     REAL(KIND=8) :: X(*)
     INTEGER(KIND=4) :: INCX
     REAL(KIND=8) :: A(LDA,*)
   END SUBROUTINE DSYR
END INTERFACE

!!! LAPACK
! dpotrf, dpotri

INTERFACE 
   SUBROUTINE DPOTRF(UPLO,N,A,LDA,INFO)
     INTEGER(KIND=4) :: LDA
     CHARACTER(LEN=1) :: UPLO
     INTEGER(KIND=4) :: N
     REAL(KIND=8) :: A(LDA,*)
     INTEGER(KIND=4) :: INFO
   END SUBROUTINE DPOTRF
END INTERFACE

INTERFACE 
   SUBROUTINE DPOTRI(UPLO,N,A,LDA,INFO)
     INTEGER(KIND=4) :: LDA
     CHARACTER(LEN=1) :: UPLO
     INTEGER(KIND=4) :: N
     REAL(KIND=8) :: A(LDA,*)
     INTEGER(KIND=4) :: INFO
   END SUBROUTINE DPOTRI
END INTERFACE

INTERFACE 
   SUBROUTINE DPOTRS( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
     implicit none
     character(len=1), intent(in) :: uplo
     integer(kind=4), intent(in) :: n, nrhs, lda, ldb
     real(kind=8), intent(in) :: a(lda,n)
     real(kind=8), intent(inout) :: b(ldb,nrhs)
     integer(kind=4), intent(out) :: info
   END SUBROUTINE DPOTRS
END INTERFACE


INTERFACE
   SUBROUTINE DSYTRS( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
     implicit none
     character(len=1), intent(in) :: uplo
     integer(kind=4), intent(in) :: n, nrhs, lda, ldb
     integer(kind=4), intent(in) :: ipiv(n)
     real(kind=8), intent(in) :: a(lda,n)
     real(kind=8), intent(inout) :: b(ldb,nrhs)
     integer(kind=4), intent(out) :: info
     END SUBROUTINE DSYTRS
END INTERFACE


interface
   subroutine dgesvd(jobu,jobvt,m,n,a,lda,s,u,ldu,vt,ldvt,work,lwork,info)
     character          jobu, jobvt
     integer            info, lda, ldu, ldvt, lwork, m, n
     double precision   a(lda,*), s(*), u(ldu,*),vt(ldvt,*),work(*)
   end subroutine dgesvd

!!! symmetric linear solver
     !! C = alpha*A*A'+beta*C
     subroutine dsyrk(uplo,trans,n,k,alpha,a,lda,beta,c,ldc)
       implicit none
       character(len=1), intent(in) :: uplo, trans
       integer(kind=4), intent(in) :: n, k, lda, ldc
       real(kind=8), intent(in) :: alpha, beta
       real(kind=8), intent(in) :: a(lda,*)
       real(kind=8), intent(inout) :: c(ldc,*)
     end subroutine dsyrk
     !! B = A\B, A general n*n
     subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
       implicit none
       integer(kind=4), intent(in) :: n, nrhs, lda, ldb
       real(kind=8), intent(inout) :: a(lda,n), b(ldb,n)
       integer(kind=4), intent(out) :: ipiv(n)
       integer(kind=4), intent(out) :: info
     end subroutine dgesv
     !! B = A\B, A symmetric n*n
     subroutine dposv(uplo, n, nrhs, a, lda, b, ldb, info) 
       implicit none
       character(len=1), intent(in) :: uplo
       integer(kind=4), intent(in) :: n, nrhs, lda, ldb
       real(kind=8), intent(inout) :: a(lda,n), b(ldb,n)
       integer(kind=4), intent(out) :: info
     end subroutine dposv

end interface

interface
   subroutine dgels(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info)
     implicit none
     character(len=1), intent(in) :: trans
     integer(kind=4), intent(in) :: m, n, nrhs, lda, ldb, lwork
     real(kind=8), intent(inout) :: a(lda,n), b(ldb,nrhs)
     real(kind=8), intent(inout) :: work(lwork)
     integer(kind=4), intent(out) :: info     
   end subroutine dgels
end interface

!DEC$ ENDIF
