!!!
!!! Calculating the Kalman filter likelihood
!!! Matlab call: dlmmex(y,F,V,x0,G,W,C0,XX)
!!! returns the -2*log(likelihood)
!!!
!!! compile: mex -O dlmmex.F90
!!!
!!! Need to have LAPACK/BLAS library mex linker option set, e.g.
!!!  -framework Accelerate in MAC.
!!!
!!! OBS: No error checking until matlab/fortran/mac/mex R2014b
!!! is fixed. Currently mexErrMsgTxt breaks the matlab session.
!!!
!!! Marko Laine <marko.laine@fmi.fi>
!!! $Revision: 0.0 $  $Date: 2014/12/29 12:00:00 $
!!!
!!! This file comes with no support, you may use it to speed up the
!!!  dlm computatios if you are fluent with fortran language mex files.
!!!

#include "fintrf.h"

#if defined(INT64)
#define INTLEN 8
#else
#define INTLEN 4
#endif

subroutine mexFunction(nlhs, plhs, nrhs, prhs)
#if defined MSWIND
  ! For Windows only!
  ! This is the library needed to set the floating point 
  ! exception word and it needs to be here in the code.  
  use msflib
#endif

  implicit none
#include "mex_inc.F90"

  ! input variables
  integer*4, intent(in) :: nlhs, nrhs
  mwPointer, intent(in) :: prhs(*)
  mwPointer, intent(out):: plhs(*)

  !! input y,F,V,x0,G,W,C0
  real(kind=8), allocatable :: y(:,:), F(:,:), V(:,:), x0(:), G(:,:), W(:,:), C0(:,:), XX(:,:)

  !! output lik
  real(kind=8) lik(1)

  !!
  mwPointer :: ptr
  mwSize :: n, m, one = 1

  interface
     function dlmlik(y,F,V,x0,G,W,C0,XX) result(lik)
       implicit none
       real(kind=8), intent(in) :: y(:,:), F(:,:), V(:,:), x0(:), G(:,:), W(:,:), C0(:,:), XX(:,:)
       real(kind=8) :: lik
     end function dlmlik
  end interface

  if (nrhs .lt. 8) then
     call mexPrintf('Error, 8 and only 8 input arguments')
     return
  elseif (nlhs .lt. 1) then
     call mexPrintf('Error, needs 1 output argument')
     return
  end if
  
  !!  y(:,:),
  ptr = prhs(1)
  m = mxGetM(ptr)
  n = mxGetN(ptr)
  allocate(y(m,n))
  call mxCopyPtrToReal8(mxGetPr(ptr),y,m*n)
  !! F(:,:),
  ptr = prhs(2)
  m = mxGetM(ptr)
  n = mxGetN(ptr)
  allocate(F(m,n))
  call mxCopyPtrToReal8(mxGetPr(ptr),F,m*n)
  !! V(:,:),
  ptr = prhs(3)
  m = mxGetM(ptr)
  n = mxGetN(ptr)
  allocate(V(m,n))
  call mxCopyPtrToReal8(mxGetPr(ptr),V,m*n)
  !! x0(:),
  ptr = prhs(4)
  m = mxGetM(ptr)
  n = mxGetN(ptr)
  allocate(x0(m*n))
  call mxCopyPtrToReal8(mxGetPr(ptr),x0,m*n)
  !! G(:,:),
  ptr = prhs(5)
  m = mxGetM(ptr)
  n = mxGetN(ptr)
  allocate(G(m,n))
  call mxCopyPtrToReal8(mxGetPr(ptr),G,m*n)
  !! W(:,:),
  ptr = prhs(6)
  m = mxGetM(ptr)
  n = mxGetN(ptr)
  allocate(W(m,n))
  call mxCopyPtrToReal8(mxGetPr(ptr),W,m*n)
  !! C0(:,:)
  ptr = prhs(7)
  m = mxGetM(ptr)
  n = mxGetN(ptr)
  allocate(C0(m,n))
  call mxCopyPtrToReal8(mxGetPr(ptr),C0,m*n)
  !! XX(:,:)
  ptr = prhs(8)
  m = mxGetM(ptr)
  n = mxGetN(ptr)
  allocate(XX(m,n))
  call mxCopyPtrToReal8(mxGetPr(ptr),XX,m*n)
  
  lik(1) = dlmlik(y,F,V,x0,G,W,C0,XX)

  plhs(1) = mxCreateDoubleMatrix(one,one, 0);
  call mxCopyReal8ToPtr(lik, mxGetPr(plhs(1)), one)

  deallocate(y, F, V, x0, G, W, C0, XX)
    
end subroutine mexFunction

!!
!! calculates the likelihood, see dlmsmo.m matlab code
!! used solvelin for Kalman filtering step
!!
function dlmlik(y,F0,V,x0,G,W,C0,XX) result(lik)
  use matlab
  implicit none
  real(kind=8), intent(in) :: y(:,:), F0(:,:), V(:,:), x0(:), G(:,:), W(:,:), C0(:,:), XX(:,:)
  ! ,Xproxy(:,:)
  real(kind=8) lik
  real(kind=8), allocatable :: rr(:), rr2(:), Cpp(:,:), FF(:,:)

  integer :: p,m,n,i,j,jj, info, nbad, ngood

  real(kind=8) :: x(size(G,1))
  real(kind=8) :: C(size(G,1),size(G,1))
  real(kind=8) :: r(size(F0,1)), r2(size(F0,1))
  real(kind=8) :: Cp(size(F0,1),size(F0,1))

  real(kind=8) :: F(size(F0,1),size(F0,2)+size(XX,2))

  logical :: igood(size(F,1))

  interface
     subroutine solvelin(Model,Yobs,Cobs,xprior,Cprior)
       implicit none
       real(8), intent(in) :: Model(:,:)
       real(8), intent(inout) :: Yobs(:), Cobs(:,:), xprior(:), Cprior(:,:)
     end subroutine solvelin
  end interface

  p = size(F0,1)
  m = size(F,2)
  n = size(y,1)
  F(:,1:size(F0,2)) = F0

  !! initialization
  x = x0
  C = C0
  r = 0.0
  lik = 0.0
  !! loop over observations
  do i=1,n
     !! add proxy variables to FX
     F(:,size(F0,2)+1:size(F,2)) = spread(XX(i,:),1,size(F,1))
     !! save obs error covariance in Cp
     Cp = 0.0 
     do j=1,p
        if (.not.isnan(y(i,j))) Cp(j,j) = V(i,j)**2
     end do
     !! observation vector in r
     r = y(i,:)
     ! need to remove missing values from y, F, and V (i.e. r,F,Cp)
     igood = .not.isnan(r)
     ngood = count(igood)
     nbad = p-ngood
     if (nbad > 0) then
        if (ngood>0) then
           allocate(rr(ngood),rr2(ngood),Cpp(ngood,ngood),FF(ngood,m))
           rr = pack(r,igood)
           jj = 1
!!!           FF = reshape(pack(F,spread(igood,p,1)),ngood,m)
!!!           Cpp = reshape(pack(Cp,spread(igood,p,1) ...),ngood,ngood)         
           do j=1,p
              if (igood(j)) then
                 FF(jj,:) = F(j,:)
                 Cpp(jj,:) = pack(Cp(j,:),igood)
                 jj = jj+1
              end if
           end do
           call solvelin(FF,rr,Cpp,x,C)
           rr2 = rr
           call dpotrs('l', ngood, 1, Cpp, ngood, rr2, ngood, info )
           if (info .ne. 0) write(*,*) 'dpotrs info:',info
           lik = lik + dot_product(rr2,rr)
           do j=1,ngood
              lik = lik + log(Cpp(j,j))*2.0
           end do
           deallocate(rr,rr2,Cpp,FF)
        else ! no observations
           ! do nothing
        end if
     else ! no missing observations
        !! solvelin replaces r with innovations,
        !! Cp with chol of Cp, the y prediction covariance, lower diagonal
        call solvelin(F,r,Cp,x,C)
        r2 = r
        !! r2 = r/Cp
        call dpotrs('l', p, 1, Cp, p, r2, p, info )
        if (info .ne. 0) write(*,*) 'dpotrs info:',info
        lik = lik + dot_product(r2,r)
        !! det(Cp) is sum of diagonals of the chol * 2
        do j=1,p
           lik = lik + log(Cp(j,j))*2.0
        end do
     end if

     !! propagate to next stage
     if (i<n) then
        x = matmul(G,x)
        C = matmul(matmul(G,C),transpose(G)) + W
        ! Not working, want to relplace x with G*x, possible with dgemv?
        !        call dgemv('n',m,m,1.0d0,G,m,x,1,0.0d0,x,1)
        !        call dsymm('r','l',m,m,1.0d0,C,m,G,m,0.0d0,C,m)
        !        call dgemm('n','t',m,m,m,1.0d0,C,m,G,m,0.0d0,C,m);
       
     end if
  end do
end function dlmlik

!!! Solve linear equation with Gaussian prior and known covariance
!!! y = Model x + e, e ~ N(0,Cobs),  with prior
!!! x ~ N(xprior,Cprior)
!!! y is n*1 and x is m*1, assume n<m and use Kalman formula.
!!!
!!! on input:
!!! Model  n*m model matrix
!!! Yobs   n*1 observations
!!! Cobs   n*n observation covariance
!!! xprior m*1 prior state
!!! Cprior m*m prior covariance
!!!
!!! on output: 
!!! Model in unchanged
!!! Yobs has model innovations Yobs-M*xprior
!!! Cobs has Cholesky factor of Model*Cprior*Model'+Cobs
!!! xprior has xest, mean estimate of the state
!!! Cprior has Cest, covariance of the state estimate
!!!
!!! marko.laine@fmi.fi 2012
subroutine solvelin(Model,Yobs,Cobs,xprior,Cprior)
  implicit none
  real(8), intent(in) :: Model(:,:)
  real(8), intent(inout) :: Yobs(:), Cobs(:,:), xprior(:), Cprior(:,:)
  !! work space
  real(8) :: Gt(size(Model,1),size(xprior,1)) ! n*m transpose of the Kalman Gain
  real(8) :: Cwork(size(Cprior,1),size(Cprior,2)) ! work space for Cest
  !!
  integer n,m,info,i
  include 'lapack_inc.f90'

  n = size(Model,1)
  m = size(Model,2)

  if (size(Yobs,1) .ne. n .or. &
       size(xprior,1) .ne. m .or. &
       any(shape(Cobs) .ne. (/n,n/)) .or. &
       any(shape(Cprior) .ne. (/m,m/))) then
     write(*,*) size(xprior)
     write(*,*) shape(Cobs)
     call mexPrintf('Error in shapes')
     return
     !   stop 'error in shapes'
  end if

  !! Calculate Kalman gain G = Cprior*Model'/(Model*Cprior*Model'+Cobs).
  ! G' = Model*Cprior
  call dsymm('r','l',n,m,1.0d0,Cprior,m,Model,n,0.0d0,Gt,n)
  ! Cobs = Model*Cprior*Model' + Cobs = G'*Model' + Cobs,
  call dgemm('n','t',n,n,m,1.0d0,Gt,n,Model,n,1.0d0,Cobs,n);
  ! G' = Cobs\G'
  ! DPOSV solves A\B = inv(A)*B, we need B/A = B*inv(A), so use (B\A)' = A'\B'
  call dposv('l',n,m,Cobs,n,Gt,n,info)
  !! Now Cobs has Cholesky of (Model*Cprior*Model'+Cobs) and 
  !! Gt is transpose of Kalman Gain.
  !! Calculate xest = xprior + G*(Yobs-M*xprior), xprior is replaced by xest.
  ! Yobs = Yobs - M*xprior, so Yobs has now the innovations.
  call dgemv('n',n,m,-1.0d0,Model,n,xprior,1,1.0d0,Yobs,1)
  ! xprior = xprior + G*Yobs
  call dgemv('t',n,m,1.0d0,Gt,n,Yobs,1,1.0d0,xprior,1)
  !! Calculate Cest = Cprior-G*Model*Cprior, Cprior is replaced by Cest.
  ! Cwork = G*Model
  call dgemm('t','n',m,m,n,1.0d0,Gt,n,Model,n,0.0d0,Cwork,m);
  ! Cprior = Cprior - Cwork*Cprior 
  call dsymm('r','l',m,m,-1.0d0,Cprior,m,Cwork,m,1.0d0,Cprior,m)
end subroutine solvelin
