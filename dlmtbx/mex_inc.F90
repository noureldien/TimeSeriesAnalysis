!!! Interfaces to Matlab mex utilities

  !! The following must be in the main program already
  !! #include "fintrf.h"
  
  interface

     function mxGetM(pm)
       mwPointer pm
       mwPointer mxGetM
     end function mxGetM

     function mxGetN(pm)
       mwPointer pm
       mwPointer mxGetN
     end function mxGetN

     function mxGetNumberOfElements(pm)
       mwPointer pm
       mwPointer mxGetNumberOfElements
     end function mxGetNumberOfElements
     
     function mxGetData(pm)
       mwPointer pm
       mwPointer mxGetData
     end function mxGetData

     function mexGetVariablePtr(varname, workspace)
       character*(*) varname, workspace
       mwPointer mexGetVariablePtr
     end function mexGetVariablePtr

     subroutine mxCopyPtrToReal8(px, y, n)
       mwPointer px
       mwSize n
       real*8 y(n)
     end subroutine mxCopyPtrToReal8

     function mxGetPr(pm)
       mwPointer pm
       mwPointer mxGetPr
      end function mxGetPr
      
     function mxIsFinite(value)
       real*8 value
       integer*4 mxIsFinite
     end function mxIsFinite

     function mxIsNaN(value)
       real*8 value
       integer*4 mxIsNaN
     end function mxIsNaN

     function mxIsDouble(pm)
       mwPointer pm
       integer*4 mxIsDouble
     end function mxIsDouble
     
     subroutine mxCopyPtrToInteger4(px, y, n)
       mwPointer px
       mwSize n
       integer*4 y(n)
     end subroutine mxCopyPtrToInteger4

     subroutine mxCopyReal8ToPtr(y, px, n)
       mwSize n 
       mwPointer px
       real*8 y(n)
     end subroutine mxCopyReal8ToPtr

     function mxCreateDoubleMatrix(m, n, ComplexFlag)
       mwSize m, n
       integer*4 ComplexFlag
       mwPointer mxCreateDoubleMatrix
     end function mxCreateDoubleMatrix

     function mxCreateDoubleScalar(value)
       real*8 value
       mwPointer mxCreateDoubleScalar
     end function mxCreateDoubleScalar

!     function mxCreateScalarDouble(value)
!       real*8 value
!       integer*4 mxCreateScalarDouble
!     end function mxCreateScalarDouble
     
     function mxCreateString(str)
       character*(*) str
       mwPointer mxCreateString
     end function mxCreateString

     subroutine mxDestroyArray(pm)
       mwPointer pm
     end subroutine mxDestroyArray

     function mxCreateStructMatrix(m, n, nfields, fieldnames)
       mwSize m, n
       integer*4 nfields
       character*(*) fieldnames(nfields)
       mwPointer mxCreateStructMatrix
     end function mxCreateStructMatrix
     
     function mxIsCell(pm)
       mwPointer pm
       integer*4 mxIsCell
     end function mxIsCell

     function mxGetCell(pm,index)
       mwPointer pm
       mwIndex index
       mwPointer mxGetCell
     end function mxGetCell

     function mxGetField(pm,index,fieldname)
       mwPointer pm
       mwIndex index
       character*(*) fieldname
       mwPointer mxGetField
     end function mxGetField

     subroutine mexErrMsgTxt(msg)
       character*(*) msg
     end subroutine mexErrMsgTxt

     subroutine mexErrMsgIdAndTxt(errorid, errormsg)
       character*(*) errorid, errormsg
     end subroutine mexErrMsgIdAndTxt
     
     subroutine mexPrintf(msg)
       character*(*) msg
     end subroutine mexPrintf

     subroutine mexSetTrapFlag(trapflag)
       integer*4 trapflag
     end subroutine mexSetTrapFlag

     function mexCallMATLAB(nlhs, plhs, nrhs, prhs, name)
       integer*4 nlhs, nrhs
       mwPointer plhs(*), prhs(*)
       character*(*) name
       integer*4 mexCallMATLAB
     end function mexCallMATLAB

     function mxAddField(pm, fieldname)
       mwPointer pm
       character*(*) fieldname
       integer*4 mxAddField
     end function mxAddField

     subroutine mxSetFieldByNumber(pm, index, fieldnumber, value)
       mwPointer pm, pvalue
       mwIndex index
       integer*4 fieldnumber
     end subroutine mxSetFieldByNumber

  end interface
