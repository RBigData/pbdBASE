#include <R.h>
#include <Rinternals.h>
#include "base_global.h"

SEXP R_PDLAPRNT(SEXP M, SEXP N, SEXP A, SEXP DESCA, SEXP CMATNM, SEXP NOUT)
{
  double work[INTEGER(DESCA)[8]];
  const int IJ = 1;
  const int SRC = 0;
  
  F77_CALL(bprnt)(INTEGER(M), INTEGER(N), REAL(A), &IJ, &IJ,
    INTEGER(DESCA), &SRC, &SRC, CHARPT(CMATNM, 0),
    INTEGER(NOUT), &work);
  
  return(R_NilValue);
} /* End of R_PDLAPRNT(). */


SEXP R_PDGEMR2D(SEXP M, SEXP N, SEXP X, SEXP DESCX, SEXP CLDIM, SEXP DESCB, SEXP CTXT)
{
  const int IJ = 1;
  SEXP B;
  
  PROTECT(B = allocMatrix(REALSXP, INTEGER(CLDIM)[0], INTEGER(CLDIM)[1]));
  
  Cpdgemr2d(INTEGER(M)[0], INTEGER(N)[0],
    REAL(X), IJ, IJ, INTEGER(DESCX),
    REAL(B), IJ, IJ, INTEGER(DESCB), INTEGER(CTXT)[0]);
  
  UNPROTECT(1);
  return(B);
} /* End of R_PDGEMR2D(). */



// next best divisor function
SEXP R_nbd(SEXP N, SEXP D)
{
  int i, test;
  
  SEXP RET;
  PROTECT(RET = allocVector(INTSXP, 1));
  INTEGER(RET)[0] = INTEGER(D)[0];
  
  for (i=INTEGER(RET)[0]; i<=INTEGER(N)[0]; i++){
    test = INTEGER(N)[0] % i;
    if (test == 0){
      INTEGER(RET)[0] = i;
      break;
    }
  }
  
  UNPROTECT(1);
  return RET;
}



