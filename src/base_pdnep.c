#include <R.h>
#include <Rinternals.h>
#include "base_global.h"

SEXP R_PDNEP(SEXP X, SEXP DESCX, SEXP XLDIM)
{
  const int N = INTEGER(DESCX)[2];
  const int IJ = 1;
  const int m = INTEGER(XLDIM)[0], n = INTEGER(XLDIM)[1];
  double *CPX;
  
  SEXP WR, WI, INFO;
  PROTECT(WR = allocVector(REALSXP, N));
  PROTECT(WI = allocVector(REALSXP, N));
  PROTECT(INFO = allocVector(INTSXP, 1));
  
  CPX = R_alloc(m*n, sizeof(double));
  memcpy(CPX, REAL(X), m*n*sizeof(double));
  
  INTEGER(INFO)[0] = 0;
  
  pdnep_(CPX, &IJ, &IJ, INTEGER(DESCX), REAL(WR), REAL(WI), INTEGER(INFO));
  
  if (INTEGER(INFO)[0] != 0)
    Rprintf("INFO = %d\n", INTEGER(INFO)[0]);
  
  UNPROTECT(3);
  return WR;
}

