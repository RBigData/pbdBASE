#include <R.h>
#include <Rinternals.h>
#include "base_global.h"

SEXP R_BCKMEANS(SEXP K, SEXP X, SEXP DESCX, SEXP IMAX, 
{

K, X, DESCX, MU, DESCMU, Z, IMAX, LDIM1, LDIM2, 
     $                    INIMHD)

{
  const int IJ = 1;
  double* cpA;
  int i, info = 0;
  int* pt_ALDIM = INTEGER(ALDIM);
  
  int INIMHD = 1;
  
  SEXP RET;
  PROTECT(RET = allocVector(REALSXP, 2));
  
  // compute inverse of condition number
  F77_CALL(condnum)(CHARPT(TYPE, 0), INTEGER(M), INTEGER(N), cpA, 
    &IJ, &IJ, INTEGER(DESCA), REAL(RET), &info);
  
  REAL(RET)[1] = (double) info;
  
  UNPROTECT(1);
  return(RET);
}
