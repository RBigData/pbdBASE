#include <R.h>
#include <Rinternals.h>
#include "../base_global.h"

SEXP R_PDCROSSPROD(SEXP TRANS, SEXP A, SEXP DESCA, SEXP CLDIM, SEXP DESCC)
{
  SEXP C;
  PROTECT(C = allocMatrix(REALSXP, INTEGER(CLDIM)[0], INTEGER(CLDIM)[1]));
  
  const double alpha = 1.0;
  const int IJ = 1;
  
  pdcrossprod_(CHARPT(TRANS, 0), &alpha, REAL(A), &IJ, &IJ, INTEGER(DESCA), 
    REAL(C), &IJ, &IJ, INTEGER(DESCC));
  
  UNPROTECT(1);
  return(C);
}

