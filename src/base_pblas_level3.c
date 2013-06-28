#include <R.h>
#include <Rinternals.h>
#include "base_global.h"

// PDTRAN
SEXP R_PDTRAN(SEXP M, SEXP N, SEXP A, SEXP DESCA, SEXP CLDIM, SEXP DESCC)
{
  SEXP C;
  PROTECT(C = allocMatrix(REALSXP, INTEGER(CLDIM)[0], INTEGER(CLDIM)[1]));
  
  const int IJ = 1;
  const double one = 1.0;
  const double zero = 0.0;
  
  F77_CALL(pdtran)(INTEGER(M), INTEGER(N), &one, 
      REAL(A), &IJ, &IJ, INTEGER(DESCA), &zero, 
      REAL(C), &IJ, &IJ, INTEGER(DESCC));
  
  UNPROTECT(1);
  return(C);
}

// PDGEMM
SEXP R_PDGEMM(SEXP TRANSA, SEXP TRANSB, SEXP M, SEXP N, SEXP K,
  SEXP A, SEXP DESCA, SEXP B, SEXP DESCB, SEXP CLDIM, SEXP DESCC)
{
  SEXP C;
  PROTECT(C = allocMatrix(REALSXP, INTEGER(CLDIM)[0], INTEGER(CLDIM)[1]));
  
  const double alpha = 1.0;
  const double beta = 0.0;
  const int IJ = 1;
  
  pdgemm_(CHARPT(TRANSA, 0), CHARPT(TRANSB, 0),
    INTEGER(M), INTEGER(N), INTEGER(K), &alpha, 
    REAL(A), &IJ, &IJ, INTEGER(DESCA),
    REAL(B), &IJ, &IJ, INTEGER(DESCB), &beta,
    REAL(C), &IJ, &IJ, INTEGER(DESCC));
  
  UNPROTECT(1);
  return(C);
}

