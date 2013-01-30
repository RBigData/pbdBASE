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
  
  F77_CALL(pdgemm)(CHARPT(TRANSA, 0), CHARPT(TRANSB, 0),
    INTEGER(M), INTEGER(N), INTEGER(K), &alpha, 
    REAL(A), &IJ, &IJ, INTEGER(DESCA),
    REAL(B), &IJ, &IJ, INTEGER(DESCB), &beta,
    REAL(C), &IJ, &IJ, INTEGER(DESCC));
  
  UNPROTECT(1);
  return(C);
}

// PDSYRK
SEXP R_PDSYRK(SEXP UPLO, SEXP TRANS, SEXP N, SEXP K, 
  SEXP A, SEXP DESCA, SEXP CLDIM, SEXP DESCC)
{
  SEXP C;
  PROTECT(C = allocMatrix(REALSXP, INTEGER(CLDIM)[0], INTEGER(CLDIM)[1]));
  
  const double alpha = 1.0;
  const double beta = 0.0;
  const int IJ = 1;
  const int size = INTEGER(N)[0];
  const char trans = 'T';
  const char uplo = 'U';
  
  pdsyrk_(CHARPT(UPLO, 0), CHARPT(TRANS, 0),
    &size, INTEGER(K), &alpha, 
    REAL(A), &IJ, &IJ, INTEGER(DESCA), &beta,
    REAL(C), &IJ, &IJ, INTEGER(DESCC));
  
  // Copy over the upper triangle to the lower
  pdmksym_(&uplo, REAL(C), &IJ, &IJ, INTEGER(DESCC));
  
  UNPROTECT(1);
  return(C);
}


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


