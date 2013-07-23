#include <R.h>
#include <Rinternals.h>
#include "../base_global.h"

SEXP R_PDCROSSPROD(SEXP UPLO, SEXP TRANS, SEXP A, SEXP DESCA, SEXP CLDIM, SEXP DESCC)
{
  SEXP C;
  PROTECT(C = allocMatrix(REALSXP, INTEGER(CLDIM)[0], INTEGER(CLDIM)[1]));
  
  const double alpha = 1.0;
  const int IJ = 1;
  
  pdcrossprod_(CHARPT(UPLO, 0), CHARPT(TRANS, 0), &alpha, 
    REAL(A), &IJ, &IJ, INTEGER(DESCA), 
    REAL(C), &IJ, &IJ, INTEGER(DESCC));
  
  UNPROTECT(1);
  return(C);
}


SEXP R_PDCHTRI(SEXP UPLO, SEXP A, SEXP ALDIM, SEXP DESCA, SEXP CLDIM, SEXP DESCC)
{
  const int IJ = 1;
  const int m = INTEGER(ALDIM)[0], n = INTEGER(ALDIM)[1];
  double *CPA;
  
  SEXP C, INFO;
  PROTECT(C = allocMatrix(REALSXP, INTEGER(CLDIM)[0], INTEGER(CLDIM)[1]));
  PROTECT(INFO = allocVector(INTSXP, 1));
  
  CPA = R_alloc(m*n, sizeof(double));
  memcpy(CPA, REAL(A), m*n*sizeof(double));
  
  INTEGER(INFO)[0] = 0;
  
  pdchtri_(CHARPT(UPLO, 0), CPA, &IJ, &IJ, INTEGER(DESCA), 
    REAL(C), &IJ, &IJ,
    INTEGER(DESCC), INTEGER(INFO));
  
  if (INTEGER(INFO)[0] != 0)
    Rprintf("INFO = %d\n", INTEGER(INFO)[0]);
  
  UNPROTECT(2);
  return(C);
}


SEXP R_PDGEEIG(SEXP X, SEXP DESCX, SEXP XLDIM)
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
  
  pdgeeig_(CPX, &IJ, &IJ, INTEGER(DESCX), REAL(WR), REAL(WI), INTEGER(INFO));
  
  if (INTEGER(INFO)[0] != 0)
    Rprintf("INFO = %d\n", INTEGER(INFO)[0]);
  
  UNPROTECT(3);
  return WR;
}

