#include "base_global.h"
#include "base/expm/matexp.h"


SEXP R_matpow_by_squaring(SEXP A, SEXP b)
{
  const int n = nrows(A);
  double *cpA;
  
  SEXP P;
  PROTECT(P = allocMatrix(REALSXP, n, n));
  
  // A is modified
  cpA = malloc(n*n*sizeof(double));
  memcpy(cpA, REAL(A), n*n*sizeof(double));
  
  matpow_by_squaring(cpA, n, INT(b,0), REAL(P));
  
  free(cpA);
  
  UNPROTECT(1);
  return(P);
}



SEXP R_matexp_pade(SEXP A)
{
  const int n = nrows(A);
  SEXP N, D;
  SEXP RET, RET_NAMES;
  
  // Allocate N and D
  PROTECT(N = allocMatrix(REALSXP, n, n));
  PROTECT(D = allocMatrix(REALSXP, n, n));
  
  // Compute N and D
  matexp_pade(n, REAL(A), REAL(N), REAL(D));
  
  // Wrangle the return
  PROTECT(RET = allocVector(VECSXP, 2));
  PROTECT(RET_NAMES = allocVector(STRSXP, 2));
  
  SET_VECTOR_ELT(RET, 0, N);
  SET_VECTOR_ELT(RET, 1, D);
  
  SET_STRING_ELT(RET_NAMES, 0, mkChar("N")); 
  SET_STRING_ELT(RET_NAMES, 1, mkChar("D")); 
  
  setAttrib(RET, R_NamesSymbol, RET_NAMES);
  
  
  UNPROTECT(4);
  return(RET);
}


