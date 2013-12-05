#include "base_global.h"
#include "base/expm/matexp.h"


SEXP R_p_matpow_by_squaring(SEXP A, SEXP desca, SEXP b)
{
  double *cpA;
  
  SEXP P;
  PROTECT(P = allocMatrix(REALSXP, nrows(A), ncols(A)));
  
  // Why did I make a copy ... ?
/*  cpA = malloc(N*N*sizeof(double));*/
/*  memcpy(cpA, REAL(A), N*N*sizeof(double));*/
  
  p_matpow_by_squaring(REAL(A), INTEGER(desca), INT(b, 0), REAL(P));
  
/*  free(cpA);*/
  
  UNPROTECT(1);
  return(P);
}




SEXP R_p_mateye(SEXP desca, SEXP ldim)
{
  SEXP RET;
  PROTECT(RET = allocMatrix(REALSXP, INT(ldim, 0), INT(ldim, 1)));
  
  p_mateye(REAL(RET), INTEGER(desca));
  
  UNPROTECT(1);
  return RET;
}




SEXP R_p_matexp_pade(SEXP A, SEXP desca)
{
  int m, n;
  SEXP N, D;
  SEXP RET, RET_NAMES;
  
  m = nrows(A);
  n = ncols(A);
  
  // Allocate N and D
  PROTECT(N = allocMatrix(REALSXP, m, n));
  PROTECT(D = allocMatrix(REALSXP, m, n));
  
  // Compute N and D
  p_matexp_pade(REAL(A), INTEGER(desca), REAL(N), REAL(D));
  
  
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
