#include "base_global.h"
#include "base/expm/matexp.h"

/*SEXP R_matexp_pade(SEXP n, SEXP A)*/
/*{*/
/*  SEXP N, D;*/
/*  SEXP RET, RET_NAMES;*/
/*  */
/*  // Allocate N and D*/
/*  PROTECT(N = allocMatrix(REALSXP, INT(n,0), INT(n,0)));*/
/*  PROTECT(D = allocMatrix(REALSXP, INT(n,0), INT(n,0)));*/
/*  */
/*  // Compute N and D*/
/*  matexp_pade(INT(n,0), REAL(A), REAL(N), REAL(D));*/
/*  */
/*  // Wrangle the return*/
/*  PROTECT(RET = allocVector(VECSXP, 2));*/
/*  PROTECT(RET_NAMES = allocVector(STRSXP, 2));*/
/*  */
/*  SET_VECTOR_ELT(RET, 0, N);*/
/*  SET_VECTOR_ELT(RET, 1, D);*/
/*  */
/*  SET_STRING_ELT(RET_NAMES, 0, mkChar("N")); */
/*  SET_STRING_ELT(RET_NAMES, 1, mkChar("D")); */
/*  */
/*  setAttrib(RET, R_NamesSymbol, RET_NAMES);*/
/*  */
/*  */
/*  UNPROTECT(4);*/
/*  return(RET);*/
/*}*/


SEXP R_matpow_by_squaring(SEXP A, SEXP desca, SEXP ldim, SEXP b)
{
  double *cpA;
  
  SEXP P;
  PROTECT(P = allocMatrix(REALSXP, INT(ldim, 0), INT(ldim, 1)));
  
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

