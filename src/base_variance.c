#include <R.h>
#include <Rinternals.h>
#include "base_global.h"

SEXP R_DDMATVAR(SEXP X, SEXP M, SEXP LCM, SEXP LCN, SEXP ICTXT)
{
  SEXP VAR;
  PROTECT(VAR = allocVector(REALSXP, INTEGER(LCN)[0]));
  
  F77_CALL(ddmatvar)(REAL(X), 
    INTEGER(M), INTEGER(LCM), INTEGER(LCN),
    INTEGER(ICTXT), REAL(VAR));
  
  UNPROTECT(1);
  return(VAR);
} 
