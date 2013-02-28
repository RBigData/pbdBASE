#include <R.h>
#include <Rinternals.h>
#include "../base_global.h"

SEXP R_PDCLVAR(SEXP X, SEXP DESCX, SEXP LSD)
{
  SEXP VAR;
  PROTECT(VAR = allocVector(REALSXP, INTEGER(LSD)[0]));
  
  pdclvar_(REAL(X), INTEGER(DESCX), REAL(VAR));
  
  UNPROTECT(1);
  return(VAR);
} 
