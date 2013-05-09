#include <R.h>
#include <Rinternals.h>
#include "base_global.h"

SEXP R_PGLM_FIT(SEXP FAMILY, SEXP LINK, SEXP INCPT, SEXP STOPRULE,
  SEXP X, SEXP DESCX, SEXP MAXITER, SEXP TOL)
{
  int i, *pt_XLDIM = INTEGER(XALDIM), *pt_YLDIM = INTEGER(YLDIM);
  
  double *pt_ORG, *pt_COPY, *pt_FT, *p_work;
  const double tmp = 0;
  double work = 0;
  
  SEXP RET, RET_NAMES, INFO, A_OUT, B_OUT, 
       EFF, FT, RSD, TAU, IPIV, RANK;
  
  /* set up return */
  PROTECT(RET = allocVector(VECSXP, 9));
  PROTECT(RET_NAMES = allocVector(STRSXP, 9));
  
  PROTECT(INFO = allocVector(INTSXP, 1));
  PROTECT(A_OUT = allocMatrix(REALSXP, pt_ALDIM[0], pt_ALDIM[1]));
  PROTECT(B_OUT = allocMatrix(REALSXP, pt_BLDIM[0], pt_BLDIM[1]));
  PROTECT(EFF = allocMatrix(REALSXP, pt_BLDIM[0], pt_BLDIM[1]));
  PROTECT(FT = allocMatrix(REALSXP, pt_BLDIM[0], pt_BLDIM[1]));
  PROTECT(RSD = allocMatrix(REALSXP, pt_BLDIM[0], pt_BLDIM[1]));
  PROTECT(TAU = allocVector(REALSXP, INTEGER(LTAU)[0]));
  PROTECT(IPIV = allocVector(INTSXP, pt_ALDIM[1]));
  PROTECT(RANK = allocVector(INTSXP, 1));
  
  SET_VECTOR_ELT(RET, 0, INFO);
  SET_VECTOR_ELT(RET, 1, A_OUT);
  SET_VECTOR_ELT(RET, 2, B_OUT);
  SET_VECTOR_ELT(RET, 3, EFF);
  SET_VECTOR_ELT(RET, 4, FT);
  SET_VECTOR_ELT(RET, 5, RSD);
  SET_VECTOR_ELT(RET, 6, TAU);
  SET_VECTOR_ELT(RET, 7, IPIV);
  SET_VECTOR_ELT(RET, 8, RANK);
  
  SET_STRING_ELT(RET_NAMES, 0, mkChar("INFO")); 
  SET_STRING_ELT(RET_NAMES, 1, mkChar("A")); 
  SET_STRING_ELT(RET_NAMES, 2, mkChar("B")); 
  SET_STRING_ELT(RET_NAMES, 3, mkChar("EFF")); 
  SET_STRING_ELT(RET_NAMES, 4, mkChar("FT")); 
  SET_STRING_ELT(RET_NAMES, 5, mkChar("RSD")); 
  SET_STRING_ELT(RET_NAMES, 6, mkChar("TAU")); 
  SET_STRING_ELT(RET_NAMES, 7, mkChar("IPIV"));
  SET_STRING_ELT(RET_NAMES, 8, mkChar("RANK")); 
  
  setAttrib(RET, R_NamesSymbol, RET_NAMES);
  
  
  SEXP BETA, WT, RESIDS, INFO;
  
  
  
  pglm_fit_(CHARPT(FAMILY, 0), CHARPT(LINK, 0), CHARPT(INCPT, 0), CHARPT(STOPRULE, 0), 
    REAL(X), INTEGER(DESCX), REAL(Y), INTEGER(DESCY), REAL(BETA), 
    REAL(WT), REAL(RESIDS), INTEGER(MAXITER), INTEGER(INFO), REAL(TOL));
  
  
}
