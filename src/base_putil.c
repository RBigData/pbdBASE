/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

// Copyright 2013, 2016 Schmidt

#include "pbdBASE.h"


SEXP R_MKSUBMAT(SEXP GBLX, SEXP LDIM, SEXP DESCX)
{
  SEXP SUBX;
  PROTECT(SUBX = allocMatrix(REALSXP, INTEGER(LDIM)[0], INTEGER(LDIM)[1]));
  
  mksubmat_(REAL(GBLX), REAL(SUBX), INTEGER(DESCX));
  
  UNPROTECT(1);
  return SUBX;
} 


SEXP R_MKGBLMAT(SEXP SUBX, SEXP DESCX, SEXP RDEST, SEXP CDEST)
{
  SEXP GBLX;
  PROTECT(GBLX = allocMatrix(REALSXP, INTEGER(DESCX)[2], INTEGER(DESCX)[3]));
  
  mkgblmat_(REAL(GBLX), REAL(SUBX), INTEGER(DESCX), INTEGER(RDEST), 
    INTEGER(CDEST));
  
  UNPROTECT(1);
  return GBLX;
} 


SEXP R_DALLREDUCE(SEXP X, SEXP LDIM, SEXP DESCX, SEXP OP, SEXP SCOPE)
{
  const int m = INTEGER(DESCX)[2], n = INTEGER(DESCX)[3];
  
  SEXP CPX;
  PROTECT(CPX = allocMatrix(REALSXP, INTEGER(LDIM)[0], INTEGER(LDIM)[1]));
  
  memcpy(REAL(CPX), REAL(X), m*n*sizeof(double));
  
  dallreduce_(REAL(CPX), INTEGER(DESCX), CHARPT(OP, 0), CHARPT(SCOPE, 0));
  
  UNPROTECT(1);
  return CPX;
} 


SEXP R_PTRI2ZERO(SEXP UPLO, SEXP DIAG, SEXP X, SEXP LDIM, SEXP DESCX)
{
  const int m = INTEGER(LDIM)[0], n = INTEGER(LDIM)[1];
  
  SEXP CPX;
  PROTECT(CPX = allocMatrix(REALSXP, m, n));
  
  memcpy(REAL(CPX), REAL(X), m*n*sizeof(double));
  
  ptri2zero_(CHARPT(UPLO, 0), CHARPT(DIAG, 0), REAL(CPX), INTEGER(DESCX));
  
  UNPROTECT(1);
  return CPX;
}


SEXP R_PDSWEEP(SEXP X, SEXP LDIM, SEXP DESCX, SEXP VEC, SEXP LVEC, SEXP MARGIN, SEXP FUN)
{
  const int m = INTEGER(LDIM)[0], n = INTEGER(LDIM)[1];
  int IJ = 1;
  
  SEXP CPX;
  PROTECT(CPX = allocMatrix(REALSXP, m, n));
  
  memcpy(REAL(CPX), REAL(X), m*n*sizeof(double));
  
  pdsweep(REAL(CPX), IJ, IJ, INTEGER(DESCX), REAL(VEC), INT(LVEC), INT(MARGIN), CHARPT(FUN, 0)[0]);
  
  UNPROTECT(1);
  return CPX;
}


SEXP R_PDGDGTK(SEXP X, SEXP LDIM, SEXP DESCX, SEXP LDIAG, SEXP RDEST, SEXP CDEST)
{
  int IJ = 1;
  
  SEXP DIAG;
  PROTECT(DIAG = allocVector(REALSXP, INTEGER(LDIAG)[0]));
  
  pdgdgtk_(REAL(X), &IJ, &IJ, INTEGER(DESCX), REAL(DIAG), 
    INTEGER(RDEST), INTEGER(CDEST));
  
  UNPROTECT(1);
  return DIAG;
}


SEXP R_PDDIAGMK(SEXP LDIM, SEXP DESCX, SEXP DIAG, SEXP LDIAG)
{
  const int m = INTEGER(LDIM)[0], n = INTEGER(LDIM)[1];
  int IJ = 1;
  
  SEXP X;
  PROTECT(X = allocMatrix(REALSXP, m, n));
  
  pddiagmk_(REAL(X), &IJ, &IJ, INTEGER(DESCX), REAL(DIAG), INTEGER(LDIAG));
  
  UNPROTECT(1);
  return X;
}



SEXP R_DHILBMK(SEXP N)
{
  int n = INTEGER(N)[0];
  
  SEXP X;
  PROTECT(X = allocMatrix(REALSXP, n, n));
  
  dhilbmk_(&n, REAL(X));
  
  UNPROTECT(1);
  return X;
}



SEXP R_PDHILBMK(SEXP LDIM, SEXP DESCX)
{
  SEXP X;
  PROTECT(X = allocMatrix(REALSXP, INTEGER(LDIM)[0], INTEGER(LDIM)[1]));
  
  pdhilbmk_(REAL(X), INTEGER(DESCX));
  
  UNPROTECT(1);
  return X;
}



SEXP R_PDMKCPN1(SEXP LDIM, SEXP DESCX, SEXP COEF)
{
  const int m = INTEGER(LDIM)[0], n = INTEGER(LDIM)[1];
  
  SEXP X;
  PROTECT(X = allocMatrix(REALSXP, m, n));
  
  pdmkcpn1_(REAL(X), INTEGER(DESCX), REAL(COEF));
  
  UNPROTECT(1);
  return X;
}
