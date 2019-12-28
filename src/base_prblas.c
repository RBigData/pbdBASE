/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

// Copyright 2013, Schmidt

#include <string.h>
#include "base/linalg/linalg.h"

// R.h and Rinternals.h needs to be included after Rconfig.h
#include "pbdBASE.h"


SEXP R_RL2BLAS(SEXP X, SEXP LDIM, SEXP DESCX, SEXP VEC, SEXP LVEC, SEXP FUN)
{
  const int m = INTEGER(LDIM)[0], n = INTEGER(LDIM)[1];
  
  SEXP CPX;
  PROTECT(CPX = allocMatrix(REALSXP, m, n));
  
  memcpy(REAL(CPX), REAL(X), m*n*sizeof(double));
  
  rl2blas_(REAL(CPX), INTEGER(DESCX), REAL(VEC), INTEGER(LVEC), INTEGER(FUN));
  
  UNPROTECT(1);
  return CPX;
}


SEXP R_RL2INSERT(SEXP X, SEXP LDIM, SEXP DESCX, SEXP VEC, SEXP LVEC, SEXP INDI, SEXP LINDI, SEXP INDJ, SEXP LINDJ)
{
  const int m = INTEGER(LDIM)[0], n = INTEGER(LDIM)[1];
  
  SEXP CPX;
  PROTECT(CPX = allocMatrix(REALSXP, m, n));
  
  memcpy(REAL(CPX), REAL(X), m*n*sizeof(double));
  
  rl2insert_(REAL(CPX), INTEGER(DESCX), REAL(VEC), INTEGER(LVEC), 
    INTEGER(INDI), INTEGER(LINDI), INTEGER(INDJ), INTEGER(LINDJ));
  
  UNPROTECT(1);
  return CPX;
}


SEXP R_RCOLCPY(SEXP X, SEXP LDIM, SEXP DESCX, SEXP XCOL, SEXP Y, SEXP DESCY, SEXP YCOL, SEXP LCOLS)
{
  const int m = INTEGER(LDIM)[0], n = INTEGER(LDIM)[1];
  
  SEXP CPX;
  PROTECT(CPX = allocMatrix(REALSXP, m, n));
  
  memcpy(REAL(CPX), REAL(X), m*n*sizeof(double));
  
  rcolcpy_(REAL(CPX), INTEGER(DESCX), INTEGER(XCOL), REAL(Y), 
    INTEGER(DESCY), INTEGER(YCOL), INTEGER(LCOLS));
  
  UNPROTECT(1);
  return CPX;
}


SEXP R_RCOLCPY2(SEXP X, SEXP LDIM, SEXP DESCX, SEXP XCOL, SEXP LXCOLS, SEXP Y, SEXP DESCY, SEXP YCOL, SEXP LYCOLS)
{
  const int m = INTEGER(LDIM)[0], n = INTEGER(LDIM)[1];
  
  SEXP CPX;
  PROTECT(CPX = allocMatrix(REALSXP, m, n));
  
  memcpy(REAL(CPX), REAL(X), m*n*sizeof(double));
  
  rcolcpy2_(REAL(CPX), INTEGER(DESCX), INTEGER(XCOL), INTEGER(LXCOLS), 
    REAL(Y), INTEGER(DESCY), INTEGER(YCOL), INTEGER(LYCOLS));
  
  UNPROTECT(1);
  return CPX;
}



SEXP R_RROWCPY(SEXP X, SEXP LDIM, SEXP DESCX, SEXP XROW, SEXP Y, SEXP DESCY, SEXP YROW, SEXP LROWS)
{
  const int m = INTEGER(LDIM)[0], n = INTEGER(LDIM)[1];
  
  SEXP CPX;
  PROTECT(CPX = allocMatrix(REALSXP, m, n));
  
  memcpy(REAL(CPX), REAL(X), m*n*sizeof(double));
  
  rrowcpy_(REAL(CPX), INTEGER(DESCX), INTEGER(XROW), REAL(Y), 
    INTEGER(DESCY), INTEGER(YROW), INTEGER(LROWS));
  
  UNPROTECT(1);
  return CPX;
}

SEXP R_RROWCPY2(SEXP X, SEXP LDIM, SEXP DESCX, SEXP XROW, SEXP LXROWS, SEXP Y, SEXP DESCY, SEXP YROW, SEXP LYROWS)
{
  const int m = INTEGER(LDIM)[0], n = INTEGER(LDIM)[1];
  
  SEXP CPX;
  PROTECT(CPX = allocMatrix(REALSXP, m, n));
  
  memcpy(REAL(CPX), REAL(X), m*n*sizeof(double));
  
  rrowcpy2_(REAL(CPX), INTEGER(DESCX), INTEGER(XROW), INTEGER(LXROWS), 
    REAL(Y), INTEGER(DESCY), INTEGER(YROW), INTEGER(LYROWS));
  
  UNPROTECT(1);
  return CPX;
}


SEXP R_PDMVSUM(SEXP X, SEXP LDIM, SEXP DESCX, SEXP Y, SEXP DESCY)
{
  const int m = INTEGER(LDIM)[0], n = INTEGER(LDIM)[1];
  
  SEXP CPX;
  PROTECT(CPX = allocMatrix(REALSXP, m, n));
  
  memcpy(REAL(CPX), REAL(X), m*n*sizeof(double));
  
  pdmvsum_(REAL(CPX), INTEGER(DESCX), REAL(Y), INTEGER(DESCY));
  
  UNPROTECT(1);
  return CPX;
}
