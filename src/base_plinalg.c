/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

// Copyright 2013, Schmidt

#include "pbdBASE.h"


SEXP R_PDCROSSPROD(SEXP UPLO, SEXP TRANS, SEXP A, SEXP DESCA, SEXP CLDIM, SEXP DESCC)
{
  R_INIT;
  double alpha = 1.0;
  int IJ = 1;
  
  SEXP C;
  newRmat(C, INT(CLDIM, 0), INT(CLDIM, 1), "dbl");
  
  pdcrossprod_(STR(UPLO, 0), STR(TRANS, 0), &alpha, 
    DBLP(A), &IJ, &IJ, INTP(DESCA), 
    DBLP(C), &IJ, &IJ, INTP(DESCC));
  
  R_END;
  return C;
}


SEXP R_PDCHTRI(SEXP UPLO, SEXP A, SEXP ALDIM, SEXP DESCA, SEXP CLDIM, SEXP DESCC)
{
  int IJ = 1;
  const int m = INTEGER(ALDIM)[0], n = INTEGER(ALDIM)[1];
  double *A_cp;
  int info = 0;
  
  SEXP C;
  PROTECT(C = allocMatrix(REALSXP, INTEGER(CLDIM)[0], INTEGER(CLDIM)[1]));
  
  A_cp = (double *) R_alloc(m*n, sizeof(double));
  memcpy(A_cp, REAL(A), m*n*sizeof(double));
  
  pdchtri_(CHARPT(UPLO, 0), A_cp, &IJ, &IJ, INTEGER(DESCA), 
    REAL(C), &IJ, &IJ,
    INTEGER(DESCC), &info);
  
  if (info != 0)
  {
    //FIXME replace with appropriate COMM_WARN
    Rprintf("INFO = %d\n", info);
  }
  
  UNPROTECT(1);
  return(C);
}

