/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

// Copyright 2013, Schmidt


#include "base/linalg/linalg.h"

// R.h and Rinternals.h needs to be included after Rconfig.h
#include "pbdBASE.h"
#include <RNACI.h>


SEXP R_PDCROSSPROD(SEXP UPLO, SEXP TRANS, SEXP A, SEXP DESCA, SEXP CLDIM, SEXP DESCC)
{
  R_INIT;
  double alpha = 1.0;
  int IJ = 1;
  
  SEXP C;
  newRmat(C, INT(CLDIM, 0), INT(CLDIM, 1), "dbl");
  
#ifdef FC_LEN_T
  pdcrossprod_(STR(UPLO, 0), STR(TRANS, 0), &alpha, 
    DBLP(A), &IJ, &IJ, INTP(DESCA), 
    DBLP(C), &IJ, &IJ, INTP(DESCC),
    (FC_LEN_T) strlen(STR(UPLO, 0)), (FC_LEN_T) strlen(STR(TRANS, 0)));
#else
  pdcrossprod_(STR(UPLO, 0), STR(TRANS, 0), &alpha, 
    DBLP(A), &IJ, &IJ, INTP(DESCA), 
    DBLP(C), &IJ, &IJ, INTP(DESCC));
#endif
  
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
  
#ifdef FC_LEN_T
  pdchtri_(CHARPT(UPLO, 0), A_cp, &IJ, &IJ, INTEGER(DESCA), 
    REAL(C), &IJ, &IJ,
    INTEGER(DESCC), &info,
    (FC_LEN_T) strlen(CHARPT(UPLO, 0)));
#else
  pdchtri_(CHARPT(UPLO, 0), A_cp, &IJ, &IJ, INTEGER(DESCA), 
    REAL(C), &IJ, &IJ,
    INTEGER(DESCC), &info);
#endif
  
  if (info != 0)
  {
    //FIXME replace with appropriate COMM_WARN
    Rprintf("INFO = %d\n", info);
  }
  
  UNPROTECT(1);
  return(C);
}
