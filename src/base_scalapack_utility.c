/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

// Copyright 2013, Schmidt

#include "pbdBASE.h"


SEXP R_PDLAPRNT(SEXP M, SEXP N, SEXP A, SEXP DESCA, SEXP CMATNM, SEXP NOUT)
{
  double work[INTEGER(DESCA)[8]];
  int IJ = 1;
  int SRC = 0;
  
  bprnt_(INTEGER(M), INTEGER(N), REAL(A), &IJ, &IJ,
         INTEGER(DESCA), &SRC, &SRC, CHARPT(CMATNM, 0),
         INTEGER(NOUT), work);
  
  return RNULL;
}



// redistributors
SEXP R_PIGEMR2D(SEXP M, SEXP N, SEXP X, SEXP DESCX, SEXP CLDIM, SEXP DESCB, SEXP CTXT)
{
  R_INIT;
  int IJ = 1;
  SEXP C;
  
  newRmat(C, INT(CLDIM, 0), INT(CLDIM, 1), "dbl");
  
  Cpigemr2d(INT(M), INT(N),
      INTEGER(X), IJ, IJ, INTEGER(DESCX),
      INTEGER(C), IJ, IJ, INTEGER(DESCB), INT(CTXT));
  
  R_END;
  return C;
}



SEXP R_PDGEMR2D(SEXP M, SEXP N, SEXP X, SEXP DESCX, SEXP CLDIM, SEXP DESCB, SEXP CTXT)
{
  R_INIT;
  int IJ = 1;
  SEXP C;
  
  newRmat(C, INT(CLDIM, 0), INT(CLDIM, 1), "dbl");
  
  Cpdgemr2d(INT(M), INT(N),
      REAL(X), IJ, IJ, INTEGER(DESCX),
      REAL(C), IJ, IJ, INTEGER(DESCB), INT(CTXT));
  
  R_END;
  return C;
}



// next best divisor function
SEXP R_nbd(SEXP N, SEXP D)
{
  R_INIT;
  int i, test;
  const int n = INT(N);
  const int d = INT(D);
  
  SEXP RET;
  newRvec(RET, 1, "int");
  INT(RET) = d;
  
  for (i=INT(RET, 0); i<=n; i++)
  {
    test = n % i;
    if (test == 0){
      INT(RET) = i;
      break;
    }
  }
  
  
  R_END;
  return RET;
}
