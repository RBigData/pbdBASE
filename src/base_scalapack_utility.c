/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

// Copyright 2013, Schmidt

#include "base_global.h"
#include <SEXPtools.h>


SEXP R_PDLAPRNT(SEXP M, SEXP N, SEXP A, SEXP DESCA, SEXP CMATNM, SEXP NOUT)
{
  double work[INTEGER(DESCA)[8]];
  const int IJ = 1;
  const int SRC = 0;
  
  bprnt_(INTEGER(M), INTEGER(N), REAL(A), &IJ, &IJ,
         INTEGER(DESCA), &SRC, &SRC, CHARPT(CMATNM, 0),
         INTEGER(NOUT), &work);
  
  return RNULL;
}



SEXP R_PDGEMR2D(SEXP M, SEXP N, SEXP X, SEXP DESCX, SEXP CLDIM, SEXP DESCB, SEXP CTXT)
{
  R_INIT;
  const int IJ = 1;
  SEXP B;
  
  newRmat(B, INT(CLDIM, 0), INT(CLDIM, 1), "dbl");
  
  
  Cpdgemr2d(INTEGER(M)[0], INTEGER(N)[0],
      REAL(X), IJ, IJ, INTEGER(DESCX),
      REAL(B), IJ, IJ, INTEGER(DESCB), INTEGER(CTXT)[0]);
  
  
  R_END;
  return B;
}



// next best divisor function
SEXP R_nbd(SEXP N, SEXP D)
{
  R_INIT;
  int i, test;
  const int n = INT(N, 0);
  const int d = INT(D, 0);
  
  SEXP RET;
  newRvec(RET, 1, "int");
  INT(RET, 0) = d;
  
  for (i=INT(RET, 0); i<=n; i++)
  {
    test = n % i;
    if (test == 0){
      INT(RET, 0) = i;
      break;
    }
  }
  
  
  R_END;
  return RET;
}



