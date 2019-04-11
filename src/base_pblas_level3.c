/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

// Copyright 2013, Schmidt and Chen

#include <RNACI.h>

#include "pblas.h"
#include "pbdBASE.h"


// Transpose
SEXP R_PDTRAN(SEXP M, SEXP N, SEXP A, SEXP DESCA, SEXP CLDIM, SEXP DESCC)
{
  R_INIT;
  int IJ = 1;
  double one = 1.0;
  double zero = 0.0;
  
  SEXP C;
  newRmat(C, INT(CLDIM, 0), INT(CLDIM, 1), "dbl");
  
  pdtran_(INTP(M), INTP(N), &one, 
          DBLP(A), &IJ, &IJ, INTP(DESCA), &zero, 
          DBLP(C), &IJ, &IJ, INTP(DESCC));
  
  
  R_END;
  return C;
}



// Mat-Mat-Mult
SEXP R_PDGEMM(SEXP TRANSA, SEXP TRANSB, SEXP M, SEXP N, SEXP K,
    SEXP A, SEXP DESCA, SEXP B, SEXP DESCB, SEXP CLDIM, SEXP DESCC)
{
  R_INIT;
  double alpha = 1.0;
  double beta = 0.0;
  int IJ = 1;
  
  SEXP C;
  newRmat(C, INT(CLDIM, 0), INT(CLDIM, 1), "dbl");
  
  
  pdgemm_(STR(TRANSA, 0), STR(TRANSB, 0),
      INTP(M), INTP(N), INTP(K), &alpha, 
      DBLP(A), &IJ, &IJ, INTP(DESCA),
      DBLP(B), &IJ, &IJ, INTP(DESCB), &beta,
      DBLP(C), &IJ, &IJ, INTP(DESCC));
  
  
  R_END;
  return C;
}
