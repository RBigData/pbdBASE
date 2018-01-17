/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

// Copyright 2013-2016 Schmidt

#include "pbdBASE.h"
#include "base/expm/matexp.h"


SEXP R_p_matpow_by_squaring(SEXP A, SEXP desca, SEXP b)
{
  R_INIT;
  const int m = nrows(A), n = ncols(A);
  double *cpA;
  
  SEXP P;
  newRmat(P, nrows(A), ncols(A), "dbl");
  
  
  // Why did I make a copy ... ? // Oh now I remember
  //FIXME check returns...
  cpA = malloc(m*n * sizeof(double));
  memcpy(cpA, REAL(A), m*n*sizeof(double));
  
  p_matpow_by_squaring(cpA, INTEGER(desca), INT(b, 0), REAL(P));
  
  free(cpA);
  
  R_END;
  return(P);
}




SEXP R_p_matexp_pade(SEXP A, SEXP desca, SEXP p)
{
  R_INIT;
  int m, n;
  SEXP N, D;
  SEXP RET, RET_NAMES;
  
  m = nrows(A);
  n = ncols(A);
  
  
  // Allocate N and D
  newRmat(N, m, n, "dbl");
  newRmat(D, m, n, "dbl");
  
  
  // Compute N and D
  p_matexp_pade(DBLP(A), INTP(desca), INT(p, 0), DBLP(N), DBLP(D));
  
  
  // Wrangle the return
  make_list_names(RET_NAMES, 2, "N", "D");
  make_list(RET, RET_NAMES, 2, N, D);
  
  R_END;
  return RET;
}
