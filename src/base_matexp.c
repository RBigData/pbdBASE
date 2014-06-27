/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

// Copyright 2013, Schmidt

#include "base_global.h"
#include "base/expm/matexp.h"

#include <SEXPtools.h>


SEXP R_matpow_by_squaring(SEXP A, SEXP b)
{
  R_INIT;
  const int n = nrows(A);
  double *cpA;
  
  SEXP P;
  newRmat(P, n, n, "dbl");
  
  // A is modified
  //FIXME check return...
  cpA = malloc(n*n*sizeof(double));
  memcpy(cpA, REAL(A), n*n*sizeof(double));
  
  matpow_by_squaring(cpA, n, INT(b,0), REAL(P));
  
  free(cpA);
  
  R_END;
  return P;
}



SEXP R_matexp_pade(SEXP A, SEXP p)
{
  R_INIT;
  const int n = nrows(A);
  SEXP N, D;
  SEXP RET, RET_NAMES;
  
  
  // Allocate N and D
  newRmat(N, n, n, "dbl");
  newRmat(D, n, n, "dbl");
  
  
  // Compute N and D
  matexp_pade(n, INT(p,0), REAL(A), REAL(N), REAL(D));
  
  
  // Wrangle the return
  RET_NAMES = make_list_names(2, "N", "D");
  RET = make_list(RET_NAMES, 2, N, D);
  
  R_END;
  return RET;
}


