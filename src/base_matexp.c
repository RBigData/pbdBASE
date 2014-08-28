/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

// Copyright 2013-2014, Schmidt


#include "pbdBASE.h"
#include "base/expm/matexp.h"


SEXP R_matexp(SEXP A, SEXP p)
{
  R_INIT;
  const int n = nrows(A);
  int i;
  double *A_cp;
  SEXP R;
  
  newRmat(R, n, n, "dbl");
  
  A_cp = malloc(n*n*sizeof(A_cp));
  
  for (i=0; i<n*n; i++)
    A_cp[i] = REAL(A)[i];
  
  
  double t = 1;
  matexp(n, INT(p), t, A_cp, REAL(R));
  
  free(A_cp);
  
  R_END;
  return R;
}


