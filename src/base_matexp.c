/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

// Copyright 2013-2015, Schmidt


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
  
  A_cp = (double *) R_alloc(n*n, sizeof(A_cp));
  
  for (i=0; i<n*n; i++)
    A_cp[i] = REAL(A)[i];
  
  
  matexp(n, INT(p), A_cp, REAL(R));
  
  R_END;
  return R;
}


