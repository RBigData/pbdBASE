/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

// Copyright 2016, Schmidt


#include "blacs.h"

// R.h and Rinternals.h needs to be included after Rconfig.h
#include "pbdBASE.h"
#include <RNACI.h>


void dmat_ldimget(int *desc, int* nrows, int* ncols);
int dmat_as_ddmatrix(int *desc, double *A_global, double *A_local);

SEXP R_redist(SEXP desc, SEXP A)
{
  R_INIT;
  SEXP ret;
  int nrows, ncols;
  
  dmat_ldimget(INTP(desc), &nrows, &ncols);
  newRmat(ret, nrows, ncols, "dbl");
  
  dmat_as_ddmatrix(INTP(desc), DBLP(A), DBLP(ret));
  
  R_END;
  return ret;
}
