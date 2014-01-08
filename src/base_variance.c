/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

// Copyright 2013, Schmidt

#include <R.h>
#include <Rinternals.h>
#include "base_global.h"

SEXP R_PDCLVAR(SEXP X, SEXP DESCX, SEXP LSD)
{
  SEXP VAR;
  PROTECT(VAR = allocVector(REALSXP, INTEGER(LSD)[0]));
  
  pdclvar_(REAL(X), INTEGER(DESCX), REAL(VAR));
  
  UNPROTECT(1);
  return(VAR);
} 
