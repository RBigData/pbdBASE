/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

// Copyright 2013, Schmidt

#include "base/stats/stats.h"

// R.h and Rinternals.h needs to be included after Rconfig.h
#include "pbdBASE.h"
#include <RNACI.h>


SEXP R_PDCLVAR(SEXP X, SEXP DESCX, SEXP LSD)
{
  R_INIT;
  SEXP VAR;
  
  newRvec(VAR, INT(LSD, 0), "dbl");
  
  pdclvar_(REAL(X), INTEGER(DESCX), REAL(VAR));
  
  R_END;
  return VAR;
} 
