/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

// Copyright 2012-2014, Schmidt

#include "base_global.h"
#include <SEXPtools.h>


SEXP R_optimal_grid(SEXP NPROCS)
{
  R_INIT;
  SEXP NPROW, NPCOL, RET, RET_NAMES;
  
  newRvec(NPROW, 1, "int", TRUE);
  newRvec(NPCOL, 1, "int", TRUE);
  
  optimalgrid_(INTP(NPROCS), INTP(NPROW), INTP(NPCOL));
  
  RET_NAMES = make_list_names(2, "nprow", "npcol");
  RET = make_list(RET_NAMES, 2, NPROW, NPCOL);
  
  R_END;
  return RET;
}



SEXP R_blacs_init(SEXP NPROW, SEXP NPCOL, SEXP ICTXT)
{
  R_INIT;
  SEXP MYROW, MYCOL, RET, RET_NAMES;
  
  newRvec(MYROW, 1, "int");
  newRvec(MYCOL, 1, "int");
  
  
  sl_init_(INTEGER(ICTXT), INTEGER(NPROW), INTEGER(NPCOL));
  blacs_gridinfo_(INTEGER(ICTXT), INTEGER(NPROW), INTEGER(NPCOL), 
      INTEGER(MYROW), INTEGER(MYCOL));
  
  
  RET_NAMES = make_list_names(5, "NPROW", "NPCOL", "ICTXT", "MYROW", "MYCOL");
  RET = make_list(RET_NAMES, 5, NPROW, NPCOL, ICTXT, MYROW, MYCOL);
  
  R_END;
  return(RET);
}


