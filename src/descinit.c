/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

// Copyright 2015, Schmidt

#include "pbdBASE.h"


SEXP R_descinit(SEXP DIM, SEXP BLDIM, SEXP ICTXT, SEXP LLD)
{
  R_INIT;
  int row_col_src = 0;
  int info = 0;
  SEXP desc;
  newRvec(desc, 9, "int");
  
  descinit_(INTP(desc), INTP(DIM), INTP(DIM)+1, INTP(BLDIM), INTP(BLDIM)+1, 
    &row_col_src, &row_col_src, INTP(ICTXT), INTP(LLD), &info);
  
  R_END;
  return desc;
}

