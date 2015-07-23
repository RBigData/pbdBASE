/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

// Copyright 2012-2015, Schmidt

#include "pbdBASE.h"

void Cblacs_get(int ConTxt, int what, int *val);
void Cblacs_gridinit(int *ConTxt, char *order, int nprow, int npcol);
void Cblacs_gridinfo(int ConTxt, int *nprow, int *npcol, int *myrow, int *mycol);
void Cblacs_exit(int NotDone);


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



SEXP R_blacs_init(SEXP NPROW_in, SEXP NPCOL_in, SEXP ICTXT_in)
{
  R_INIT;
  SEXP NPROW, NPCOL, ICTXT, MYROW, MYCOL, RET, RET_NAMES;
  
  newRvec(NPROW, 1, "int");
  newRvec(NPCOL, 1, "int");
  newRvec(ICTXT, 1, "int");
  newRvec(MYROW, 1, "int");
  newRvec(MYCOL, 1, "int");
  
  INT(NPROW) = INT(NPROW_in);
  INT(NPCOL) = INT(NPCOL_in);
  INT(ICTXT) = INT(ICTXT_in);
  
  char order = 'R';
  
/*  sl_init_(INTP(ICTXT), INTP(NPROW), INTP(NPCOL));*/
  Cblacs_get(INT(ICTXT_in), 0, INTP(ICTXT));
  Cblacs_gridinit(INTP(ICTXT), &order, INT(NPROW), INT(NPCOL));
  Cblacs_gridinfo(INT(ICTXT), INTP(NPROW), INTP(NPCOL), INTP(MYROW), INTP(MYCOL));
  
  RET_NAMES = make_list_names(5, "NPROW", "NPCOL", "ICTXT", "MYROW", "MYCOL");
  RET = make_list(RET_NAMES, 5, NPROW, NPCOL, ICTXT, MYROW, MYCOL);
  
  R_END;
  return(RET);
}



SEXP R_blacs_exit(SEXP CONT)
{
  Cblacs_exit(INT(CONT));
  
  return R_NilValue;
}


