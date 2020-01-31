/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

// Copyright 2015, Schmidt

#include "base/utils/utils.h"
#include "scalapack.h"

// R.h and Rinternals.h needs to be included after Rconfig.h
#include "pbdBASE.h"
#include <RNACI.h>


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



SEXP R_NUMROC(SEXP N, SEXP NB, SEXP IPROC, SEXP NPROCS)
{
  R_INIT;
  SEXP NUM;
  newRvec(NUM, 1, "int");
  
  numrocwrap_(INTP(N), INTP(NB), INTP(IPROC), INTP(NPROCS), INTP(NUM));
  
  R_END;
  return NUM;
}



static void l2g_coord(int* ret, int i, int j, int* bldim, int* procs, int* myproc)
{
  const int nprocs = procs[0] * procs[1];
  
  ret[0] = nprocs*bldim[0] * (i-1)/bldim[0] + (i-1)%bldim[0] + ((nprocs+myproc[0])%nprocs)*bldim[0] + 1;
  ret[1] = nprocs*bldim[1] * (j-1)/bldim[1] + (j-1)%bldim[1] + ((nprocs+myproc[1])%nprocs)*bldim[1] + 1;
}

SEXP l2g_coords(SEXP ind, SEXP bldim, SEXP procs, SEXP myproc)
{
  SEXP ret;
  PROTECT(ret = allocVector(INTSXP, 6));
  
  l2g_coord(INTEGER(ret), INTEGER(ind)[0], INTEGER(ind)[1], INTEGER(bldim), INTEGER(procs), INTEGER(myproc));
  
  UNPROTECT(1);
  return ret;
}



static void g2l_coord(int* ret, int i, int j, int* bldim, int* procs, int* src)
{
  // matrix block position
  ret[0] = i / (procs[0] * bldim[0]);
  ret[1] = j / (procs[1] * bldim[1]);
  
  // process grid block
  ret[2] = (src[0] + i / bldim[0]) % procs[0];
  ret[3] = (src[1] + j / bldim[1]) % procs[1];
  
  // local coordinates
  ret[4] = i % bldim[0] + bldim[0] * ret[0];
  ret[5] = j % bldim[1] + bldim[1] * ret[1];
}

SEXP g2l_coords(SEXP ind, SEXP bldim, SEXP procs, SEXP src)
{
  SEXP ret;
  PROTECT(ret = allocVector(INTSXP, 6));
  
  g2l_coord(INTEGER(ret), INTEGER(ind)[0], INTEGER(ind)[1], INTEGER(bldim), INTEGER(procs), INTEGER(src));
  
  UNPROTECT(1);
  return ret;
}
