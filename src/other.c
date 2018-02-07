/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

// Copyright 2013-2014, Schmidt and Chen

#include "pbdBASE.h"


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
