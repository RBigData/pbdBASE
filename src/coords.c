/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

// Copyright 2015, Schmidt

#include "pbdBASE.h"


// C translation of indxg2l
static inline int indxg2l(const int INDXGLOB, const int NB, const int NPROCS)
{
  return NB * ((INDXGLOB - 1)/(NB*NPROCS)) + ((INDXGLOB - 1) % NB) + 1;
}

// C translation of indxg2p
static inline int indxg2p(const int INDXGLOB, const int NB, const int NPROCS)
{
  return ((INDXGLOB - 1) / NB) % NPROCS;
}




static inline int matcoord(const int *restrict dim, const int *restrict bldim, 
  const int gi, const int gj, int *restrict i, int *restrict j,
  const int nprow, const int npcol, const int myrow, const int mycol)
{
  if (gi < 1 || gj < 1 || gi > dim[0] || gj > dim[1])
    return NA_INTEGER;
  
  *i = indxg2l(gi, bldim[0], nprow);
  *j = indxg2l(gj, bldim[1], npcol);
  //indxg2lsub_(&gi, bldim, &nprow, i);
  //indxg2lsub_(&gj, bldim+1, &npcol, j);
  
  const int pi = indxg2p(gi, bldim[0], nprow);
  const int pj = indxg2p(gj, bldim[1], npcol);
  
  if (pi == myrow && pj == mycol)
    return 0;
  else
    return -1;
}



static inline int infog1l(const int global, const int bldim, const int np, const int myp)
{
  int block = (global-1)/bldim;
  int src = block % np;
  int local = (block / np + 1) * bldim + 1;
  
  if (((myp + np) % np) >= block)
  {
    if (myp == src)
      local += ((global-1) % bldim);
    
    local -= bldim;
  }
  
  if (myp == src)
    return local;
  else
    return -1;
}



#define intvecelt(vec,i) INT(VECTOR_ELT(vec,i))
SEXP R_g2lcoord(SEXP dim, SEXP bldim, SEXP gi, SEXP gj, SEXP gridinfo)
{
  R_INIT;
  SEXP ret;
  
  int i = NA_INTEGER;
  int j = NA_INTEGER;
  
  // gridinfo = list(nprow, npcol, ictxt, myrow, mycol)
  int check = matcoord(INTP(dim), INTP(bldim), INT(gi), INT(gj), &i, &j,
    intvecelt(gridinfo, 0), intvecelt(gridinfo, 1), intvecelt(gridinfo, 3), intvecelt(gridinfo, 4));
  
  //if (i == -1 || j == -1)
  if (check == NA_INTEGER || check == -1)
  {
    newRvec(ret, 1, "int");
    INT(ret) = check;
  }
  else
  {
    newRvec(ret, 2, "int");
    INT(ret, 0) = i;
    INT(ret, 1) = j;
  }
  
  R_END;
  return ret;
}

