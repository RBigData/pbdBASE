// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// Copyright 2013-2014, 2016 Schmidt


// sweep array out of distributed matrix
// inputs/outputs
//   x = submatrix of data which should globally be "swept"
// inputs
//   ix/jx = 
//   descx = descriptor array for x
//   vec = vector to "sweep" through x
//   lvec = length of vec
//   margin = 1 for row sweeping, 2 for column sweeping
//   fun = char with 4 possibilities, describing the type of sweep to perform:
//     "+", "-", "*", "/"

#include "utils.h"


#define UNUSED(x) (void)(x)

#define ROWS 1
#define COLS 2

static inline int ind(const int i, const int m)
{
  int ind = i%m;
  if (ind == 0)
    ind = m;
  
  return ind;
}


static inline int indxl2g(const int indxloc, const int nb, const int iproc, const int isrcproc, const int nprocs)
{
  return nprocs*nb*((indxloc - 1)/nb) + ((indxloc - 1) % nb) + ((nprocs + iproc - isrcproc) % nprocs)*nb + 1;
}


static inline void l2gpair(const int i, const int j, int* restrict gi, int* restrict gj, const int *restrict desc, const int *restrict blacs)
{
  *gi = indxl2g(i, desc[5-1], blacs[4-1], 0, blacs[2-1]);
  *gj = indxl2g(j, desc[6-1], blacs[5-1], 0, blacs[3-1]);
}


#define SWEEP(OP) \
  { \
    if (margin == ROWS) \
    { \
      k = descx[2]; \
      for (int j=1; j<=n; j++) \
      { \
        for (int i=1; i<=m; i++) \
        { \
          l2gpair(i, j, &gi, &gj, descx, blacs); \
          pos = ind(gi + k*(gj-1), lvec) -1; \
          x[i-1 + m*(j-1)] = x[i-1 + m*(j-1)] OP vec[pos]; \
        } \
      } \
    } \
    else if (margin == COLS) \
    { \
      k = descx[3]; \
      for (int j=1; j<=n; j++) \
      { \
        for (int i=1; i<=m; i++) \
        { \
          l2gpair(i, j, &gi, &gj, descx, blacs); \
          pos = ind(gj + k*(gi-1), lvec) - 1; \
          x[i-1 + m*(j-1)] = x[i-1 + m*(j-1)] OP vec[pos]; \
        } \
      } \
    } \
  }



void pdsweep(double *restrict x, const int ix, const int jx, int *restrict descx, double *restrict vec, const int lvec, const int margin, const char fun)
{
  UNUSED(ix);
  UNUSED(jx);
  
  int k, m, n, pos, gi, gj;
  int ldm[2], blacs[5];
  
  // get local and proc grid info
  pdims_(descx, ldm, blacs);
  
  m = ldm[0];
  n = ldm[1];
  
  // only do work if we own any local pieces
  if (m == 0 || n == 0)
    return;
  
  // addition
  if (fun == '+')
    SWEEP(+)
  else if (fun == '-')
    SWEEP(-)
  else if (fun == '*')
    SWEEP(*)
  else if (fun == '/')
    SWEEP(/)
  
}
