/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

// Copyright 2013, Schmidt

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "matexp.h"



// c = a * b
static inline void p_matprod(double *a, int *desca, double *b, int *descb, double *c, int *descc)
{
  const char trans = 'N';
  const int ij = 1;
  const double one = 1.0, zero = 0.0;
  int m, n, k;
  
  
  m = desca[2];
  n = descb[3];
  k = desca[3];
  
  
  pdgemm_(&trans, &trans, &m, &n, &k, &one, a, &ij, &ij, desca, b, &ij, &ij, descb, &zero, c, &ij, &ij, descc);
}


// Copy a ONTO b, i.e. b = a
static inline void p_matcopy(double *a, int *desca, double *b, int *descb)
{
  const char uplo = 'A';
  const int n = desca[2];
  const int ij = 1;
  
  pdlacpy_(&uplo, &n, &n, a, &ij, &ij, desca, b, &ij, &ij, descb);
}


// Identity matrix
/*static inline */
void p_mateye(double *a, int *desca)
{
  const int ij = 1;
  double diag = 1.0;
  int ldiag = 1;
  pddiagmk_(a, &ij, &ij, desca, &diag, &ldiag);
}
// Fix this later, this is so stupid
#if 0
{
  int i, j, gi, gj, ti, tj;
  int mb_a, nb_a, minb_a;
  int ldm[2];
  int blacs[5];
  int m, n, gm;
  
  // Get local dim and context grid info
  pdims_(desca, ldm, blacs);
  
  m = ldm[0];
  n = ldm[1];
  
  gm = desca[2];
  
  mb_a = desca[4];
  nb_a = desca[5];
  
  mb_a = MIN(mb_a, gm);
  nb_a = MIN(nb_a, gm);
  
  // Initialize
  for (i=0; i<m*n; i++)
    a[i] = 0.0;
  
  // Fill diagonal with 1's
  for (j=0; j<n; j++)//+=nb_a)
  {
    for (i=0; i<m; i++)//+=mb_a)
    {
      l2gpair_(&i, &j, &gi, &gj, desca, blacs);
      
/*      if (m-i+1 > 0)*/
        mb_a = MIN(mb_a, m-i+1);
/*      if (n-j+1 > 0)*/
        nb_a = MIN(nb_a, n-j+1);
      
      
      
      minb_a = mb_a<=nb_a?mb_a:nb_a;
      
      if (gi == gj)
        a[i + m*j] = 1.0;
      // only try if a diagonal entry is in this sub-block
/*      if (abs(gi-gj) <= minb_a)*/
/*      {*/
/*        for (tj=0; tj<nb_a; tj++)*/
/*        {*/
/*          for (ti=0; ti<mb_a; ti++)*/
/*          {*/
/*            if (gi+ti == gj+tj)*/
/*            {*/
/*              a[i+ti + m*(j+tj)] = 1.0;*/
/*            }*/
/*          }*/
/*        }*/
/*      }*/
    }
  }
  
}
#endif




// Exponentiation by squaring
// P = A^b
void p_matpow_by_squaring(double *A, int *desca, int b, double *P)
{
  int n, m;
  int ldm[2];
  int blacs[5];
  int i, j;
  double tmp, tmpj;
  double *TMP;
  
  p_mateye(P, desca);
  
  
  // Trivial cases
  if (b == 0)
    return;
  
  if (b == 1)
  {
    p_matcopy(A, desca, P, desca);
    return;
  }
  
  
  // Get local dim and context grid info
  pdims_(desca, ldm, blacs);
  
  m = ldm[0];
  n = ldm[1];
  
  
  // General case
  TMP = malloc(m*n*sizeof(double));
  
  while (b)
  {
    if (b&1)
    {
      p_matprod(P, desca, A, desca, TMP, desca);
      p_matcopy(TMP, desca, P, desca);
    }
    
    b >>=1;
    p_matprod(A, desca, A, desca, TMP, desca);
    p_matcopy(TMP, desca, A, desca);
  }
  
  free(TMP);
}




// Matrix exponentiation using Pade' approximations
// p==q==13
void matexp_pade_fillmats(const unsigned int m, const unsigned int n, const unsigned int i, double *N, double *D, double *B, double *C);

void p_matexp_pade(double *A, int *desca, double *N, double *D)
{
  int m, n;
  int i;
  int ldm[2], blacs[5];
  double *B, *C;
  
  // Get local dim and context grid info
  pdims_(desca, ldm, blacs);
  
  m = ldm[0];
  n = ldm[1];
  
  m = m?m:1;
  n = n?n:1;
  
  // Power of A
  B = calloc(m*n, sizeof(double));
  // Temporary storage for matrix multiplication
  C = malloc(m*n*sizeof(double));
  
  assert(B != NULL);
  assert(C != NULL);
  
  p_matcopy(A, desca, C, desca);
  
  // Initialize
  p_mateye(D, desca);
  memcpy(N, D, m*n*sizeof(double));
  
  // Fill N and D
  for (i=1; i<=13; i++)
  {
    // C = A*B
    if (i > 1)
      p_matprod(A, desca, B, desca, C, desca);
    
    // Update matrices
    matexp_pade_fillmats(m, n, i, N, D, B, C);
  }
  
  free(B);
  free(C);
}


