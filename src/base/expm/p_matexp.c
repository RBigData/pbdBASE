/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

// Copyright 2013, Schmidt

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define EXPSGN(x,pow) (x==0?(pow==0?1:0):(x>0?1:(pow%2==0?1:(-1))))
#define MIN(a,b) (a<b?a:b)

// ScaLAPACK functions
void pdgemm_(char *transa, char *transb, int *m, int *n, int *k, double *alpha, double *a, int *ia, int *ja, int *desca, double *b, int *ib, int *jb, int *descb, double *beta, double *c, int *ic, int *jc, int *descc);
void pdlacpy_(char *uplo, int *m, int *n, double *a, int *ia, int *ja, int *desca, double *b, int *ib, int *jb, int *descb);

// pbdBASE functions
pdims_(int *desc, int *ldm, int *blacs);
l2gpair_(int *i, int *j, int *gi, int *gj, int *desc, int *blacs);


// c = a * b for square matrices
static inline void p_matprod(double *a, int *desca, double *b, int *descb, double *c, int *descc)
{
  const char trans = 'N';
  const int ij = 1;
  const int n = desca[2];
  const double one = 1.0, zero = 0.0;
  
  pdgemm_(&trans, &trans, &n, &n, &n, &one, a, &ij, &ij, desca, b, &ij, &ij, descb, &zero, c, &ij, &ij, descc);
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
#if defined FIXTHISMESSLATERPLEASE
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


/*// Matrix exponentiation using Pade' approximations*/
/*// p==q==13*/
/*void p_matexp_pade(const unsigned int n, double *A, double *N, double *D)*/
/*{*/
/*  int i, j;*/
/*  int itmp;*/
/*  double tmp, tmpj;*/
/*  double *B, *C;*/
/*  */
/*  B = malloc(n*n*sizeof(double));*/
/*  C = malloc(n*n*sizeof(double));*/
/*  */
/*  // Initialize*/
/*  #pragma omp for simd*/
/*  {*/
/*    for (i=0; i<n*n; i++)*/
/*    {*/
/*      N[i] = 0.0;*/
/*      D[i] = 0.0;*/
/*    }*/
/*    */
/*    // Fill diagonal with 1*/
/*    j = 0;*/
/*    for (i=0; i<n*n; i+=n)*/
/*    {*/
/*      i += j;*/
/*      */
/*      N[i] = 1;*/
/*      D[i] = 1;*/
/*      */
/*      j = 1;*/
/*    }*/
/*  }*/
/*  */
/*  // Fill N and D*/
/*  for (i=1; i<=13; i++)*/
/*  {*/
/*    // C = A*B*/
/*    if (i>1)*/
/*      matprod(n, A, B, C);*/
/*    else*/
/*    {*/
/*      #pragma omp for simd*/
/*      {*/
/*        for (j=0; j<n*n; j++)*/
/*          C[j] = A[j];*/
/*      }*/
/*    }*/
/*    */
/*    #pragma omp for simd*/
/*    {*/
/*      // B = C*/
/*      for (j=0; j<n*n; j++)*/
/*        B[j] = C[j];*/
/*      */
/*      // N = pade_coef[i] * C*/
/*      // D = (-1)^j * pade_coef[i] * C*/
/*      tmp = dmat_pade_coefs[i];*/
/*      itmp = EXPSGN(-1, i);*/
/*      */
/*      if (itmp == 1)*/
/*      {*/
/*        for (j=0; j<n*n; j++)*/
/*        {*/
/*          tmpj = tmp * C[j];*/
/*          N[j] += tmpj;*/
/*          D[j] += tmpj;*/
/*        }*/
/*      }*/
/*      else*/
/*      {*/
/*        for (j=0; j<n*n; j++)*/
/*        {*/
/*          tmpj = tmp * C[j];*/
/*          N[j] += tmpj;*/
/*          D[j] -= tmpj;*/
/*        }*/
/*      }*/
/*    }*/
/*  }*/
/*  */
/*  free(B);*/
/*  free(C);*/
/*}*/



// Exponentiation by squaring
// P = A^b
void p_matpow_by_squaring(double *A, int *desca, int b, double *P)
{
  int n, m;
  int ldm[2];
  int blacs[5];
  int i, j;
  int itmp;
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

