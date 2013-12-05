  /* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

// Copyright 2013, Schmidt

#include <stdio.h>
#include <stdlib.h>

#include "matexp.h"

void dgemm_(char *transa, char *transb, int *m, int *n, int *k, double *alpha, double *a, int *lda, double *b, int *ldb, double *beta, double *c, int *ldc);

// C = A * B for square matrices
static inline void matprod(const unsigned int n, double *a, double *b, double *c)
{
  char trans = 'N';
  const double one = 1.0, zero = 0.0;

  dgemm_(&trans, &trans, &n, &n, &n, &one, a, &n, b, &n, &zero, c, &n);
}

// Copy A ONTO B, i.e. B = A
static inline void matcopy(const unsigned int n, double *A, double *B)
{
  int i;

  #pragma omp for simd
  {
    for (i=0; i<n*n; i++)
      B[i] = A[i];
  }
}

// Zero matrix
static inline void matzero(const unsigned int n, double *a)
{
  int i;

  #pragma omp for simd
  {
    for (i=0; i<n*n; i++)
      a[i] = 0.0;
  }
}

// Identity matrix
static inline void mateye(const unsigned int n, double *a)
{
  int i, j;

  #pragma omp for simd
  {
    matzero(n, a);

    // Fill diagonal with 1
    j = 0;
    for (i=0; i<n*n; i+=n)
    {
      i += j;

      a[i] = 1.0;

      j = 1.0;
    }
  }
}


/* r_m(x) = p_m(x) / q_m(x), where
   p_m(x) = sum_{j=0}^m (2m-j)!m!/(2m)!/(m-j)!/j! * x^j

   and q_m(x) = p_m(-x)
*/


// Matrix exponentiation using Pade' approximations
// p==q==13
void matexp_pade(const unsigned int n, double *A, double *N, double *D)
{
  int i, j;
  int sign;
  double tmp, tmpj;
  double *B, *C;

  B = calloc(n*n, sizeof(double));
  C = calloc(n*n, sizeof(double));

  // Initialize
  #pragma omp for simd
  {
    // Fill diagonal with 1
    j = 0;
    for (i=0; i<n*n; i+=n)
    {
      i += j;
      
      N[i] = 1;
      D[i] = 1;
      
      j = 1;
    }
  }
  
  // Fill N and D
  for (i=1; i<=13; i++)
  {
    // C = A*B
    if (i == 1)
      matcopy(n, A, C);
    else
      matprod(n, A, B, C);
    
    
    // N = pade_coef[i] * C
    // D = (-1)^j * pade_coef[i] * C
    tmp = matexp_pade_coefs[i];
    sign = SGNEXP(-1, i);
    
    if (sign == 1)
    {
      for (j=0; j<n*n; j++)
      {
        B[j] = C[j];
        tmpj = tmp * C[j];
        N[j] += tmpj;
        D[j] += tmpj;
      }
    }
    else
    {
      for (j=0; j<n*n; j++)
      {
        B[j] = C[j];
        tmpj = tmp * C[j];
        N[j] += tmpj;
        D[j] -= tmpj;
      }
    }
  }
  
  free(B);
  free(C);
}



// Exponentiation by squaring
// P = A^b
void matpow_by_squaring(double *A, int n, int b, double *P)
{
  int i, j;
  double tmp, tmpj;
  double *TMP;
  
  mateye(n, P);
  
  
  // Trivial cases
  if (b == 0)
    return;
  
  if (b == 1)
  {
    matcopy(n, A, P);
    return;
  }


  // General case
  TMP = malloc(n*n*sizeof(double));
  
  while (b)
  {
    if (b&1)
    {
      matprod(n, P, A, TMP);
      matcopy(n, TMP, P);
    }
    
    b >>=1;
    matprod(n, A, A, TMP);
    matcopy(n, TMP, A);
  }
  
  free(TMP);
}

