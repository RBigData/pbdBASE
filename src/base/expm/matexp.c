  /* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

// Copyright 2013, Schmidt

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#ifdef _OPENMP
  #include <omp.h>
  #if _OPENMP >= 201307
    #define _OPENMP_SUPPORT_SIMD
  #endif
#endif

#include "matexp.h"


void dgemm_(char *transa, char *transb, int *m, int *n, int *k, double *alpha, double *a, int *lda, double *b, int *ldb, double *beta, double *c, int *ldc);
void dlacpy_(char *uplo, int *m, int *n, double *a, int *lda, double *b, int *ldb);

// C = A * B for square matrices
static inline void matprod(int n, double *a, double *b, double *c)
{
  char trans = 'N';
  double one = 1.0, zero = 0.0;
  
  dgemm_(&trans, &trans, &n, &n, &n, &one, a, &n, b, &n, &zero, c, &n);
}

// Copy A ONTO B, i.e. B = A
static inline void matcopy(int n, double *A, double *B)
{
  char uplo = 'A';
  
  dlacpy_(&uplo, &n, &n, A, &n, B, &n);
}

// Zero matrix
static inline void matzero(const unsigned int n, double *a)
{
  int i;
  
  #if defined(_OPENMP_SUPPORT_SIMD)
  #pragma omp for simd
  #endif
    for (i=0; i<n*n; i++)
      a[i] = 0.0;
}

// Identity matrix
static inline void mateye(const unsigned int n, double *a)
{
  int i;
  
  matzero(n, a);
  
  // Fill diagonal with 1
  i = 0;
  while (i < n*n)
  {
    a[i] = 1.0;
    
    i += n+1;
  }
}




// Exponentiation by squaring --- A is modified
// P = A^b
void matpow_by_squaring(double *A, int n, int b, double *P)
{
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





/* r_m(x) = p_m(x) / q_m(x), where
   p_m(x) = sum_{j=0}^m (2m-j)!m!/(2m)!/(m-j)!/j! * x^j

   and q_m(x) = p_m(-x)
*/


// Matrix exponentiation using Pade' approximations
// p==q==13

void matexp_pade_fillmats(const unsigned int m, const unsigned int n, const unsigned int i, double *N, double *D, double *B, double *C)
{
  int j;
  const double tmp = matexp_pade_coefs[i];
  double tmpj;
  
  
  if (SGNEXP(-1, i) == 1)
  {
    for (j=0; j<m*n; j++)
    {
      // B = C
      tmpj = C[j];
      B[j] = tmpj;
      
      tmpj *= tmp;
      // N = pade_coef[i] * C
      N[j] += tmpj;
      // D = (-1)^j * pade_coef[i] * C
      D[j] += tmpj;
    }
  }
  else
  {
    for (j=0; j<m*n; j++)
    {
      // B = C
      tmpj = C[j];
      B[j] = tmpj;
      
      tmpj *= tmp;
      // N = pade_coef[i] * C
      N[j] += tmpj;
      // D = (-1)^j * pade_coef[i] * C
      D[j] -= tmpj;
    }
  }
}

void matexp_pade(const unsigned int n, double *A, double *N, double *D)
{
  int i;
  double *B, *C;

  // Power of A
  B = calloc(n*n, sizeof(double));
  // Temporary storage for matrix multiplication
  C = malloc(n*n * sizeof(double));

  assert(B != NULL);
  assert(C != NULL);

  matcopy(n, A, C);

  for (i=0; i<n*n; i++)
  {
    N[i] = 0.0;
    D[i] = 0.0;
  }

  // Initialize N and D
  // Fill diagonal with 1
  i = 0;
  while (i < n*n)
  {
    N[i] = 1.0;
    D[i] = 1.0;
    
    i += n+1;
  }
  
  
  // Fill N and D
  for (i=1; i<=PADE_PQ; i++)
  {
    // C = A*B
    if (i > 1)
      matprod(n, A, B, C);
      
    // Update matrices
    matexp_pade_fillmats(n, n, i, N, D, B, C);
  }
  
  free(B);
  free(C);
}


