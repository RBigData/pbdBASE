/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

// Copyright 2013-2014, Drew Schmidt

// Matrix exponentiation algorithm from:
// "New Scaling and Squaring Algorithm for the Matrix Exponential",
// Awad H. Al-Mohy and Nicholas J. Higham, August 2009


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


// -------------------------------------------------------- 
// Utilities
// -------------------------------------------------------- 

// C = A * B for square matrices
static void matprod(int n, double *a, double *b, double *c)
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



// Identity matrix
static inline void mateye(const unsigned int n, double *a)
{
  int i;
  
  for (i=0; i<n*n; i++)
    a[i] = 0.0;
  
  // Fill diagonal with 1
  i = 0;
  while (i < n*n)
  {
    a[i] = 1.0;
    
    i += n+1;
  }
}



// 1-norm for a square matrix
static double matnorm_1(const double *x, const int n)
{
  int i, j;
  double norm = 0;
  double tmp;
  
  // max(colSums(abs(x))) 
  for (j=0; j<n; j++)
  {
    tmp = 0;
    
    for (i=0; i<n; i++)
      tmp += fabs(x[i*n + j]);
    
    if (tmp > norm)
      norm = tmp;
  }
  
  return norm;
}



#define NTHETA 5

static int matexp_scale_factor(const double *x, const int n)
{
  int i;
  const double theta[] = {1.5e-2, 2.5e-1, 9.5e-1, 2.1e0, 5.4e0};
  
  const double x_1 = matnorm_1(x, n);
  
  for (i=0; i<NTHETA; i++)
  {
    if (x_1 <= theta[i])
      return 0;
  }
  
  i = (int) ceil(log2(x_1/theta[5]));
  
  return 1 << i;
}



// Matrix power by squaring: P = A^b --- A is modified
static void matpow_by_squaring(double *A, int n, int b, double *P)
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



// -------------------------------------------------------- 
// Matrix Exponentiation via Pade' Approximations
// -------------------------------------------------------- 



/* r_m(x) = p_m(x) / q_m(x), where
   p_m(x) = sum_{j=0}^m (2m-j)!m!/(2m)!/(m-j)!/j! * x^j
  
   and q_m(x) = p_m(-x)
*/

// (p==q) <= 13
void matexp_pade_fillmats(const int m, const int n, const int i, double *N, double *D, double *B, double *C)
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



static void matexp_pade(const int n, const int p, double *A, double *N)
{
  int i, info = 0;
  int *ipiv;
  double *B, *C, *D;
  
  // Power of A
  B = calloc(n*n, sizeof(double));
  assert(B != NULL);
  
  // Temporary storage for matrix multiplication
  C = malloc(n*n * sizeof(double));
  assert(C != NULL);
  
  D = malloc(n*n * sizeof(double));
  assert(D != NULL);
  
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
  for (i=1; i<=p; i++)
  {
    // C = A*B
    if (i > 1)
      matprod(n, A, B, C);
      
    // Update matrices
    matexp_pade_fillmats(n, n, i, N, D, B, C);
  }
  
  // R <- solve(D) %*% N
  ipiv = calloc(n, sizeof(double));
  assert(ipiv != NULL);
  
  dgesv_(&n, &n, D, &n, ipiv, N, &n, &info);
  
  
  free(B);
  free(C);
  free(D);
  free(ipiv);
}



void matexp(const int n, const int p, const double t, double *x, double *ret)
{
  const int m = matexp_scale_factor(x, n);
  const int one = 1;
  double tmp;
  
  if (m == 0)
  {
    dscal_(&n, &t, x, &one);
    return matexp_pade(n, p, x, ret);
  }
  
  tmp = t / ((double) m);
  dscal_(&n, &tmp, x, &one);
  
  
  matexp_pade(n, p, x, ret);
  
  matcopy(n, ret, x);
  
  matpow_by_squaring(x, n, m, ret);
}

