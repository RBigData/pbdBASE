/* Copyright (C) 2013-2015 Drew Schmidt. 
 *               2015 Wei-Chen Chen
 * All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

  * Redistributions of source code must retain the above copyright notice, this
    list of conditions and the following disclaimer.
  * Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions and the following disclaimer in the documentation
    and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. */


/* Matrix exponentiation algorithm from:
   "New Scaling and Squaring Algorithm for the Matrix Exponential", by
   Awad H. Al-Mohy and Nicholas J. Higham, August 2009
*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

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
  unsigned int i;
  
  for (i=0; i<n*n; i++)
    a[i] = 0.0;
  
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
      tmp += fabs(x[i + j*n]);
    
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
  
  i = (int) ceil(log2(x_1/theta[4]));
  
  return 1 << i;
}



// Matrix power by squaring: P = A^b (A is garbage on exit)
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
    
    b >>= 1;
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

// Workhorse for matexp_pade
void matexp_pade_fillmats(const int m, const int n, const int i, double *restrict N, double *restrict D, double *restrict B, const double *restrict C)
{
  int j;
  const double tmp = matexp_pade_coefs[i];
  double tmpj;
  const int sgn = SGNEXP(-1, i);
  
    /* Performs the following actions:
        B = C
        N = pade_coef[i] * C
        D = (-1)^j * pade_coef[i] * C
    */
    for (j=0; j<m*n; j++)
    {
      tmpj = C[j];
      B[j] = tmpj;
      
      tmpj *= tmp;
      
      N[j] += tmpj;
      D[j] += sgn*tmpj;
    }
}



// Exponentiation via Pade' expansion
static void matexp_pade(int n, const int p, double *A, double *N)
{
  int i, info = 0;
  int *ipiv;
  double *B, *C, *D;
  
  // Power of A
  B = calloc(n*n, sizeof(*B));
  assert(B != NULL);
  
  // Temporary storage for matrix multiplication
  C = calloc(n*n, sizeof(*C));
  assert(C != NULL);
  
  D = calloc(n*n, sizeof(*D));
  assert(D != NULL);
  
  matcopy(n, A, C);
  
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
  
  // R <- inverse(D) %*% N
  ipiv = calloc(n, sizeof(*ipiv));
  assert(ipiv != NULL);
  
  dgesv_(&n, &n, D, &n, ipiv, N, &n, &info);
  
  
  free(B);
  free(C);
  free(D);
  free(ipiv);
}



/*
  n       Number of rows/cols of (square) matrix x.
  
  p       Order of the Pade' approximation. 0 < p <= 13.
  
  t       Scaling factor for x (t=1 canonical).
  
  x       Input (square) matrix.  On function exit, the values
          in x are garbage.
  
  ret     On exit, ret = expm(x).
*/

void matexp(int n, const int p, double *x, double *ret)
{
  int m;
  int nn = n*n;
  int one = 1;
  double tmp;
  
  m = matexp_scale_factor(x, n);
  
  if (m == 0)
    matexp_pade(n, p, x, ret);
  else
  {
    tmp = 1. / ((double) m);
    dscal_(&nn, &tmp, x, &one);
  
    matexp_pade(n, p, x, ret);
  
    matcopy(n, ret, x);
  
    matpow_by_squaring(x, n, m, ret);
  }
}
