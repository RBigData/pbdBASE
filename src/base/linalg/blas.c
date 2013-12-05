/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

// Copyright 2013, Schmidt

#include <ctype.h>

#include "scalapack.h"
#include "errors.h"


void dgemm_(char *transa, char *transb, int *m, int *n, int *k, double *alpha, double *a, int *lda, double *b, int *ldb, double *beta, double *c, int *ldc);
void dlacpy_(char *uplo, int *m, int *n, double *a, int *ia, int *ja, int *desca, double *b, int *ib, int *jb, int *descb);


// Copy a ONTO b, i.e. b = a
int p_matcopy(ddmatrix *a, ddmatrix *b)
{
  const char uplo = 'A';
  const int ij = 1;
  
  //FIXME change matrix attributes as well
  
  dlacpy_(&uplo, DIM(a), &DIM(a)[1],
           DATA(a), &LDA(a), 
           DATA(b), &LDA(b));
  
  return NOPROBLEMO;
}



// c = op(a)*op(b)
int p_matprod(char opa, matrix *a, char opb, matrix *b, matrix *c)
{
  const int ij = 1;
  const double one = 1.0, zero = 0.0;
  int m, n, k, k_check;
  
  
  // Quick return if possible
  if (DIM(a)[0] < 1 || DIM(a)[1] < 1 || DIM(b)[0] < 1 || DIM(b)[1] < 1)
    return DEGENERATE;
  
  
  opa = toupper(opa);
  opb = toupper(opb);
  
  if (opa == 'N')
  {
    m = DIM(a)[0];
    k = DIM(a)[1];
  }
  else
  {
    m = DIM(a)[1];
    k = DIM(a)[0];
  }
  
  if (opb == 'N')
  {
    n = DIM(b)[1];
    k_check = DIM(b)[0];
  }
  else
  {
    n = DIM(b)[0];
    k_check = DIM(b)[1];
  }
  
  
  // Check conformality
  if (k != k_check)
    return NONCONFORMABLE;
  
  
  // Compute product
  dgemm_(&opa, &opb, &m, &n, &k, &one, 
          DATA(a), LDA(a), 
          DATA(b), LDA(b), &zero, 
          DATA(c), LDA(c));
  
  
  return NOPROBLEMO;
}


