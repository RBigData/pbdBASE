/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

// Copyright 2013, Schmidt

#include <ctype.h>

#include "scalapach.h"
#include "errors.h"


void pdgemm_(char *transa, char *transb, int *m, int *n, int *k, double *alpha, double *a, int *ia, int *ja, int *desca, double *b, int *ib, int *jb, int *descb, double *beta, double *c, int *ic, int *jc, int *descc);
void pdlacpy_(char *uplo, int *m, int *n, double *a, int *ia, int *ja, int *desca, double *b, int *ib, int *jb, int *descb);


ddmatrix *p_newmat(double *Data, int *desc, int *ldim)
{
  ddmatrix *a = malloc(sizeof(struct ddmatrix));
  
  assert(a != NULL);
  
  P_DATA(a) = Data;
  P_DESC(a) = desc;
  P_LDIM(a) = ldim;
  
  return a;
}



// Copy a ONTO b, i.e. b = a
int p_matcopy(ddmatrix *a, ddmatrix *b)
{
  int i;
  const char uplo = 'A';
  const int ij = 1;
  
  // Data
  pdlacpy_(&uplo, P_DESC(a)[2], P_DESC(a)[3], 
           P_DATA(a), &ij, &ij, P_DESC(a), 
           P_DATA(b), &ij, &ij, P_DESC(b));
  
  // Descriptor array
  for (i=0; i<9; i++)
    P_DESC(b)[i] = P_DESC(a)[i];
  
  // ldim
  P_LDIM(b)[0] = P_LDIM(a)[0];
  P_LDIM(b)[1] = P_LDIM(a)[1];
  
  
  return NOPROBLEMO;
}



int p_delmat(ddmatrix *a)
{
  assert(a != NULL);
  
  free(P_DATA(a));
  free(a);
  
  return NOPROBLEMO;
}



void p_printmat(ddmatrix *x)
{
  // If comm.rank() == 0
  PRINT("    Dense Distributed Matrix\n");
  PRINT("--------------------------------\n");
/*  PRINT("Process grid:          %dx%d\n", blacs[2], blacs[3]);*/
  PRINT("Global dimension:       %dx%d\n", P_DIM(x)[0], P_DIM(x)[1]);
  PRINT("(max) Local dimension:  %dx%d\n", P_LDIM(x)[0], P_LDIM(x)[1]);
  PRINT("Blocking:               %dx%d\n", P_DESC(x)[4], P_DESC(x)[5]);
  PRINT("BLACS ICTXT:            %d\n", P_DESC(x)[1]);
}



// c = op(a)*op(b)
int p_matprod(char opa, ddmatrix *a, char opb, ddmatrix *b, ddmatrix *c)
{
  const int ij = 1;
  const double one = 1.0, zero = 0.0;
  int m, n, k, k_check;
  
  
  // Quick return if possible
  if (P_DIM(a)[0] < 1 || P_DIM(a)[1] < 1 || P_DIM(b)[0] < 1 || P_DIM(b)[1] < 1)
    return DEGENERATE;
  
  
  opa = toupper(opa);
  opb = toupper(opb);
  
  if (opa == 'N')
  {
    m = P_DIM(a)[0];
    k = P_DIM(a)[1];
  }
  else
  {
    m = P_DIM(a)[1];
    k = P_DIM(a)[0];
  }
  
  if (opb == 'N')
  {
    n = P_DIM(b)[1];
    k_check = P_DIM(b)[0];
  }
  else
  {
    n = P_DIM(b)[0];
    k_check = P_DIM(b)[1];
  }
  
  
  // Check conformality
  if (k != k_check)
    return NONCONFORMABLE;
  
  
  // Compute product
  pdgemm_(&opa, &opb, &m, &n, &k, &one, 
          P_DATA(a), &ij, &ij, P_DESC(a), 
          P_DATA(b), &ij, &ij, P_DESC(descb), &zero, 
          P_DATA(c), &ij, &ij, P_DESC(descc));
  
  
  return NOPROBLEMO;
}


