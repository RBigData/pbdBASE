/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

// Copyright 2013-2015, Schmidt and Chen

#include "pbdBASE.h"


SEXP R_NUMROC(SEXP N, SEXP NB, SEXP IPROC, SEXP NPROCS)
{
  R_INIT;
  SEXP NUM;
  newRvec(NUM, 1, "int");
  
  numrocwrap_(INTP(N), INTP(NB), INTP(IPROC), INTP(NPROCS), INTP(NUM));
  
  R_END;
  return NUM;
}


// -------------------------------------------------------- 
// Linear equations 
// -------------------------------------------------------- 


/* Solving systems of linear equations */
SEXP R_PDGESV(SEXP N, SEXP NRHS, SEXP MXLDIMS, SEXP A, SEXP DESCA, SEXP B, SEXP DESCB)
{
  R_INIT;
  int IJ = 1;
  int * ipiv;
  double *A_cp;
  
  SEXP RET, RET_NAMES, INFO, B_OUT;
  newRvec(INFO, 1, "int");
  newRmat(B_OUT, nrows(B), ncols(B), "dbl");
  
  
  // Copy A and B since pdgesv writes in place
  A_cp = (double *) R_alloc(nrows(A)*ncols(A), sizeof(double));
  //FIXME check returns...
  memcpy(A_cp, DBLP(A), nrows(A)*ncols(A)*sizeof(double));
  memcpy(DBLP(B_OUT), DBLP(B), nrows(B)*ncols(B)*sizeof(double));
  
  
  // Call pdgesv
    ipiv = (int *) R_alloc(INT(MXLDIMS, 0) + INT(DESCA, 5), sizeof(int));
/*  ipiv = (int *) R_alloc(nrows(B) + INT(DESCA, 5), sizeof(int));*/
  
  
  INT(INFO, 0) = 0;
  
  pdgesv_(INTP(N), INTP(NRHS),
    A_cp, &IJ, &IJ, INTP(DESCA), ipiv,
    DBLP(B_OUT), &IJ, &IJ, INTP(DESCB), INTP(INFO));
  
  
  // Manage return
  make_list_names(RET_NAMES, 2, "info", "B");
  make_list(RET, RET_NAMES, 2, INFO, B_OUT);
  
  R_END;
  return RET;
}


/* Matrix inverse */
SEXP R_PDGETRI(SEXP A, SEXP DESCA)
{
  R_INIT;
  int IJ = 1;
  
  SEXP RET, RET_NAMES, INFO, INV;
  
  newRvec(INFO, 1, "int");
  newRmat(INV, nrows(A), ncols(A), "dbl");
  
  
  // Compute inverse
  pdinv_(DBLP(A), &IJ, &IJ, INTP(DESCA), DBLP(INV), INTP(INFO));
  
  
  // Manage return
  make_list_names(RET_NAMES, 2, "info", "A");
  make_list(RET, RET_NAMES, 2, INFO, INV);
  
  R_END;
  return RET;
}



// -------------------------------------------------------- 
// Factorizations 
// -------------------------------------------------------- 

/* Symmetric Eigen */
SEXP R_PDSYEVR(SEXP JOBZ, SEXP UPLO, SEXP N, SEXP A, SEXP DESCA, SEXP DESCZ)
{
  R_INIT;
  SEXP RET, RET_NAMES, INFO, W, Z;
  char range = 'A';
  int IJ = 1;
  int lwork = -1;
  int *iwork;
  int liwork = -1;
  double temp_work = 0;
  double *work;
  double *A_cp;
  double tmp = 0;
  int itmp = 0;
  int m, nz;
  
  newRvec(INFO, 1, "int");
  INT(INFO, 0) = 0;
  newRvec(W, INT(N, 0), "dbl");
  newRmat(Z, nrows(A), ncols(A), "dbl");
  
  /* Query size of workspace */
  // pdsyev_(CHARPT(JOBZ, 0), CHARPT(UPLO, 0), INTP(N),
  //     &tmp, &IJ, &IJ, INTP(DESCA),
  //     &tmp, &tmp, &IJ, &IJ, INTP(DESCZ),
  //     &temp_work, &lwork, INTP(INFO));
  
  pdsyevr_(CHARPT(JOBZ, 0), &range, CHARPT(UPLO, 0), INTP(N),
    &tmp, &IJ, &IJ, INTP(DESCA), &tmp, &tmp, &itmp, &itmp,
    &m, &nz, DBLP(W), DBLP(Z), &IJ, &IJ, INTP(DESCZ),
    &temp_work, &lwork, &liwork, &liwork, INTP(INFO));
  
  /* Allocate workspace and calculate */
  const size_t size = nrows(A)*ncols(A);
  A_cp = (double *) R_alloc(size, sizeof(*A_cp));
  memcpy(A_cp, DBLP(A), size*sizeof(*A_cp));
  
  lwork = (int) temp_work;
  lwork = nonzero(lwork);
  work = (double *) R_alloc(lwork, sizeof(*work));
  
  liwork = nonzero(liwork);
  iwork = (int *) R_alloc(liwork, sizeof(*iwork));
  
  pdsyevr_(CHARPT(JOBZ, 0), &range, CHARPT(UPLO, 0), INTP(N),
    A_cp, &IJ, &IJ, INTP(DESCA), &tmp, &tmp, &itmp, &itmp,
    &m, &nz, DBLP(W), DBLP(Z), &IJ, &IJ, INTP(DESCZ),
    work, &lwork, iwork, &liwork, INTP(INFO));
  
  // pdsyev_(CHARPT(JOBZ, 0), CHARPT(UPLO, 0), INTP(N),
  //     A_cp, &IJ, &IJ, INTP(DESCA),
  //     DBLP(W), DBLP(Z), &IJ, &IJ, INTP(DESCZ),
  //     work, &lwork, INTP(INFO));
  
  
  // Manage return
  make_list_names(RET_NAMES, 3, "values", "vectors", "info");
  make_list(RET, RET_NAMES, 3, W, Z, INFO);
  
  R_END;
  return RET;
} 


/* Non-Symmetric Eigen */
// :[


/* LU factorization */
SEXP R_PDGETRF(SEXP M, SEXP N, SEXP A, SEXP CLDIM, SEXP DESCA, SEXP LIPIV)
{
  R_INIT;
  int *ipiv;
  int IJ = 1;
  SEXP RET, RET_NAMES, INFO, C;
  
  newRvec(INFO, 1, "int");
  newRmat(C, INT(CLDIM, 0), INT(CLDIM, 1), "dbl");
  
  
  // A = LU
  memcpy(DBLP(C), DBLP(A), nrows(A)*ncols(A)*sizeof(double));
  
  INT(INFO, 0) = 0;
  
  INT(LIPIV) = nonzero(INT(LIPIV));
  ipiv = (int*) R_alloc(INT(LIPIV), sizeof(int));
  
  pdgetrf_(INTP(M), INTP(N), DBLP(C), &IJ, &IJ, INTP(DESCA), ipiv, INTP(INFO));
  
  // Manage return
  make_list_names(RET_NAMES, 2, "info", "A");
  make_list(RET, RET_NAMES, 2, INFO, C);
  
  R_END;
  return RET;
}



/* Cholesky */
SEXP R_PDPOTRF(SEXP N, SEXP A, SEXP DESCA, SEXP UPLO)
{
  R_INIT;
  int IJ = 1;
  SEXP RET, RET_NAMES, INFO, C;
  
  newRvec(INFO, 1, "int");
  newRmat(C, nrows(A), ncols(A), "dbl");
  
  // Compute chol
  memcpy(DBLP(C), DBLP(A), nrows(A)*ncols(A)*sizeof(double));
  
  INT(INFO, 0) = 0;
  
  pdpotrf_(STR(UPLO, 0), INTP(N), DBLP(C), &IJ, &IJ, INTP(DESCA), INTP(INFO));
  
  // Manage return
  make_list_names(RET_NAMES, 2, "info", "A");
  make_list(RET, RET_NAMES, 2, INFO, C);
  
  R_END;
  return(RET);
}



// The beast
SEXP R_PDSYEVX(SEXP JOBZ, SEXP RANGE, SEXP N, SEXP A, SEXP DESCA, SEXP VL, SEXP VU, SEXP IL, SEXP IU, SEXP ABSTOL, SEXP ORFAC)
{
  R_INIT;
  char uplo = 'U';
  int IJ = 1;
  int i;
  int m, nz;
  int lwork, liwork, info;
  int descz[9], ldm[2], blacs[5];
  int tmp_liwork;
  int *iwork, *ifail, *iclustr;
  
  double tmp_lwork;
  double *work;
  double *w, *z, *gap;
  double *a;
  
  SEXP RET, RET_NAMES, W, Z, M;
  
  
  // grid and local information
  pdims_(INTEGER(DESCA), ldm, blacs);
  
  ldm[0] = nrows(A);//nonzero(ldm[0]);
  ldm[1] = ncols(A);//nonzero(ldm[1]);
  
  
  // Setup for the setup
  for (i=0; i<9; i++)
    descz[i] = INT(DESCA, i);
  
  w = (double*) R_alloc(INT(N), sizeof(double));
  z = (double*) R_alloc(ldm[0]*ldm[1], sizeof(double));
  gap = (double*) R_alloc(blacs[1]*blacs[2], sizeof(double));
  
  
  a = (double*) R_alloc(ldm[0]*ldm[1], sizeof(double));
  
  memcpy(a, DBLP(A), nrows(A)*ncols(A)*sizeof(double));
  
  ifail = (int*) R_alloc(INT(N, 0), sizeof(int));
  iclustr = (int*) R_alloc(2*blacs[1]*blacs[2], sizeof(int));
  
  
  // Allocate local workspace
  lwork = -1;
  liwork = -1;
  info = 0;
  
  pdsyevx_(CHARPT(JOBZ, 0), CHARPT(RANGE, 0), &uplo, 
    INTP(N), a, &IJ, &IJ, INTP(DESCA), 
    DBLP(VL), DBLP(VU), INTP(IL), INTP(IU), 
    DBLP(ABSTOL), &m, &nz, w, 
    DBLP(ORFAC), z, &IJ, &IJ, descz, 
    &tmp_lwork, &lwork, &tmp_liwork, &liwork, 
    ifail, iclustr, gap, &info);
  
  lwork = nonzero( ((int) tmp_lwork) );
  work = (double*) R_alloc(lwork, sizeof(double));
  
  liwork = nonzero(tmp_liwork);
  iwork = (int*) R_alloc(liwork, sizeof(int));
  
  // Compute eigenvalues
  m = 0;
  info = 0;
  
  pdsyevx_(CHARPT(JOBZ, 0), CHARPT(RANGE, 0), &uplo, 
    INTP(N), a, &IJ, &IJ, INTP(DESCA), 
    DBLP(VL), DBLP(VU), INTP(IL), INTP(IU), 
    DBLP(ABSTOL), &m, &nz, w, 
    DBLP(ORFAC), z, &IJ, &IJ, descz, 
    work, &lwork, iwork, &liwork, 
    ifail, iclustr, gap, &info);
  
  
  newRvec(W, m, "dbl");
  
  for (i=0; i<m; i++)
    DBL(W, i) = w[i];
  
/*  SEXP IFAIL;*/
/*    PROTECT(IFAIL = allocVector(INTSXP, m));*/
/*    for (i=0; i<m; i++)*/
/*        INTEGER(IFAIL)[0] = ifail[i];*/
  
  
  // Manage the return
  if (CHARPT(JOBZ, 0)[0] == 'N') // Only eigenvalues are computed
  {
    make_list_names(RET_NAMES, 1, "values");
    make_list(RET, RET_NAMES, 1, W);
  }
  else // eigenvalues + eigenvectors
  {
    newRmat(Z, ldm[0], ldm[1], "dbl");
    
    for (i=0; i<ldm[0]*ldm[1]; i++)
        DBL(Z, i) = z[i];
    
    newRvec(M, 1, "int");
    INT(M, 0) = m;
    
    make_list_names(RET_NAMES, 3, "values", "vectors", "m");
    make_list(RET, RET_NAMES, 3, W, Z, M);
  }
  
  
  R_END;
  return RET;
}


// -------------------------------------------------------- 
// Auxillary 
// -------------------------------------------------------- 


// Matrix norms
SEXP R_PDLANGE(SEXP TYPE, SEXP M, SEXP N, SEXP A, SEXP DESCA)
{
  R_INIT;
  int IJ = 1;
  
  SEXP VAL;
  newRvec(VAL, 1, "dbl");
  
  matnorm_(DBLP(VAL), STR(TYPE, 0), INTP(M), INTP(N), DBLP(A), &IJ, &IJ, INTP(DESCA));
  
  R_END;
  return VAL;
}



// Condition # estimator for general matrix
SEXP R_PDGECON(SEXP TYPE, SEXP M, SEXP N, SEXP A, SEXP DESCA)
{
  R_INIT;
  int IJ = 1;
  double* cpA;
  int info = 0;
  const int m = nrows(A);
  const int n = ncols(A);
  SEXP RET;
  newRvec(RET, 2, "dbl"); // RET = {cond_num, info}
  
  cpA = malloc(m*n * sizeof(*cpA));
  memcpy(cpA, DBLP(A), m*n*sizeof(*cpA));
  
  // compute inverse of condition number
  condnum_(CHARPT(TYPE, 0), INTP(M), INTP(N), cpA, &IJ, &IJ, INTP(DESCA), DBLP(RET), &info);
  
  DBL(RET, 1) = (double) info;
  
  free(cpA);
  
  R_END;
  return RET;
}



// Condition # estimator for triangular matrix
SEXP R_PDTRCON(SEXP TYPE, SEXP UPLO, SEXP DIAG, SEXP N, SEXP A, SEXP DESCA)
{
  R_INIT;
  double* work;
  double tmp;
  int* iwork;
  int lwork, liwork, info = 0;
  int IJ = 1, in1 = -1;
  
  SEXP RET;
  newRvec(RET, 2, "dbl");
  
  // workspace query and allocate work vectors
  pdtrcon_(CHARPT(TYPE, 0), CHARPT(UPLO, 0), CHARPT(DIAG, 0),
    INTP(N), DBLP(A), &IJ, &IJ, INTP(DESCA), DBLP(RET), 
    &tmp, &in1, &liwork, &in1, &info);
  
  lwork = (int) tmp;
  work = malloc(lwork * sizeof(*work));
  iwork = malloc(liwork * sizeof(*iwork));
  
  // compute inverse of condition number
  info = 0;
  pdtrcon_(CHARPT(TYPE, 0), CHARPT(UPLO, 0), CHARPT(DIAG, 0),
    INTP(N), DBLP(A), &IJ, &IJ, INTP(DESCA), DBLP(RET), 
    work, &lwork, iwork, &liwork, &info);
  
  DBL(RET, 1) = (double) info;
  
  free(work);
  free(iwork);
  
  R_END;
  return RET;
}
