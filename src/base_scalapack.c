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
  RET_NAMES = make_list_names(2, "info", "B");
  RET = make_list(RET_NAMES, 2, INFO, B_OUT);
  
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
  RET_NAMES = make_list_names(2, "info", "A");
  RET = make_list(RET_NAMES, 2, INFO, INV);
  
  R_END;
  return RET;
}



// -------------------------------------------------------- 
// Factorizations 
// -------------------------------------------------------- 


/* SVD */
SEXP R_PDGESVD(SEXP M, SEXP N, SEXP ASIZE, SEXP A, SEXP DESCA, 
    SEXP ULDIM, SEXP DESCU, SEXP VTLDIM, SEXP DESCVT, SEXP JOBU, SEXP JOBVT, 
    SEXP INPLACE)
{
  R_INIT;
  double *A_OUT;
  int temp_IJ = 1, temp_lwork = -1;
  double temp_A = 0, temp_work = 0, *WORK;
  SEXP RET, RET_NAMES, INFO, D, U, VT;
  
  newRvec(INFO, 1, "int");
  newRvec(D, INT(ASIZE, 0), "dbl");
  newRmat(U, INT(ULDIM, 0), INT(ULDIM, 1), "dbl");
  newRmat(VT, INT(VTLDIM, 0), INT(VTLDIM, 1), "dbl");
  
  
  // Query size of workspace
  INT(INFO, 0) = 0;
  
  pdgesvd_(STR(JOBU, 0), STR(JOBVT, 0),
    INTP(M), INTP(N),
    &temp_A, &temp_IJ, &temp_IJ, INTP(DESCA),
    &temp_A, &temp_A, &temp_IJ, &temp_IJ, INTP(DESCU),
    &temp_A, &temp_IJ, &temp_IJ, INTP(DESCVT),
    &temp_work, &temp_lwork, INTP(INFO));
      
  // Allocate work vector and calculate svd
  temp_lwork = (int) temp_work;
  temp_lwork = nonzero(temp_lwork);
  
  WORK = (double *) R_alloc(temp_lwork, sizeof(double));
  
  A_OUT = (double *) R_alloc(nrows(A)*ncols(A), sizeof(double));
  memcpy(A_OUT, REAL(A), nrows(A)*ncols(A)*sizeof(double));
  
  pdgesvd_(STR(JOBU, 0), STR(JOBVT, 0),
    INTP(M), INTP(N),
    A_OUT, &temp_IJ, &temp_IJ, INTP(DESCA),
    REAL(D), REAL(U), &temp_IJ, &temp_IJ, INTP(DESCU),
    REAL(VT), &temp_IJ, &temp_IJ, INTP(DESCVT),
    WORK, &temp_lwork, INTP(INFO));
  
  // Manage return
  RET_NAMES = make_list_names(4, "info", "d", "u", "vt");
  RET = make_list(RET_NAMES, 4, INFO, D, U, VT);
  
  R_END;
  return RET;
} 


/* Symmetric Eigen */
SEXP R_PDSYEV(SEXP JOBZ, SEXP UPLO, SEXP N, SEXP A, SEXP DESCA, SEXP ZLDIM, SEXP DESCZ)
{
  R_INIT;
  SEXP RET, RET_NAMES, INFO, W, Z;
  int temp_IJ = 1, temp_lwork = -1;
  double temp_A = 0, temp_work = 0, *WORK;
  
  newRvec(INFO, 1, "int");
  newRvec(W, INT(N, 0), "dbl");
  newRmat(Z, INT(ZLDIM, 0), INT(ZLDIM, 1), "dbl");
  
  
  /* Query size of workspace */
  INT(INFO, 0) = 0;
  
  pdsyev_(CHARPT(JOBZ, 0), CHARPT(UPLO, 0), INTP(N),
      &temp_A, &temp_IJ, &temp_IJ, INTP(DESCA),
      &temp_A, &temp_A, &temp_IJ, &temp_IJ, INTP(DESCZ),
      &temp_work, &temp_lwork, INTP(INFO));
      
  /* Allocate work vector and calculate svd */
  temp_lwork = (int) temp_work;
  temp_lwork = nonzero(temp_lwork);
  
  WORK = (double *) R_alloc(temp_lwork, sizeof(double));
  
  pdsyev_(CHARPT(JOBZ, 0), CHARPT(UPLO, 0), INTP(N),
      DBLP(A), &temp_IJ, &temp_IJ, INTP(DESCA),
      DBLP(W), DBLP(Z), &temp_IJ, &temp_IJ, INTP(DESCZ),
      WORK, &temp_lwork, INTP(INFO));
  
  
  // Manage return
  RET_NAMES = make_list_names(3, "values", "vectors", "info");
  RET = make_list(RET_NAMES, 3, W, Z, INFO);
  
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
  ipiv = R_alloc(INT(LIPIV), sizeof(int));
  
  pdgetrf_(INTP(M), INTP(N), DBLP(C), &IJ, &IJ, INTP(DESCA), ipiv, INTP(INFO));
  
  // Manage return
  RET_NAMES = make_list_names(2, "info", "A");
  RET = make_list(RET_NAMES, 2, INFO, C);
  
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
  RET_NAMES = make_list_names(2, "info", "A");
  RET = make_list(RET_NAMES, 2, INFO, C);
  
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
  
  w = R_alloc(INT(N), sizeof(double));
  z = R_alloc(ldm[0]*ldm[1], sizeof(double));
  gap = R_alloc(blacs[1]*blacs[2], sizeof(double));
  
  
  a = R_alloc(ldm[0]*ldm[1], sizeof(double));
  
  memcpy(a, DBLP(A), nrows(A)*ncols(A)*sizeof(double));
  
  ifail = R_alloc(INT(N, 0), sizeof(int));
  iclustr = R_alloc(2*blacs[1]*blacs[2], sizeof(int));
  
  
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
  work = R_alloc(lwork, sizeof(double));
  
  liwork = nonzero(tmp_liwork);
  iwork = R_alloc(liwork, sizeof(int));
  
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
    RET_NAMES = make_list_names(1, "values");
    RET = make_list(RET_NAMES, 1, W);
  }
  else // eigenvalues + eigenvectors
  {
    newRmat(Z, ldm[0], ldm[1], "dbl");
    
    for (i=0; i<ldm[0]*ldm[1]; i++)
        DBL(Z, i) = z[i];
    
    newRvec(M, 1, "int");
    INT(M, 0) = m;
    
    RET_NAMES = make_list_names(3, "values", "vectors", "m");
    RET = make_list(RET_NAMES, 3, W, Z, M);
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
  
  SEXP RET;
  newRvec(RET, 2, "dbl"); // RET = {cond_num, info}
  
  
  // Copy A
  cpA = R_alloc(nrows(A)*ncols(A), sizeof(double));
  memcpy(cpA, REAL(A), nrows(A)*ncols(A)*sizeof(double));
  
  
  // compute inverse of condition number
  condnum_(CHARPT(TYPE, 0), INTP(M), INTP(N), cpA, 
    &IJ, &IJ, INTP(DESCA), DBLP(RET), &info);
  
  
  DBL(RET, 1) = (double) info;
  
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
  work = (double *) R_alloc(lwork, sizeof(double));
  iwork = (int *) R_alloc(liwork, sizeof(int));
  
  // compute inverse of condition number
  info = 0;
  pdtrcon_(CHARPT(TYPE, 0), CHARPT(UPLO, 0), CHARPT(DIAG, 0),
    INTP(N), DBLP(A), &IJ, &IJ, INTP(DESCA), DBLP(RET), 
    work, &lwork, iwork, &liwork, &info);
  
  DBL(RET, 1) = (double) info;
  
  R_END;
  return RET;
}


