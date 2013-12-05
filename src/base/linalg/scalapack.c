/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

// Copyright 2013, Schmidt



// -------------------------------------------------------- 
// Linear equations 
// -------------------------------------------------------- 

pdgesv_(int *n, int *nrhs, double *a, int *ia, int *ja, int *desca, int *ipiv, double *b, int *ib, int *jb, int *descb, int *info);
void pdinv_(double *x, int *descx, double *inv, int *info);


// Solve A*x = b
int p_solve(ddmatrix *A, ddmatrix *b, ddmatrix *x)
{
  int info;
  const int ij = 1;
  ddmatrix *A_;
  int * ipiv;
  
  ipiv = (int *) malloc( + P_DESC(A)[5] * sizeof(int)); // FIXME
  
  p_matcopy(A, A_);
  p_matcopy(b, x);
  
  pdgesv(P_NROWS(A), P_COLS(x), P_DATA(A), &ij, &ij, P_DESC(a), 
         &ipiv, P_DATA(x), &ij, &ij, P_DESC(x), &info);
  
  
  p_delmat(A_);
  free(ipiv);
  
  return info;
}



// Matrix inverse 
int p_matinv(ddmatrix *x, ddmatrix *inv)
{
  int info;
  const int ij = 1;
  
  pdinv_(P_DATA(x), &ij, &ij, P_DESC(x), P_DATA(inv), &info);
  
  return info;
}

#if 0

// -------------------------------------------------------- 
// Factorizations 
// -------------------------------------------------------- 


/* SVD */
PDGESVD


/* Symmetric Eigen */
PDSYEV

/* Non-Symmetric Eigen */
// :[


/* LU factorization */
PDGETRF


/* Cholesky */
PDPOTRF



// The beast
PDSYEVX


// -------------------------------------------------------- 
// Auxillary 
// -------------------------------------------------------- 


// Matrix norms
PDLANGE(SEXP TYPE, SEXP M, SEXP N, SEXP A, SEXP DESCA)
{
    const int IJ = 1;
    double *work;
    
    SEXP VAL;
    PROTECT(VAL = allocVector(REALSXP, 1));
    
    F77_CALL(matnorm)(REAL(VAL), CHARPT(TYPE, 0), INTEGER(M),
        INTEGER(N), REAL(A), &IJ, &IJ, INTEGER(DESCA));
    
    UNPROTECT(1);
    return(VAL);
}


// Condition # estimator for general matrix
SEXP R_PDGECON(SEXP TYPE, SEXP M, SEXP N, SEXP A, SEXP DESCA, SEXP ALDIM)
{
    const int IJ = 1;
    double* cpA;
    int i, info = 0;
    int* pt_ALDIM = INTEGER(ALDIM);
    
    SEXP RET;
    PROTECT(RET = allocVector(REALSXP, 2));
    
    // Copy A
    i = pt_ALDIM[0] * pt_ALDIM[1];
    cpA = R_alloc(i, sizeof(double));
    memcpy(cpA, REAL(A), i * sizeof(double));
    
    // compute inverse of condition number
    F77_CALL(condnum)(CHARPT(TYPE, 0), INTEGER(M), INTEGER(N), cpA, 
        &IJ, &IJ, INTEGER(DESCA), REAL(RET), &info);
    
    REAL(RET)[1] = (double) info;
    
    UNPROTECT(1);
    return(RET);
}



// Condition # estimator for triangular matrix
SEXP R_PDTRCON(SEXP TYPE, SEXP UPLO, SEXP DIAG, SEXP N, SEXP A, SEXP DESCA)
{
    double* work;
    double tmp;
    int* iwork;
    int i, lwork, liwork, info = 0;
    // consts 
    const int IJ = 1, in1 = -1;
    // R objects
    SEXP RET;
    PROTECT(RET = allocVector(REALSXP, 2));
    
    // workspace query and allocate work vectors
    F77_CALL(pdtrcon)(CHARPT(TYPE, 0), CHARPT(UPLO, 0), CHARPT(DIAG, 0),
        INTEGER(N), REAL(A), &IJ, &IJ, INTEGER(DESCA), REAL(RET), 
        &tmp, &in1, &liwork, &in1, &info);
    
    lwork = (int) tmp;
    work = (double *) R_alloc(lwork, sizeof(double));
    iwork = (int *) R_alloc(liwork, sizeof(int));
    
    // compute inverse of condition number
    info = 0;
    F77_CALL(pdtrcon)(CHARPT(TYPE, 0), CHARPT(UPLO, 0), CHARPT(DIAG, 0),
        INTEGER(N), REAL(A), &IJ, &IJ, INTEGER(DESCA), REAL(RET), 
        work, &lwork, iwork, &liwork, &info);
    
    REAL(RET)[1] = (double) info;
    
    UNPROTECT(1);
    return(RET);
}

#endif

