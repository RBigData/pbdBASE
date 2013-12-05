#include <R.h>
#include <Rinternals.h>
#include "base_global.h"

#define TRUE 1
#define FALSE 0


SEXP R_NUMROC(SEXP N, SEXP NB, SEXP IPROC, SEXP NPROCS)
{
    SEXP NUM;
    PROTECT(NUM = allocVector(INTSXP, 1));
    
    numrocwrap_(INTEGER(N), INTEGER(NB), INTEGER(IPROC), INTEGER(NPROCS), INTEGER(NUM));
    
    UNPROTECT(1);
    return NUM;
}


// -------------------------------------------------------- 
// Linear equations 
// -------------------------------------------------------- 


/* Solving systems of linear equations */
SEXP R_PDGESV(SEXP N, SEXP NRHS, SEXP MXLDIMS, SEXP A, SEXP ALDIM, SEXP DESCA,
    SEXP B, SEXP BLDIM, SEXP DESCB)
{
    int i, *pt_ALDIM = INTEGER(ALDIM), *pt_BLDIM = INTEGER(BLDIM);
    double *pt_ORG, *pt_COPY, *A_OUT;
    SEXP RET, RET_NAMES, INFO, B_OUT;
    
    /* Protect R objects. */
    PROTECT(RET = allocVector(VECSXP, 2));
    PROTECT(RET_NAMES = allocVector(STRSXP, 2));
    PROTECT(INFO = allocVector(INTSXP, 1));
    PROTECT(B_OUT = allocMatrix(REALSXP, pt_BLDIM[0], pt_BLDIM[1]));
    
    SET_VECTOR_ELT(RET, 0, INFO);
    SET_VECTOR_ELT(RET, 1, B_OUT);
    SET_STRING_ELT(RET_NAMES, 0, mkChar("info")); 
    SET_STRING_ELT(RET_NAMES, 1, mkChar("B")); 
    setAttrib(RET, R_NamesSymbol, RET_NAMES);
    
    /* Copy A and B since pdgesv writes in place */
    A_OUT = (double *) R_alloc(pt_ALDIM[0] * pt_ALDIM[1], sizeof(double));
    pt_ORG = REAL(A);
    pt_COPY = A_OUT;
    for(i = 0; i < pt_ALDIM[0] * pt_ALDIM[1]; i++){
        *pt_COPY = *pt_ORG;
        pt_ORG++;
        pt_COPY++;
    }
    
    pt_ORG = REAL(B);
    pt_COPY = REAL(B_OUT);
    for(i = 0; i < pt_BLDIM[0] * pt_BLDIM[1]; i++){
        *pt_COPY = *pt_ORG;
        pt_ORG++;
        pt_COPY++;
    }
    
    const int IJ = 1;
    int * ipiv;
    ipiv = (int *) R_alloc(INTEGER(MXLDIMS)[0] + INTEGER(DESCA)[5], sizeof(int));
    
    /* Set info */
    INTEGER(INFO)[0] = 0;
    
    F77_CALL(pdgesv)(INTEGER(N), INTEGER(NRHS),
        A_OUT, &IJ, &IJ, INTEGER(DESCA), ipiv,
        REAL(B_OUT), &IJ, &IJ, INTEGER(DESCB), INTEGER(INFO));
    
    /* Return. */
    UNPROTECT(4);
    return(RET);
} /* End of R_PDGESV(). */


/* Matrix inverse */
SEXP R_PDGETRI(SEXP A, SEXP DESCA)
{
    const int ij = 1;
    
    /* R objects. */
    SEXP RET, RET_NAMES, INFO, INV;
    
    PROTECT(RET = allocVector(VECSXP, 2));
    PROTECT(RET_NAMES = allocVector(STRSXP, 2));
    PROTECT(INFO = allocVector(INTSXP, 1));
    PROTECT(INV = allocMatrix(REALSXP, nrows(A), ncols(A)));
    
    setAttrib(RET, R_NamesSymbol, RET_NAMES);
    
    
    // Compute inverse
    pdinv_(REAL(A), &ij, &ij, INTEGER(DESCA), REAL(INV), INTEGER(INFO));
    
    
    // Return
    SET_VECTOR_ELT(RET, 0, INFO);
    SET_VECTOR_ELT(RET, 1, INV);
    SET_STRING_ELT(RET_NAMES, 0, mkChar("info")); 
    SET_STRING_ELT(RET_NAMES, 1, mkChar("A")); 
    
    UNPROTECT(4);
    return(RET);
} /* End of R_PDGETRI(). */



// -------------------------------------------------------- 
// Factorizations 
// -------------------------------------------------------- 


/* SVD */
SEXP R_PDGESVD(SEXP M, SEXP N, SEXP ASIZE, SEXP A, SEXP DESCA, SEXP ALDIM, 
    SEXP ULDIM, SEXP DESCU, SEXP VTLDIM, SEXP DESCVT, SEXP JOBU, SEXP JOBVT, 
    SEXP INPLACE)
{
    int i, *pt_ALDIM = INTEGER(ALDIM);
    double *A_OUT;
    SEXP RET, RET_NAMES, INFO, D, U, VT;

    /* Extra needed. */
    int temp_IJ = 1, temp_lwork = -1;
    double temp_A = 0, temp_work = 0, *WORK;

    /* Protect R objects. */
    PROTECT(A);
    
    PROTECT(RET = allocVector(VECSXP, 4));
    PROTECT(RET_NAMES = allocVector(STRSXP, 4));
    
    PROTECT(INFO = allocVector(INTSXP, 1));
    PROTECT(D = allocVector(REALSXP, INTEGER(ASIZE)[0]));
    PROTECT(U = allocMatrix(REALSXP, INTEGER(ULDIM)[0], INTEGER(ULDIM)[1]));
    PROTECT(VT = allocMatrix(REALSXP,
            INTEGER(VTLDIM)[0], INTEGER(VTLDIM)[1]));
    
    SET_VECTOR_ELT(RET, 0, INFO);
    SET_VECTOR_ELT(RET, 1, D);
    SET_VECTOR_ELT(RET, 2, U);
    SET_VECTOR_ELT(RET, 3, VT);
    
    SET_STRING_ELT(RET_NAMES, 0, mkChar("info")); 
    SET_STRING_ELT(RET_NAMES, 1, mkChar("d")); 
    SET_STRING_ELT(RET_NAMES, 2, mkChar("u")); 
    SET_STRING_ELT(RET_NAMES, 3, mkChar("vt")); 
    setAttrib(RET, R_NamesSymbol, RET_NAMES);
    
    
    /* Query size of workspace */
    INTEGER(INFO)[0] = 0;
    F77_CALL(pdgesvd)(CHARPT(JOBU, 0), CHARPT(JOBVT, 0),
        INTEGER(M), INTEGER(N),
        &temp_A, &temp_IJ, &temp_IJ, INTEGER(DESCA),
        &temp_A, &temp_A, &temp_IJ, &temp_IJ, INTEGER(DESCU),
        &temp_A, &temp_IJ, &temp_IJ, INTEGER(DESCVT),
        &temp_work, &temp_lwork, INTEGER(INFO));
        
    /* Allocate work vector and calculate svd */
    temp_lwork = (int) temp_work;
    temp_lwork = nonzero(temp_lwork);
    
    WORK = (double *) R_alloc(temp_lwork, sizeof(double));
    
    INTEGER(INFO)[0] = 0;
    
/*    if (CHARPT(INPLACE, 0)[0] == 'N'){*/
        /* Make copy of original data, since pdgesvd destroys it */
        i = pt_ALDIM[0] * pt_ALDIM[1];
        A_OUT = (double *) R_alloc(i, sizeof(double));
        memcpy(A_OUT, REAL(A), i * sizeof(double));
        
        pdgesvd_(CHARPT(JOBU, 0), CHARPT(JOBVT, 0),
            INTEGER(M), INTEGER(N),
            A_OUT, &temp_IJ, &temp_IJ, INTEGER(DESCA),
            REAL(D), REAL(U), &temp_IJ, &temp_IJ, INTEGER(DESCU),
            REAL(VT), &temp_IJ, &temp_IJ, INTEGER(DESCVT),
            WORK, &temp_lwork, INTEGER(INFO));
/*    }*/
/*    else {*/
/*        pdgesvd_(CHARPT(JOBU, 0), CHARPT(JOBVT, 0),*/
/*            INTEGER(M), INTEGER(N),*/
/*            REAL(A), &temp_IJ, &temp_IJ, INTEGER(DESCA),*/
/*            REAL(D), REAL(U), &temp_IJ, &temp_IJ, INTEGER(DESCU),*/
/*            REAL(VT), &temp_IJ, &temp_IJ, INTEGER(DESCVT),*/
/*            WORK, &temp_lwork, INTEGER(INFO));*/
/*    }*/
    
    UNPROTECT(7);
    
    return(RET);
} 


/* Symmetric Eigen */
SEXP R_PDSYEV(SEXP JOBZ, SEXP UPLO, SEXP N, SEXP A, SEXP DESCA, SEXP ALDIM, SEXP ZLDIM, SEXP DESCZ)
{
    int i, *pt_ALDIM = INTEGER(ALDIM);
    double *A_OUT;
    SEXP RET, RET_NAMES, INFO, W, Z;
    
    /* Extra needed. */
    int temp_IJ = 1, temp_lwork = -1;
    double temp_A = 0, temp_work = 0, *WORK;
    
    /* Protect R objects. */
    PROTECT(A);
    
    PROTECT(RET = allocVector(VECSXP, 3));
    PROTECT(RET_NAMES = allocVector(STRSXP, 3));
    
    PROTECT(W = allocVector(REALSXP, INTEGER(N)[0]));
    PROTECT(Z = allocMatrix(REALSXP, INTEGER(ZLDIM)[0], INTEGER(ZLDIM)[1]));
    PROTECT(INFO = allocVector(INTSXP, 1));
    
    SET_VECTOR_ELT(RET, 0, W);
    SET_VECTOR_ELT(RET, 1, Z);
    SET_VECTOR_ELT(RET, 2, INFO);
    
    SET_STRING_ELT(RET_NAMES, 0, mkChar("values")); 
    SET_STRING_ELT(RET_NAMES, 1, mkChar("vectors")); 
    SET_STRING_ELT(RET_NAMES, 2, mkChar("info")); 
    setAttrib(RET, R_NamesSymbol, RET_NAMES);
    
    /* Query size of workspace */
    INTEGER(INFO)[0] = 0;
    pdsyev_(CHARPT(JOBZ, 0), CHARPT(UPLO, 0), INTEGER(N),
        &temp_A, &temp_IJ, &temp_IJ, INTEGER(DESCA),
        &temp_A, &temp_A, &temp_IJ, &temp_IJ, INTEGER(DESCZ),
        &temp_work, &temp_lwork, INTEGER(INFO));
        
    /* Allocate work vector and calculate svd */
    temp_lwork = (int) temp_work;
    temp_lwork = nonzero(temp_lwork);
    
    WORK = (double *) R_alloc(temp_lwork, sizeof(double));
    
    INTEGER(INFO)[0] = 0;
    
    pdsyev_(CHARPT(JOBZ, 0), CHARPT(UPLO, 0), INTEGER(N),
        REAL(A), &temp_IJ, &temp_IJ, INTEGER(DESCA),
        REAL(W), REAL(Z), &temp_IJ, &temp_IJ, INTEGER(DESCZ),
        WORK, &temp_lwork, INTEGER(INFO));
    
    UNPROTECT(6);
    
    return(RET);
} 


/* Non-Symmetric Eigen */
// :[


/* LU factorization */
SEXP R_PDGETRF(SEXP M, SEXP N, SEXP A, SEXP CLDIM, SEXP DESCA, SEXP LIPIV)
{
    int i, *pt_CLDIM = INTEGER(CLDIM), *ipiv;
    double *pt_A, *pt_C;
    const int IJ = 1;
    SEXP RET, RET_NAMES, INFO, C;
    
    /* Protect R objects. */
    PROTECT(RET = allocVector(VECSXP, 2));
    PROTECT(RET_NAMES = allocVector(STRSXP, 2));
    PROTECT(INFO = allocVector(INTSXP, 1));
    PROTECT(C = allocMatrix(REALSXP, pt_CLDIM[0], pt_CLDIM[1]));
    
    SET_VECTOR_ELT(RET, 0, INFO);
    SET_VECTOR_ELT(RET, 1, C);
    SET_STRING_ELT(RET_NAMES, 0, mkChar("info")); 
    SET_STRING_ELT(RET_NAMES, 1, mkChar("A")); 
    setAttrib(RET, R_NamesSymbol, RET_NAMES);
    
    /* Set INFO and Copy A -> C. */
    INTEGER(INFO)[0] = 0;
    pt_A = REAL(A);
    pt_C = REAL(C);
    for(i = 0; i < pt_CLDIM[0] * pt_CLDIM[1]; i++){
        *pt_C = *pt_A;
        pt_A++;
        pt_C++;
    }
    
    LIPIV = nonzero(LIPIV);
    ipiv = (int *) R_alloc(INTEGER(LIPIV), sizeof(int));
    
    INTEGER(INFO)[0] = 0;
    F77_CALL(pdgetrf)(INTEGER(M), INTEGER(N), REAL(C), 
        &IJ, &IJ, INTEGER(DESCA), ipiv, INTEGER(INFO));
    
    /* Return. */
    UNPROTECT(4);
    return(RET);
} /* End of R_PDGETRF(). */



/* Cholesky */
SEXP R_PDPOTRF(SEXP N, SEXP A, SEXP CLDIM, SEXP DESCA, SEXP UPLO)
{
    int i, *pt_CLDIM = INTEGER(CLDIM);
    const int IJ = 1;
    double *pt_A, *pt_C;
    SEXP RET, RET_NAMES, INFO, C;
    
    /* Protect R objects. */
    PROTECT(RET = allocVector(VECSXP, 2));
    PROTECT(RET_NAMES = allocVector(STRSXP, 2));
    PROTECT(INFO = allocVector(INTSXP, 1));
    PROTECT(C = allocMatrix(REALSXP, pt_CLDIM[0], pt_CLDIM[1]));
    
    SET_VECTOR_ELT(RET, 0, INFO);
    SET_VECTOR_ELT(RET, 1, C);
    SET_STRING_ELT(RET_NAMES, 0, mkChar("info")); 
    SET_STRING_ELT(RET_NAMES, 1, mkChar("A")); 
    setAttrib(RET, R_NamesSymbol, RET_NAMES);
    
    /* Copy A -> C and set INFO and return R objects. */
    INTEGER(INFO)[0] = 0;
    pt_A = REAL(A);
    pt_C = REAL(C);
    for(i = 0; i < pt_CLDIM[0] * pt_CLDIM[1]; i++){
        *pt_C = *pt_A;
        pt_A++;
        pt_C++;
    }
    
    // Call Fortran.
    F77_CALL(pdpotrf)(CHARPT(UPLO, 0), INTEGER(N),
        REAL(C), &IJ, &IJ, INTEGER(DESCA), INTEGER(INFO));
    
    // Return. 
    UNPROTECT(4);
                return(RET);
}



// The beast
SEXP R_PDSYEVX(SEXP JOBZ, SEXP RANGE, SEXP N, SEXP A, SEXP DESCA, SEXP VL, SEXP VU, SEXP IL, SEXP IU, SEXP ABSTOL, SEXP ORFAC)
{
    char uplo = 'U';
    
    const int IJ = 1;
    int i, j;
    int m, nz;
    int lwork, liwork, info;
    int descz[9], ldm[2], blacs[5];
    int tmp_liwork;
    int ownany;
    int unpt;
    int *iwork, *ifail, *iclustr;
    
    double tmp_lwork;
    double *work;
    double *w, *z, *gap;
    double *a;
    
    SEXP RET, RET_NAMES, W, Z, IFAIL, M;
    
    
    // grid and local information
    pdims_(INTEGER(DESCA), ldm, blacs);
    if (ldm[0] < 1 || ldm[1] < 1)
        ownany = FALSE;
    else
        ownany = TRUE;
    
    ldm[0] = nonzero(ldm[0]);
    ldm[1] = nonzero(ldm[1]);
    
    
    // Setup for the setup
    for (i=0; i<9; i++)
        descz[i] = INTEGER(DESCA)[i];
    
    w = R_alloc(INTEGER(N)[0], sizeof(double));
    z = R_alloc(ldm[0]*ldm[1], sizeof(double));
    gap = R_alloc(blacs[1]*blacs[2], sizeof(double));
    
    PROTECT(A);
    a = R_alloc(ldm[0]*ldm[1], sizeof(double));
    for (i=0; i<ldm[0]*ldm[1]; i++)
        a[i] = REAL(A)[i];
    
    ifail = R_alloc(INTEGER(N)[0], sizeof(int));
    iclustr = R_alloc(2*blacs[1]*blacs[2], sizeof(int));
    
    
    // Allocate local workspace
    lwork = -1;
    liwork = -1;
    info = 0;
    
    pdsyevx_(CHARPT(JOBZ, 0), CHARPT(RANGE, 0), &uplo, 
        INTEGER(N), a, &IJ, &IJ, INTEGER(DESCA), 
        REAL(VL), REAL(VU), INTEGER(IL), INTEGER(IU), 
        REAL(ABSTOL), &m, &nz, w, 
        REAL(ORFAC), z, &IJ, &IJ, descz, 
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
        INTEGER(N), a, &IJ, &IJ, INTEGER(DESCA), 
        REAL(VL), REAL(VU), INTEGER(IL), INTEGER(IU), 
        REAL(ABSTOL), &m, &nz, w, 
        REAL(ORFAC), z, &IJ, &IJ, descz, 
        work, &lwork, iwork, &liwork, 
        ifail, iclustr, gap, &info);
    
    
    PROTECT(W = allocVector(REALSXP, m));
    for (i=0; i<m; i++)
        REAL(W)[i] = w[i];
    
    
/*    PROTECT(IFAIL = allocVector(INTSXP, m));*/
/*    for (i=0; i<m; i++)*/
/*        INTEGER(IFAIL)[0] = ifail[i];*/
    
    
    // Manage the return
    if (CHARPT(JOBZ, 0)[0] == 'N') // Only eigenvalues are computed
    {
        PROTECT(RET = allocVector(VECSXP, 1));
        PROTECT(RET_NAMES = allocVector(STRSXP, 1));
        
        SET_VECTOR_ELT(RET, 0, W);
        
        SET_STRING_ELT(RET_NAMES, 0, mkChar("values")); 
        
        setAttrib(RET, R_NamesSymbol, RET_NAMES);
        
        unpt = 4;
    }
    else // eigenvalues + eigenvectors
    {
        PROTECT(Z = allocMatrix(REALSXP, ldm[0], ldm[1]));
        for (i=0; i<ldm[0]*ldm[1]; i++)
            REAL(Z)[i] = z[i];
        
        PROTECT(M = allocVector(INTSXP, 1));
        INTEGER(M)[0] = m;
        
        PROTECT(RET = allocVector(VECSXP, 3));
        PROTECT(RET_NAMES = allocVector(STRSXP, 3));
        
        SET_VECTOR_ELT(RET, 0, W);
        SET_VECTOR_ELT(RET, 1, Z);
        SET_VECTOR_ELT(RET, 2, M);
        
        SET_STRING_ELT(RET_NAMES, 0, mkChar("values")); 
        SET_STRING_ELT(RET_NAMES, 1, mkChar("vectors")); 
        SET_STRING_ELT(RET_NAMES, 2, mkChar("m")); 
        
        setAttrib(RET, R_NamesSymbol, RET_NAMES);
        
        unpt = 6;
    }
    
    
    UNPROTECT(unpt);
    return RET;
}


// -------------------------------------------------------- 
// Auxillary 
// -------------------------------------------------------- 


// Matrix norms
SEXP R_PDLANGE(SEXP TYPE, SEXP M, SEXP N, SEXP A, SEXP DESCA)
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


