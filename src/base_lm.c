/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

// Copyright 2013, Schmidt

#include <R.h>
#include <Rinternals.h>

#include "base_global.h"
#include "Rtools/Rtools.h"


/* For computing LLS solution, either over or    under-determined. */
/* In the case that A is rank deficient, the 'limited pivoting    */
/* strategy from R's dqrls.f is used. I don't think this is         */
/* numerically stable, but it's the cost of preserving the            */
/* order of the model matrix, which has important interpretive    */
/* value sometimes.                                                                                         */
SEXP R_PDGELS(SEXP TOL, SEXP M, SEXP N, SEXP NRHS,
    SEXP A, SEXP DESCA, SEXP B, SEXP DESCB,
    SEXP LTAU)
{
    int ptct = 0;
    int lwork = -1;
    const int IJ = 1;
    int i;
    double *pt_ORG, *pt_COPY, *pt_EFF, *pt_FT, *pt_RSD, *p_work;
    const double tmp = 0.0;
    double work = 0.0;
    
    char trans = 'N'; // If trans='T', expect all hell to break loose
    
    SEXP RET, RET_NAMES, INFO, A_OUT, B_OUT, EFF, FT, RSD, TAU, IPIV, RANK;
    
    /* set up return */
    PT(RET = allocVector(VECSXP, 9), ptct);
    PT(RET_NAMES = allocVector(STRSXP, 9), ptct);
    
    PT(INFO = allocVector(INTSXP, 1), ptct);
    PT(A_OUT = allocMatrix(REALSXP, nrows(A), ncols(A)), ptct);
    PT(B_OUT = allocMatrix(REALSXP, nrows(B), ncols(B)), ptct);
    PT(EFF = allocMatrix(REALSXP, nrows(B), ncols(B)), ptct);
    PT(FT = allocMatrix(REALSXP, nrows(B), ncols(B)), ptct);
    PT(RSD = allocMatrix(REALSXP, nrows(B), ncols(B)), ptct);
    PT(TAU = allocVector(REALSXP, INT(LTAU, 0)), ptct);
    PT(IPIV = allocVector(INTSXP, ncols(A)), ptct);
    PT(RANK = allocVector(INTSXP, 1), ptct);
    
    
    /* Copy A and B since pdgels writes in place, also initialize */
    memcpy(DBLP(A_OUT), DBLP(A), nrows(A)*ncols(A)*sizeof(double));
    
    /* Set FT, RSD, and EFF to 0 */
    pt_ORG = REAL(B);
    pt_COPY = REAL(B_OUT);
    pt_EFF = REAL(EFF);
    pt_FT = REAL(FT);
    pt_RSD = REAL(RSD);
    
    for(i = 0; i < nrows(B)*ncols(B); i++){
        *pt_COPY = *pt_ORG;
        *pt_FT = 0.0;
        *pt_RSD = 0.0;
        
        pt_EFF++;
        pt_FT++;
        pt_RSD++;
        
        pt_ORG++;
        pt_COPY++;
    }
    
    
    /* workspace query */
    INT(INFO, 0) = 0;
    
    rpdgels_(REAL(TOL), &trans,
        INTP(M), INTP(N), INTP(NRHS),
        &tmp, &IJ, &IJ, INTP(DESCA),
        &tmp, &IJ, &IJ, INTP(DESCB),
        &tmp, &tmp, &tmp,
        &tmp, &work, &lwork,
        &IJ, &IJ, INTP(INFO));
    
    
    /* allocate work vector */
    lwork = (int) work;
    lwork = nonzero(lwork);
    p_work = (double *) R_alloc(lwork, sizeof(double));
    
    
    
    /*    and compute LLS solution */
    rpdgels_(REAL(TOL), &trans,
        INTP(M), INTP(N), INTP(NRHS),
        REAL(A_OUT), &IJ, &IJ, INTP(DESCA),
        REAL(B_OUT), &IJ, &IJ, INTP(DESCB),
        REAL(EFF), REAL(FT), REAL(RSD),
        REAL(TAU), p_work, &lwork,
        INTP(IPIV), INTP(RANK), INTP(INFO));
    
    
    // Manage return
    RET_NAMES = make_list_names(9, "INFO", "A", "B", "EFF", "FT", "RSD", "TAU", "IPIV", "RANK");
    RET = make_list(RET_NAMES, INFO, A_OUT, B_OUT, EFF, FT, RSD, TAU, IPIV, RANK);
    
    UNPT(ptct);
    return(RET);
}


/* ----------------------------------------------------- */
/*                     QR functions no one will ever use                     */
/* ----------------------------------------------------- */

/* Computing QR */
SEXP R_PDGEQPF(SEXP TOL, SEXP M, SEXP N,
    SEXP A, SEXP ALDIM, SEXP DESCA,
    SEXP LTAU)
{
    int i, *pt_ALDIM = INTEGER(ALDIM);
    int lwork = -1;
    const int IJ = 1;

    double *pt_ORG, *pt_COPY;
    double work = 0.0;
    const double tmp = 0.0;

    double *p_work;

    SEXP RET, RET_NAMES, INFO, A_OUT, TAU, IPIV, RANK;

    /* Protect R objects. */
    PROTECT(INFO = allocVector(INTSXP, 1));
    PROTECT(A_OUT = allocMatrix(REALSXP, pt_ALDIM[0], pt_ALDIM[1]));
    PROTECT(TAU = allocVector(REALSXP, INTEGER(LTAU)[0]));
    PROTECT(IPIV = allocVector(INTSXP, pt_ALDIM[1]));
    PROTECT(RANK = allocVector(INTSXP, 1));

    /* Manage return */
    PROTECT(RET = allocVector(VECSXP, 5));
    PROTECT(RET_NAMES = allocVector(STRSXP, 5));

    SET_VECTOR_ELT(RET, 0, A_OUT);
    SET_VECTOR_ELT(RET, 1, RANK);
    SET_VECTOR_ELT(RET, 2, TAU);
    SET_VECTOR_ELT(RET, 3, IPIV);
    SET_VECTOR_ELT(RET, 4, INFO);

    SET_STRING_ELT(RET_NAMES, 0, mkChar("qr"));
    SET_STRING_ELT(RET_NAMES, 1, mkChar("rank"));
    SET_STRING_ELT(RET_NAMES, 2, mkChar("tau"));
    SET_STRING_ELT(RET_NAMES, 3, mkChar("pivot"));
    SET_STRING_ELT(RET_NAMES, 4, mkChar("INFO"));

    setAttrib(RET, R_NamesSymbol, RET_NAMES);

    /* Copy A since pdorgqr writes in place */
    pt_ORG = REAL(A);
    pt_COPY = REAL(A_OUT);
    for(i = 0; i < pt_ALDIM[0] * pt_ALDIM[1]; i++){
        *pt_COPY = *pt_ORG;
        pt_ORG++;
        pt_COPY++;
    }

    /* workspace query */
    INTEGER(INFO)[0] = 0;
    rpdgeqpf_(REAL(TOL), INTEGER(M), INTEGER(N),
        &tmp, &IJ, &IJ, INTEGER(DESCA),
        &IJ, &tmp,
        &work, &lwork, &IJ, INTEGER(INFO));

    /* allocate work vector and factor A=QR */
    lwork = (int) work;
    lwork = nonzero(lwork);
    p_work = (double *) R_alloc(lwork, sizeof(double));

    INTEGER(INFO)[0] = 0;
    rpdgeqpf_(REAL(TOL), INTEGER(M), INTEGER(N),
        REAL(A_OUT), &IJ, &IJ, INTEGER(DESCA),
        INTEGER(IPIV), REAL(TAU),
        p_work, &lwork, INTEGER(RANK), INTEGER(INFO));

    /* Return. */
    UNPROTECT(7);
    return(RET);
}




/* For computing Q*y or Q^T*y */
SEXP R_PDORMQR(SEXP SIDE, SEXP TRANS, SEXP M, SEXP N, SEXP K,
    SEXP A, SEXP ALDIM, SEXP DESCA,
    SEXP TAU,
    SEXP B, SEXP BLDIM, SEXP DESCB)
{
    int i, *pt_ALDIM = INTEGER(ALDIM), *pt_BLDIM = INTEGER(BLDIM);
    int lwork = -1;
    const int IJ = 1;

    double *pt_ORG, *pt_COPY, *A_CPY;
    double work = 0.0;
    const double tmp = 0.0;

    double *p_work;

    SEXP RET, RET_NAMES, INFO, B_OUT;

    /* Protect R objects. */
    PROTECT(RET = allocVector(VECSXP, 2));
    PROTECT(RET_NAMES = allocVector(STRSXP, 2));
    PROTECT(INFO = allocVector(INTSXP, 1));
    PROTECT(B_OUT = allocMatrix(REALSXP, pt_BLDIM[0], pt_BLDIM[1]));

    SET_VECTOR_ELT(RET, 0, INFO);
    SET_VECTOR_ELT(RET, 1, B_OUT);
    SET_STRING_ELT(RET_NAMES, 0, mkChar("INFO"));
    SET_STRING_ELT(RET_NAMES, 1, mkChar("B"));
    setAttrib(RET, R_NamesSymbol, RET_NAMES);

    /* Copy A and B since pdormqr writes in place */
    A_CPY = (double *) R_alloc(pt_ALDIM[0] * pt_ALDIM[1], sizeof(double));
    pt_ORG = REAL(A);
    pt_COPY = A_CPY;
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

    /* workspace query */
    INTEGER(INFO)[0] = 0;
    F77_CALL(pdormqr)(CHARPT(SIDE, 0), CHARPT(TRANS, 0),
        INTEGER(M), INTEGER(N), INTEGER(K),
        &tmp, &IJ, &IJ, INTEGER(DESCA),
        &tmp,
        &tmp, &IJ, &IJ, INTEGER(DESCB),
        &work, &lwork, INTEGER(INFO));

    /* allocate work vector and compute Q*y or Q^T*y */
    lwork = (int) work;
    lwork = nonzero(lwork);
    p_work = (double *) R_alloc(lwork, sizeof(double));

    INTEGER(INFO)[0] = 0;
    F77_CALL(pdormqr)(CHARPT(SIDE, 0), CHARPT(TRANS, 0),
        INTEGER(M), INTEGER(N), INTEGER(K),
        A_CPY, &IJ, &IJ, INTEGER(DESCA),
        REAL(TAU),
        REAL(B_OUT), &IJ, &IJ, INTEGER(DESCB),
        p_work, &lwork, INTEGER(INFO));

    /* Return. */
    UNPROTECT(4);
    return(RET);
}







/* recovering Q from a QR */
SEXP R_PDORGQR(SEXP M, SEXP N, SEXP K, SEXP A, SEXP ALDIM, SEXP DESCA, SEXP TAU)
{
    int i, *pt_ALDIM = INTEGER(ALDIM);
    int lwork = -1;
    const int IJ = 1;

    double *pt_ORG, *pt_COPY;
    double work = 0.0;
    const double tmp = 0.0;

    double *p_work;

    SEXP RET, RET_NAMES, INFO, A_OUT;

    /* Protect R objects. */
    PROTECT(RET = allocVector(VECSXP, 2));
    PROTECT(RET_NAMES = allocVector(STRSXP, 2));
    PROTECT(INFO = allocVector(INTSXP, 1));
    PROTECT(A_OUT = allocMatrix(REALSXP, pt_ALDIM[0], pt_ALDIM[1]));

    SET_VECTOR_ELT(RET, 0, INFO);
    SET_VECTOR_ELT(RET, 1, A_OUT);
    SET_STRING_ELT(RET_NAMES, 0, mkChar("INFO"));
    SET_STRING_ELT(RET_NAMES, 1, mkChar("A"));

    setAttrib(RET, R_NamesSymbol, RET_NAMES);

    /* Copy A since pdorgqr writes in place */
    pt_ORG = REAL(A);
    pt_COPY = REAL(A_OUT);
    for(i = 0; i < pt_ALDIM[0] * pt_ALDIM[1]; i++){
        *pt_COPY = *pt_ORG;
        pt_ORG++;
        pt_COPY++;
    }

    /* workspace query */
    INTEGER(INFO)[0] = 0;
    pdorgqr_(INTEGER(M), INTEGER(N), INTEGER(K),
        &tmp, &IJ, &IJ, INTEGER(DESCA),
        &tmp,
        &work, &lwork, INTEGER(INFO));

    /* allocate work vector and recover Q */
    lwork = (int) work;
    lwork = nonzero(lwork);
    p_work = (double *) R_alloc(lwork, sizeof(double));

    INTEGER(INFO)[0] = 0;
    pdorgqr_(INTEGER(M), INTEGER(N), INTEGER(K),
        REAL(A_OUT), &IJ, &IJ, INTEGER(DESCA),
        REAL(TAU),
        p_work, &lwork, INTEGER(INFO));

    /* Return. */
    UNPROTECT(4);
    return(RET);
}

