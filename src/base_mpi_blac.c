#include <R.h>
#include <Rinternals.h>
#include "base_global.h"


SEXP R_optimal_grid(SEXP NPROCS)
{
    SEXP NPROW, NPCOL, RET, RET_NAMES;
    
    PROTECT(RET = allocVector(VECSXP, 5));
    PROTECT(RET_NAMES = allocVector(STRSXP, 5));
    
    PROTECT(NPROW = allocVector(INTSXP, 1));
    PROTECT(NPCOL = allocVector(INTSXP, 1));
    
    SET_VECTOR_ELT(RET, 0, NPROW);
    SET_VECTOR_ELT(RET, 1, NPCOL);
    
    SET_STRING_ELT(RET_NAMES, 0, mkChar("nprow")); 
    SET_STRING_ELT(RET_NAMES, 1, mkChar("npcol")); 
    
    setAttrib(RET, R_NamesSymbol, RET_NAMES);
    
    optimalgrid_(INTEGER(NPROCS), INTEGER(NPROW), INTEGER(NPCOL));
    
    UNPROTECT(4);
    return(RET);
}

SEXP R_blacs_init(SEXP NPROW, SEXP NPCOL, SEXP ICTXT)
{
    SEXP MYROW, MYCOL, RET, RET_NAMES;
    
    PROTECT(RET = allocVector(VECSXP, 5));
    PROTECT(RET_NAMES = allocVector(STRSXP, 5));
    
    PROTECT(MYROW = allocVector(INTSXP, 1));
    PROTECT(MYCOL = allocVector(INTSXP, 1));
    
    SET_VECTOR_ELT(RET, 0, NPROW);
    SET_VECTOR_ELT(RET, 1, NPCOL);
    SET_VECTOR_ELT(RET, 2, ICTXT);
    SET_VECTOR_ELT(RET, 3, MYROW);
    SET_VECTOR_ELT(RET, 4, MYCOL);
    
    SET_STRING_ELT(RET_NAMES, 0, mkChar("NPROW")); 
    SET_STRING_ELT(RET_NAMES, 1, mkChar("NPCOL")); 
    SET_STRING_ELT(RET_NAMES, 2, mkChar("ICTXT")); 
    SET_STRING_ELT(RET_NAMES, 3, mkChar("MYROW")); 
    SET_STRING_ELT(RET_NAMES, 4, mkChar("MYCOL")); 
    
    setAttrib(RET, R_NamesSymbol, RET_NAMES);
    
    /* initialize */
    sl_init_(INTEGER(ICTXT), INTEGER(NPROW), INTEGER(NPCOL));
    blacs_gridinfo_(INTEGER(ICTXT), INTEGER(NPROW), INTEGER(NPCOL), 
        INTEGER(MYROW), INTEGER(MYCOL));
    
    UNPROTECT(4);
    return(RET);
}


/* Reductions */

SEXP R_igsum2d1(SEXP ICTXT, SEXP SCOPE, SEXP M, SEXP N, SEXP A, SEXP LDA, SEXP RDEST, SEXP CDEST)
{
    int i;
    const int m = INTEGER(M)[0], n = INTEGER(N)[0];
    char top = ' ';
    
    SEXP OUT;
    PROTECT(OUT = allocMatrix(INTSXP, m, n));
    
    memcpy(INTEGER(OUT), INTEGER(A), m*n*sizeof(int));
    
    Cigsum2d(INTEGER(ICTXT)[0], CHARPT(SCOPE, 0), &top, m, n, INTEGER(OUT), 
        INTEGER(LDA)[0], INTEGER(RDEST)[0], INTEGER(CDEST)[0]);
    
    UNPROTECT(1);
    return(OUT);
}

SEXP R_dgsum2d1(SEXP ICTXT, SEXP SCOPE, SEXP M, SEXP N, SEXP A, SEXP LDA, SEXP RDEST, SEXP CDEST)
{
    int i;
    const int m = INTEGER(M)[0], n = INTEGER(N)[0];
    char top = ' ';
    
    SEXP OUT;
    PROTECT(OUT = allocMatrix(REALSXP, m, n));
    
    memcpy(REAL(OUT), REAL(A), m*n*sizeof(double));
    
    Cdgsum2d(INTEGER(ICTXT)[0], CHARPT(SCOPE, 0), &top, m, n, REAL(OUT), 
        INTEGER(LDA)[0], INTEGER(RDEST)[0], INTEGER(CDEST)[0]);
    
    UNPROTECT(1);
    return(OUT);
}


SEXP R_igamx2d1(SEXP ICTXT, SEXP SCOPE, SEXP M, SEXP N, SEXP A, SEXP LDA, SEXP RDEST, SEXP CDEST)
{
    int i;
    const int m = INTEGER(M)[0], n = INTEGER(N)[0];
    char top = ' ';
    const int rcflag = -1;
    
    SEXP OUT;
    PROTECT(OUT = allocMatrix(INTSXP, m, n));
    
    memcpy(INTEGER(OUT), INTEGER(A), m*n*sizeof(int));
    
    Cdgamx2d(INTEGER(ICTXT)[0], CHARPT(SCOPE, 0), &top, m, n, INTEGER(OUT), 
        INTEGER(LDA)[0], rcflag, rcflag, rcflag, INTEGER(RDEST)[0], INTEGER(CDEST)[0]);
    
    UNPROTECT(1);
    return(OUT);
}

SEXP R_dgamx2d1(SEXP ICTXT, SEXP SCOPE, SEXP M, SEXP N, SEXP A, SEXP LDA, SEXP RDEST, SEXP CDEST)
{
    int i;
    const int m = INTEGER(M)[0], n = INTEGER(N)[0];
    char top = ' ';
    const int rcflag = -1;
    
    SEXP OUT;
    PROTECT(OUT = allocMatrix(REALSXP, m, n));
    
    memcpy(REAL(OUT), REAL(A), m*n*sizeof(double));
    
    Cdgamx2d(INTEGER(ICTXT)[0], CHARPT(SCOPE, 0), &top, m, n, REAL(OUT), 
        INTEGER(LDA)[0], rcflag, rcflag, rcflag, INTEGER(RDEST)[0], INTEGER(CDEST)[0]);
    
    UNPROTECT(1);
    return(OUT);
}


SEXP R_igamn2d1(SEXP ICTXT, SEXP SCOPE, SEXP M, SEXP N, SEXP A, SEXP LDA, SEXP RDEST, SEXP CDEST)
{
    int i;
    const int m = INTEGER(M)[0], n = INTEGER(N)[0];
    char top = ' ';
    const int rcflag = -1;
    
    SEXP OUT;
    PROTECT(OUT = allocMatrix(INTSXP, m, n));
    
    memcpy(INTEGER(OUT), INTEGER(A), m*n*sizeof(int));
    
    Cdgamn2d(INTEGER(ICTXT)[0], CHARPT(SCOPE, 0), &top, m, n, INTEGER(OUT), 
        INTEGER(LDA)[0], rcflag, rcflag, rcflag, INTEGER(RDEST)[0], INTEGER(CDEST)[0]);
    
    UNPROTECT(1);
    return(OUT);
}


SEXP R_dgamn2d1(SEXP ICTXT, SEXP SCOPE, SEXP M, SEXP N, SEXP A, SEXP LDA, SEXP RDEST, SEXP CDEST)
{
    int i;
    const int m = INTEGER(M)[0], n = INTEGER(N)[0];
    char top = ' ';
    const int rcflag = -1;
    
    SEXP OUT;
    PROTECT(OUT = allocMatrix(REALSXP, m, n));
    
    memcpy(REAL(OUT), REAL(A), m*n*sizeof(double));
    
    Cdgamn2d(INTEGER(ICTXT)[0], CHARPT(SCOPE, 0), &top, m, n, REAL(OUT), 
        INTEGER(LDA)[0], rcflag, rcflag, rcflag, INTEGER(RDEST)[0], INTEGER(CDEST)[0]);
    
    UNPROTECT(1);
    return(OUT);
}

// Point to point send/receive
SEXP R_dgesd2d1(SEXP ICTXT, SEXP M, SEXP N, SEXP A, SEXP LDA, SEXP RDEST, SEXP CDEST)
{
    int i;
    const int m = INTEGER(M)[0], n = INTEGER(N)[0];
    
    SEXP OUT;
    PROTECT(OUT = allocMatrix(REALSXP, m, n));
    
    memcpy(REAL(OUT), REAL(A), m*n*sizeof(double));
    
    Cdgesd2d(INTEGER(ICTXT)[0], m, n, REAL(OUT), INTEGER(LDA)[0], 
        INTEGER(RDEST)[0], INTEGER(CDEST)[0]);
    
    UNPROTECT(1);
    return(OUT);
}

SEXP R_dgerv2d1(SEXP ICTXT, SEXP M, SEXP N, SEXP A, SEXP LDA, SEXP RDEST, SEXP CDEST)
{
    int i;
    const int m = INTEGER(M)[0], n = INTEGER(N)[0];
    
    SEXP OUT;
    PROTECT(OUT = allocMatrix(REALSXP, m, n));
    
    memcpy(REAL(OUT), REAL(A), m*n*sizeof(double));
    
    Cdgerv2d(INTEGER(ICTXT)[0], m, n, REAL(OUT), INTEGER(LDA)[0], 
        INTEGER(RDEST)[0], INTEGER(CDEST)[0]);
    
    UNPROTECT(1);
    return(OUT);
}


