#include <R.h>
#include <Rinternals.h>
#include "base_global.h"

/*SEXP R_mpi_blacs_initialize(SEXP NPROW, SEXP NPCOL, SEXP ICTXT, SEXP MYROW,*/
/*    SEXP MYCOL)*/
/*{*/
/*  F77_CALL(mpi_blacs_initialize)(INTEGER(NPROW), INTEGER(NPCOL),*/
/*      INTEGER(ICTXT), INTEGER(MYROW), INTEGER(MYCOL));*/
/*  return(R_NilValue);*/
/*}*/

SEXP R_dgsum2d(SEXP ICTXT, SEXP SCOPE, SEXP M, SEXP A, SEXP LDA)
{
  int i;
/*  double *pt_A = REAL(A), *pt_OUT;*/
  const int m = INTEGER(M)[0];
  char top = ' ';
  const int dest = -1;
  
  SEXP OUT;
  PROTECT(OUT = allocVector(REALSXP, m));
  
  memcpy(REAL(OUT), REAL(A), m*sizeof(double));
  
/*  pt_OUT = REAL(OUT);*/
/*  for(i = 0; i < INTEGER(M)[0]; i++){*/
/*    *pt_OUT = *pt_A;*/
/*    pt_OUT++;*/
/*    pt_A++;*/
/*  }*/
  
/*  Cdgsum2d(int ConTxt, char *scope, char *top, int m, int n, double *A,*/
/*    int lda, int rdest, int cdest)*/
  
  Cdgsum2d(INTEGER(ICTXT)[0], CHARPT(SCOPE, 0), &top, m, 
    1, REAL(OUT), INTEGER(LDA)[0], dest, dest);
  
  UNPROTECT(1);
  return(OUT);
}



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


