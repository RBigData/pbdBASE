#include <R.h>
#include <Rinternals.h>
#include "base_global.h"

SEXP R_mpi_blacs_initialize(SEXP NPROW, SEXP NPCOL, SEXP ICTXT, SEXP MYROW,
    SEXP MYCOL)
{
  F77_CALL(mpi_blacs_initialize)(INTEGER(NPROW), INTEGER(NPCOL),
      INTEGER(ICTXT), INTEGER(MYROW), INTEGER(MYCOL));
  return(R_NilValue);
}

SEXP R_dgsum2d(SEXP ICTXT, SEXP SCOPE, SEXP M, SEXP N, SEXP A, SEXP LDA)
{
  int i;
  double *pt_A = REAL(A), *pt_OUT;
  char top = ' ';
  const int dest = -1;
  SEXP OUT;

  PROTECT(OUT = allocMatrix(REALSXP, INTEGER(M)[0], INTEGER(N)[0]));
  pt_OUT = REAL(OUT);
  for(i = 0; i < INTEGER(M)[0] * INTEGER(N)[0]; i++){
    *pt_OUT = *pt_A;
    pt_OUT++;
    pt_A++;
  }
  
/*  Cdgsum2d(int ConTxt, char *scope, char *top, int m, int n, double *A,*/
/*    int lda, int rdest, int cdest)*/
  
  Cdgsum2d(INTEGER(ICTXT)[0], CHARPT(SCOPE, 0), &top, INTEGER(M)[0], 
    INTEGER(N)[0], REAL(OUT), INTEGER(LDA)[0], dest, dest);
  
  UNPROTECT(1);
  return(OUT);
}




/*SEXP R_dims(SEXP DESC)*/
/*{*/
/*  SEXP LDM, BLACS;*/
/*  PROTECT(LDM = allocVector(INTSXP, 2));*/
/*  PROTECT(BLACS = allocVector(INTSXP, 4));*/
/*  */
/*  dims_(INTEGER(DESC), INTEGER(LDM), INTEGER(BLACS));*/
/*  */
/*  Rprintf("%d %d\n", INTEGER(LDM)[0], INTEGER(LDM)[1]);*/
/*  */
/*  UNPROTECT(2);*/
/*  return(R_NilValue);*/
/*}*/




