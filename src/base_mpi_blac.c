/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

// Copyright 2013, Schmidt

#include "base_global.h"
#include <SEXPtools.h>


SEXP R_optimal_grid(SEXP NPROCS)
{
  R_INIT;
  SEXP NPROW, NPCOL, RET, RET_NAMES;
  
  newRvec(NPROW, 1, "int");
  newRvec(NPCOL, 1, "int");
  
  optimalgrid_(INTP(NPROCS), INTP(NPROW), INTP(NPCOL));
  
  RET_NAMES = make_list_names(2, "nprow", "npcol");
  RET = make_list(RET_NAMES, NPROW, NPCOL);
  
  R_END;
  return RET;
}



SEXP R_blacs_init(SEXP NPROW, SEXP NPCOL, SEXP ICTXT)
{
  R_INIT;
  SEXP MYROW, MYCOL, RET, RET_NAMES;
  
  newRvec(MYROW, 1, "int");
  newRvec(MYCOL, 1, "int");
  
  
  sl_init_(INTEGER(ICTXT), INTEGER(NPROW), INTEGER(NPCOL));
  blacs_gridinfo_(INTEGER(ICTXT), INTEGER(NPROW), INTEGER(NPCOL), 
      INTEGER(MYROW), INTEGER(MYCOL));
  
  
  RET_NAMES = make_list_names(5, "NPROW", "NPCOL", "ICTXT", "MYROW", "MYCOL");
  RET = make_list(RET_NAMES, NPROW, NPCOL, ICTXT, MYROW, MYCOL);
  
  R_END;
  return(RET);
}



/* Reductions */
SEXP R_igsum2d1(SEXP ICTXT, SEXP SCOPE, SEXP M, SEXP N, SEXP A, SEXP LDA, SEXP RDEST, SEXP CDEST)
{
  R_INIT;
  int i;
  const int m = INT(M, 0), n = INT(N, 0);
  char top = ' ';
  
  SEXP OUT;
  newRmat(OUT, m, n, "int");
  
  memcpy(INTEGER(OUT), INTEGER(A), m*n*sizeof(int));
  
  Cigsum2d(INTEGER(ICTXT)[0], CHARPT(SCOPE, 0), &top, m, n, INTEGER(OUT), 
      INTEGER(LDA)[0], INTEGER(RDEST)[0], INTEGER(CDEST)[0]);
  
  R_END;
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


