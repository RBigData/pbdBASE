/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

// Copyright 2013, Schmidt

#include "pbdBASE.h"


/* Computing QR */
SEXP R_PDGEQPF(SEXP TOL, SEXP M, SEXP N,
    SEXP A, SEXP DESCA)
{
  R_INIT;
  int lwork = -1;
  int IJ = 1;
  double work = 0.0;
  double tmp = 0.0;
  double *p_work;
  const int ltau = MIN(INT(M, 0), INT(N, 0));
  SEXP RET, RET_NAMES, INFO, A_OUT, TAU, IPIV, RANK;
  
  newRvec(INFO, 1, "int", true);
  newRmat(A_OUT, nrows(A), ncols(A), "dbl");
  newRvec(TAU, ltau, "dbl");
  newRvec(IPIV, ncols(A), "int");
  newRvec(RANK, 1, "int");
  
  
  /* Copy A since pdorgqr writes in place */
  memcpy(DBLP(A_OUT), DBLP(A), nrows(A)*ncols(A)*sizeof(double));
  
  /* workspace query */
  rpdgeqpf_(DBLP(TOL), INTP(M), INTP(N),
      &tmp, &IJ, &IJ, INTP(DESCA),
      &IJ, &tmp,
      &work, &lwork, &IJ, INTP(INFO));
  
  /* allocate work vector and factor A=QR */
  lwork = (int) work;
  lwork = nonzero(lwork);
  p_work = (double *) R_alloc(lwork, sizeof(double));
  
  rpdgeqpf_(DBLP(TOL), INTP(M), INTP(N),
      DBLP(A_OUT), &IJ, &IJ, INTP(DESCA),
      INTP(IPIV), DBLP(TAU),
      p_work, &lwork, INTP(RANK), INTP(INFO));
  
  
  // Manage return
  RET_NAMES = make_list_names(5, "qr", "rank", "tau", "pivot", "INFO");
  RET = make_list(RET_NAMES, 5, A_OUT, RANK, TAU, IPIV, INFO);
  
  R_END;
  return RET;
}



/* For computing Q*y or Q^T*y */
SEXP R_PDORMQR(SEXP SIDE, SEXP TRANS, SEXP M, SEXP N, SEXP K,
    SEXP A, SEXP ALDIM, SEXP DESCA,
    SEXP TAU,
    SEXP B, SEXP BLDIM, SEXP DESCB)
{
  R_INIT;
  int i, *pt_ALDIM = INTEGER(ALDIM), *pt_BLDIM = INTEGER(BLDIM);
  int lwork = -1;
  int IJ = 1;
  double *pt_ORG, *pt_COPY, *A_CPY;
  double work = 0.0;
  double tmp = 0.0;
  double *p_work;
  SEXP RET, RET_NAMES, INFO, B_OUT;
  
  /* Protect R objects. */
  newRvec(INFO, 1, "int", true);
  newRmat(B_OUT, pt_BLDIM[0], pt_BLDIM[1], "dbl");
  
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
  pdormqr_(STR(SIDE, 0), STR(TRANS, 0),
      INTEGER(M), INTEGER(N), INTEGER(K),
      &tmp, &IJ, &IJ, INTEGER(DESCA),
      &tmp,
      &tmp, &IJ, &IJ, INTEGER(DESCB),
      &work, &lwork, INTEGER(INFO));
  
  /* allocate work vector and compute Q*y or Q^T*y */
  lwork = (int) work;
  lwork = nonzero(lwork);
  p_work = (double *) R_alloc(lwork, sizeof(double));
  
  pdormqr_(STR(SIDE, 0), STR(TRANS, 0),
      INTEGER(M), INTEGER(N), INTEGER(K),
      A_CPY, &IJ, &IJ, INTEGER(DESCA),
      REAL(TAU),
      REAL(B_OUT), &IJ, &IJ, INTEGER(DESCB),
      p_work, &lwork, INTEGER(INFO));
  
  /* Return. */
  RET_NAMES = make_list_names(2, "INFO", "B");
  RET = make_list(RET_NAMES, 2, INFO, B_OUT);
  
  R_END;
  return RET;
}



/* recovering Q from a QR */
SEXP R_PDORGQR(SEXP M, SEXP N, SEXP K, SEXP A, SEXP ALDIM, SEXP DESCA, SEXP TAU)
{
  R_INIT;
  int i, *pt_ALDIM = INTEGER(ALDIM);
  int lwork = -1;
  int IJ = 1;
  double *pt_ORG, *pt_COPY;
  double work = 0.0;
  double tmp = 0.0;
  double *p_work;
  SEXP RET, RET_NAMES, INFO, A_OUT;
  
  /* Protect R objects. */
  newRvec(INFO, 1, "int", true);
  newRmat(A_OUT, pt_ALDIM[0], pt_ALDIM[1], "dbl");
  
  /* Copy A since pdorgqr writes in place */
  pt_ORG = REAL(A);
  pt_COPY = REAL(A_OUT);
  for(i = 0; i < pt_ALDIM[0] * pt_ALDIM[1]; i++){
      *pt_COPY = *pt_ORG;
      pt_ORG++;
      pt_COPY++;
  }
  
  /* workspace query */
  pdorgqr_(INTEGER(M), INTEGER(N), INTEGER(K),
      &tmp, &IJ, &IJ, INTEGER(DESCA),
      &tmp,
      &work, &lwork, INTEGER(INFO));
  
  /* allocate work vector and recover Q */
  lwork = (int) work;
  lwork = nonzero(lwork);
  p_work = (double *) R_alloc(lwork, sizeof(double));
  
  pdorgqr_(INTEGER(M), INTEGER(N), INTEGER(K),
      REAL(A_OUT), &IJ, &IJ, INTEGER(DESCA),
      REAL(TAU),
      p_work, &lwork, INTEGER(INFO));
  
  /* Return. */
  RET_NAMES = make_list_names(2, "INFO", "A");
  RET = make_list(RET_NAMES, 2, INFO, A_OUT);
  
  R_END;
  return RET;
}

