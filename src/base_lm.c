/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

// Copyright 2013, Schmidt

#include "pbdBASE.h"


/* For computing LLS solution, either over or    under-determined. */
/* In the case that A is rank deficient, the 'limited pivoting    */
/* strategy from R's dqrls.f is used. I don't think this is         */
/* numerically stable, but it's the cost of preserving the            */
/* order of the model matrix, which has important interpretive    */
/* value sometimes.                                                                                         */
SEXP R_PDGELS(SEXP TOL, SEXP M, SEXP N, SEXP NRHS,
    SEXP A, SEXP DESCA, SEXP B, SEXP DESCB, SEXP LTAU)
{
  R_INIT;
  int lwork = -1;
  int IJ = 1;
  int i;
  double *pt_ORG, *pt_COPY, *pt_EFF, *pt_FT, *pt_RSD, *p_work;
  double tmp = 0.0;
  double work = 0.0;
  int NN = INT(N);
  int DESCA_CP[9];
  
  for (i=0; i<9; i++)
    DESCA_CP[i] = INT(DESCA, i);
  
  char trans = 'N'; // If trans='T', expect all hell to break loose
  
  SEXP RET, RET_NAMES, INFO, A_OUT, B_OUT, EFF, FT, RSD, TAU, IPIV, RANK;
  
  /* set up return */
  newRvec(INFO, 1, "int", true);
  newRmat(A_OUT, nrows(A), ncols(A), "dbl");
  newRmat(B_OUT, nrows(B), ncols(B), "dbl");
  newRmat(EFF, nrows(B), ncols(B), "dbl");
  newRmat(FT, nrows(B), ncols(B), "dbl");
  newRmat(RSD, nrows(B), ncols(B), "dbl");
  newRvec(TAU, INT(LTAU, 0), "dbl");
  newRvec(IPIV, ncols(A), "int");
  newRvec(RANK, 1, "int");
  
  
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
  rpdgels_(REAL(TOL), &trans,
    INTP(M), &NN, INTP(NRHS),
    &tmp, &IJ, &IJ, DESCA_CP,
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
    INTP(M), &NN, INTP(NRHS),
    REAL(A_OUT), &IJ, &IJ, DESCA_CP,
    REAL(B_OUT), &IJ, &IJ, INTP(DESCB),
    REAL(EFF), REAL(FT), REAL(RSD),
    REAL(TAU), p_work, &lwork,
    INTP(IPIV), INTP(RANK), INTP(INFO));
  
  
  // Manage return
  RET_NAMES = make_list_names(9, "INFO", "A", "B", "EFF", "FT", "RSD", "TAU", "IPIV", "RANK");
  RET = make_list(RET_NAMES, 9, INFO, A_OUT, B_OUT, EFF, FT, RSD, TAU, IPIV, RANK);
  
  R_END;
  return RET;
}
