#include <R.h>
#include <Rinternals.h>
#include "base_global.h"

SEXP R_PDGESV(SEXP N, SEXP NRHS, SEXP MXLDIMS, 
  SEXP A, SEXP ALDIM, SEXP DESCA,
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
/*  IPIV ( MXLDIMS(1) + DESCA(6) )*/
  int * ipiv;
  ipiv = (int *) R_alloc(INTEGER(MXLDIMS)[0] + INTEGER(DESCA)[5], sizeof(int));
  
  /* Set info */
  INTEGER(INFO)[0] = 0;
  
  F77_CALL(pdgesv)(INTEGER(N), INTEGER(NRHS),
    A_OUT, &IJ, &IJ, INTEGER(DESCA), ipiv,
    REAL(B_OUT), &IJ, &IJ, INTEGER(DESCB), INTEGER(INFO));
  
  /* Return. */
  UNPROTECT(4);
  // Free(A_OUT);    // .Call takes care of this.
  return(RET);
} /* End of R_PDGESV(). */



SEXP R_PDGESVD(SEXP M, SEXP N, SEXP ASIZE, 
                SEXP A, SEXP DESCA, SEXP ALDIM, SEXP ULDIM,
                SEXP DESCU, SEXP VTLDIM, SEXP DESCVT, SEXP JOBU, SEXP JOBVT)
{
  int i, *pt_ALDIM = INTEGER(ALDIM);
  double *A_OUT;
  SEXP RET, RET_NAMES, INFO, D, U, VT;

  /* Extra needed. */
  int temp_IJ = 1, temp_lwork = -1;
  double temp_A = 0, temp_work = 0, *WORK;

  /* Protect R objects. */
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

  /* Make copy of original data, since pdgesvd destroys it */
  i = pt_ALDIM[0] * pt_ALDIM[1];
  A_OUT = (double *) R_alloc(i, sizeof(double));
  memcpy(A_OUT, REAL(A), i * sizeof(double));

  /* Query size of workspace */
    INTEGER(INFO)[0] = 0;
    F77_CALL(pdgesvd)(CHARPT(JOBU, 0), CHARPT(JOBVT, 0),
      INTEGER(M), INTEGER(N),
      &temp_A, &temp_IJ, &temp_IJ, INTEGER(DESCA),
      &temp_A, &temp_A, &temp_IJ, &temp_IJ, INTEGER(DESCU),
      &temp_A, &temp_IJ, &temp_IJ, INTEGER(DESCVT),
      &temp_work, &temp_lwork, INTEGER(INFO));
    temp_lwork = (int) temp_work;

    /* Allocate work vector and calculate svd */
    WORK = (double *) R_alloc(temp_lwork, sizeof(double));
    
    INTEGER(INFO)[0] = 0;
    F77_CALL(pdgesvd)(CHARPT(JOBU, 0), CHARPT(JOBVT, 0),
      INTEGER(M), INTEGER(N),
      A_OUT, &temp_IJ, &temp_IJ, INTEGER(DESCA),
      REAL(D), REAL(U), &temp_IJ, &temp_IJ, INTEGER(DESCU),
      REAL(VT), &temp_IJ, &temp_IJ, INTEGER(DESCVT),
      WORK, &temp_lwork, INTEGER(INFO));

  /* Return. */
  UNPROTECT(6);

  return(RET);
} /* End of R_PDGESVD(). */



SEXP R_PDGETRI(SEXP A, SEXP CLDIM, SEXP DESCA, SEXP N)
{
  int *pt_CLDIM = INTEGER(CLDIM), *ipiv, *iwork;
  int i, lwork, liwork;
  double *pt_A, *pt_C, *work;
  
  const int query = -1, IJ = 1;
  double tmp_A = 0, tmp1 = 0;
  
  /* R objects. */
  SEXP RET, RET_NAMES, INFO, C;
  
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
  
  /* IPIV ( N + DESCA(6) ) */
  ipiv = (int *) R_alloc(INTEGER(N)[0] + INTEGER(DESCA)[5], sizeof(int));
  
  /* LU decomposition */
  INTEGER(INFO)[0] = 0;
  F77_CALL(pdgetrf)(INTEGER(N), INTEGER(N), REAL(C), &IJ, &IJ, 
    INTEGER(DESCA), ipiv, INTEGER(INFO));
  
  /* workspace query for inverse */
  INTEGER(INFO)[0] = 0;
  F77_CALL(pdgetri)(INTEGER(N), &tmp_A, &IJ, &IJ, INTEGER(DESCA),
    &IJ, &tmp1, &query, &liwork, &query, INTEGER(INFO));
  
  /* allocate work arrays and invert A */
  lwork = (int) tmp1;
  work = (double *) R_alloc(lwork, sizeof(double));
  
  iwork = (int *) R_alloc(liwork, sizeof(int));
  
  INTEGER(INFO)[0] = 0;
  F77_CALL(pdgetri)(INTEGER(N), REAL(C), &IJ, &IJ, INTEGER(DESCA),
    ipiv, work, &lwork, iwork, &liwork, INTEGER(INFO));
  
  /* Return. */
  UNPROTECT(4);
        return(RET);
} /* End of R_PDGETRI(). */


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
  
  ipiv = (int *) R_alloc(INTEGER(LIPIV), sizeof(int));
  
  INTEGER(INFO)[0] = 0;
  F77_CALL(pdgetrf)(INTEGER(M), INTEGER(N), REAL(C), 
    &IJ, &IJ, INTEGER(DESCA), ipiv, INTEGER(INFO));
  
  /* Return. */
  UNPROTECT(4);
        return(RET);
} /* End of R_PDGETRF(). */



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
  
  /* Call Fortran. */
  F77_CALL(pdpotrf)(CHARPT(UPLO, 0), INTEGER(N),
    REAL(C), &IJ, &IJ, INTEGER(DESCA), INTEGER(INFO));
  
  /* Return. */
  UNPROTECT(4);
        return(RET);
} /* End of R_PDPOTRF(). */
