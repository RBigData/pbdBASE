/*  Copyright (c) 2012-2013, 2017 by Drew Schmidt <wrathematics@gmail.com>
    All rights reserved.
    
    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:
    
    1. Redistributions of source code must retain the above copyright notice,
    this list of conditions and the following disclaimer.
    
    2. Redistributions in binary form must reproduce the above copyright
    notice, this list of conditions and the following disclaimer in the
    documentation and/or other materials provided with the distribution.
    
    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
    PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
    CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
    PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
    PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
    LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

// We improve the base ScaLAPACK performance as follows. For m>n:
//   * Compute A = QR via pdgeqrf()
//   * Compute SVD(R) = U_R * Sigma * V^T
//   * If needed, U = Q * U_R
// For m<n we can use LQ, but this is not yet complete.

#include <stdint.h>
#include "scalapack.h"

// R.h and Rinternals.h needs to be included after Rconfig.h
#include "pbdBASE.h"
#include <RNACI.h>


typedef int32_t len_t;

typedef struct dvector
{
  len_t len;
  double *data;
} dvector_t;

typedef struct dmatrix
{
  len_t nrows;
  len_t ncols;
  double *restrict data;
} dmatrix_t;

typedef struct ddmatrix
{
  len_t nrows;
  len_t ncols;
  dmatrix_t *restrict data;
  int *restrict desc;
} ddmatrix_t;



typedef struct svdparam
{
  bool inplace;
  bool retu;
  bool retvt;
  bool descending;
} svdparam_t;

typedef struct dsvdstruct
{
  dvector_t *restrict d;
  ddmatrix_t *restrict u;
  ddmatrix_t *restrict vt;
  int info;
} dsvd_t;



static inline int svd_noqr(const char *jobu, const char *jobvt, ddmatrix_t *A, dsvd_t *svd)
{
  // const char jobu = retu ? 'V' : 'N';
  // const char jobvt = retvt ? 'V' : 'N';
  
  double tmp;
  int lwork = -1;
  double *work;
  
#ifdef FC_LEN_T
  pdgesvd_(jobu, jobvt, &A->nrows, &A->ncols, A->data->data, &(int){1}, &(int){1}, A->desc,
    svd->d->data, svd->u->data->data, &(int){1}, &(int){1}, svd->u->desc,
    svd->vt->data->data, &(int){1}, &(int){1}, svd->vt->desc,
    &tmp, &lwork, &svd->info,
    (FC_LEN_T) strlen(jobu), (FC_LEN_T) strlen(jobvt));
#else
  pdgesvd_(jobu, jobvt, &A->nrows, &A->ncols, A->data->data, &(int){1}, &(int){1}, A->desc,
    svd->d->data, svd->u->data->data, &(int){1}, &(int){1}, svd->u->desc,
    svd->vt->data->data, &(int){1}, &(int){1}, svd->vt->desc,
    &tmp, &lwork, &svd->info);
#endif
  
  lwork = (int) tmp;
  work = malloc(lwork * sizeof(*work));
  if (work == NULL)
    return -1;
  
#ifdef FC_LEN_T
  pdgesvd_(jobu, jobvt, &A->nrows, &A->ncols, A->data->data, &(int){1}, &(int){1}, A->desc,
    svd->d->data, svd->u->data->data, &(int){1}, &(int){1}, svd->u->desc,
    svd->vt->data->data, &(int){1}, &(int){1}, svd->vt->desc,
    work, &lwork, &svd->info,
    (FC_LEN_T) strlen(jobu), (FC_LEN_T) strlen(jobvt));
#else
  pdgesvd_(jobu, jobvt, &A->nrows, &A->ncols, A->data->data, &(int){1}, &(int){1}, A->desc,
    svd->d->data, svd->u->data->data, &(int){1}, &(int){1}, svd->u->desc,
    svd->vt->data->data, &(int){1}, &(int){1}, svd->vt->desc,
    work, &lwork, &svd->info);
#endif
  
  free(work);
  
  return 0;
}



// svd with no QR
SEXP R_PDGESVD(SEXP M, SEXP N, SEXP ASIZE, SEXP A, SEXP DESCA, 
    SEXP ULDIM, SEXP DESCU, SEXP VTLDIM, SEXP DESCVT, SEXP JOBU, SEXP JOBVT, 
    SEXP INPLACE)
{
  UNUSED(INPLACE);
  
  R_INIT;
  SEXP RET, RET_NAMES, INFO, D, U, VT;
  double *A_cp;
  dsvd_t svd;
  dvector_t d;
  dmatrix_t a_local, u_local, vt_local;
  ddmatrix_t a, u, vt;
  
  newRvec(INFO, 1, "int");
  newRvec(D, INT(ASIZE, 0), "dbl");
  newRmat(U, INT(ULDIM, 0), INT(ULDIM, 1), "dbl");
  newRmat(VT, INT(VTLDIM, 0), INT(VTLDIM, 1), "dbl");
  
  a_local.nrows = nrows(A);
  a_local.ncols = ncols(A);
  
  a.data = &a_local;
  a.nrows = INT(M);
  a.ncols = INT(N);
  a.desc = INTP(DESCA);
  
  d.data = DBLP(D);
  
  u_local.data = DBLP(U);
  u.data = &u_local;
  u.desc = INTP(DESCU);
  
  vt_local.data = DBLP(VT);
  vt.data = &vt_local;
  vt.desc = INTP(DESCVT);
  
  svd.d = &d;
  svd.u = &u;
  svd.vt = &vt;
  
  
  A_cp = (double *) malloc(nrows(A)*ncols(A) * sizeof(double));
  if (A_cp == NULL)
    goto oom;
  memcpy(A_cp, REAL(A), nrows(A)*ncols(A)*sizeof(double));
  
  a_local.data = A_cp;
  
  int check = svd_noqr(STR(JOBU), STR(JOBVT), &a, &svd);
  
  free(A_cp);
  if (check)
    goto oom;
  
  // Manage return
  INT(INFO) = svd.info;
  make_list_names(RET_NAMES, 4, "info", "d", "u", "vt");
  make_list(RET, RET_NAMES, 4, INFO, D, U, VT);
  
  R_END;
  return RET;
  
  
  
oom:
  R_END;
  error("out of memory");
  
  return R_NilValue; // for static analysis
}
