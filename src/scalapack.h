#ifndef __PBDBASE_SCALAPACK__
#define __PBDBASE_SCALAPACK__


// BLACS
void sl_init_(int *ictxt, int *nprow, int *npcol);
void blacs_gridinfo_(int *ictxt, int *nprow, int *npcol, int *myrow,
  int *mycol);
void Cigsum2d(int ConTxt, char *scope, char *top, int m, int n, int *A,
  int lda, int rdest, int cdest);
void Cdgsum2d(int ConTxt, char *scope, char *top, int m, int n, double *A,
  int lda, int rdest, int cdest);
void Cigamx2d(int ConTxt, char *scope, char *top, int m, int n, int *A,
  int lda, int *rA, int *cA, int ldia, int rdest, int cdest);
void Cdgamx2d(int ConTxt, char *scope, char *top, int m, int n, double *A,
  int lda, int *rA, int *cA, int ldia, int rdest, int cdest);
void Cigamn2d(int ConTxt, char *scope, char *top, int m, int n, int *A,
  int lda, int *rA, int *cA, int ldia, int rdest, int cdest);
void Cdgamn2d(int ConTxt, char *scope, char *top, int m, int n, double *A,
  int lda, int *rA, int *cA, int ldia, int rdest, int cdest);
void Cdgesd2d(int ConTxt, int m, int n, double *A, int lda,
  int rdest, int cdest);
void Cdgerv2d(int ConTxt, int m, int n, double *A, int lda, int rsrc, int csrc);


// PBLAS
void pdtran_(int *m, int *n, double *alpha, double *a, int *ia, int *ja, 
  int *desca, double *beta, double *c, int *ic, int *jc, int *descc);
void pdgemm_(char *transa, char *transb, int *m, int *n, int *k, double *alpha, 
  double *a, int *ia, int *ja, int *desca, double *b, int *ib, int *jb, 
  int *descb, double *beta, double *c, int *ic, int *jc, int *descc);


// SCALAPACK
void pdgesv_(int *n, int *nrhs, double *a, int *ia, int *ja, int *desca, 
  int *ipiv, double *b, int *ib, int *jb, int *descb, int *info);
void pdgesvd_(const char *restrict jobu, const char *restrict jobvt,
  const int *restrict m, const int *restrict n, double *restrict a,
  const int *restrict ia, const int *restrict ja, const int *restrict desca,
  double *restrict s, double *restrict u, const int *restrict iu,
  const int *restrict ju, const int *restrict descu, double *restrict vt,
  const int *restrict ivt, const int *restrict jvt, const int *restrict descvt,
  double *restrict work, const int *restrict lwork, int *restrict info);
void pdsyev_(char *jobz, char *uplo, int *n, double *a, int *ia, int *ja, 
  int *desca, double *w, double *z, int *iz, int *jz, int *descz, double *work, 
  int *lwork, int *info);
void pdsyevr_(char *jobz, char *range, char *uplo, int  *n, double *a, int *ia,
  int *ja, int *desca, double *vl, double *vu, int *il, int *iu, int *m, int *nz,
  double *w, double *z, int *iz, int *jz, int *descz, double *work, int *lwork,
  int *iwork, int  *liwork, int *info);
void pdgetrf_(int *m, int *n, double *a, int *ia, int *ja, int *desca, 
  int *ipiv, int *info);
void pdpotrf_(char *uplo, int *n, double *a, int *ia, int *ja, int *desca, 
  int *info);
void pdsyevx_(char *jobz, char *range, char *uplo, int *n, double *a, int *ia, 
  int *ja, int *desca, double *vl, double *vu, int *il, int *iu, 
  double *abstol, int *m, int *nz, double *w, double *orfac, double *z, int *iz,
  int *jz, int *descz, double *work, int *lwork, int *iwork, int *liwork, 
  int *ifail, int *iclustr, double *gap, int *info);
void pdtrcon_(char *norm, char *uplo, char *diag, int *n, double *a, int *ia, 
  int *ja, int *desca, double *rcond, double *work, int *lwork, int *iwork, 
  int *liwork, int *info);
void pdormqr_(char *side, char *trans, int *m, int *n, int *k, double *a, 
  int *ia, int *ja, int *desca, double *tau, double *c, int *ic, int *jc, 
  int *descc, double *work, int *lwork, int *info);
void pdorgqr_(int *m, int *n, int *k, double *a, int *ia, int *ja, int *desca, 
  double *tau, double *work, int *lwork, int *info);
void pdgelqf_(int *m, int *n, double *a, int *ia, int *ja, int *desca,
  double *tau, double *work, int *lwork, int *info);
void pdorglq_(int *m, int *n, int *k, double *a, int *ia, int *ja, int *desca,
  double *tau, double *work, int *lwork, int *info);


// TOOLS
void bprnt_(int *m, int *n, double *a, int *ia, int *ja, int *desca,
  int *irprnt, int *icprnt, char *chatnm, int *nout, double *work);
void descinit_(int *desc, int *m, int *n, int *mb, int *nb, int *isrc,
  int *icsrc, int *ictxt, int *lld, int *info);


// REDIST
void pigemr2d_(int *m, int *n, int *a, int *ia, int *ja, int *desca,
  int *b, int *ib, int *jb, int *descb, int *ictxt);
void Cpigemr2d(int m, int n, int *a, int ia, int ja, int *desca,
  int *b, int ib, int jb, int *descb, int ictxt);
void pdgemr2d_(int *m, int *n, double *a, int *ia, int *ja, int *desca,
  double *b, int *ib, int *jb, int *descb, int *ictxt);
void Cpdgemr2d(int m, int n, double *a, int ia, int ja, int *desca,
  double *b, int ib, int jb, int *descb, int ictxt);


// ???
int numroc_(const int* n, const int* nb, const int* iproc, const int* isrcproc, const int* nprocs);
void Cblacs_gridinfo(int ConTxt, int *nprow, int *npcol, int *myrow, int *mycol);


#endif
