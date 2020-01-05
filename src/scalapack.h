#ifndef __PBDBASE_SCALAPACK__
#define __PBDBASE_SCALAPACK__


// For C/Fortran char* string lengths using size_t
#ifdef USE_FC_LEN_T
  #include <stddef.h>
  #include <Rconfig.h>    // this defines FC_LEN_T
  #include <string.h>
#endif


// PBLAS
void pdtran_(int *m, int *n, double *alpha, double *a, int *ia, int *ja, 
  int *desca, double *beta, double *c, int *ic, int *jc, int *descc);
#ifdef FC_LEN_T
  void pdgemm_(char *transa, char *transb, int *m, int *n, int *k, double *alpha, 
    double *a, int *ia, int *ja, int *desca, double *b, int *ib, int *jb, 
    int *descb, double *beta, double *c, int *ic, int *jc, int *descc,
    FC_LEN_T transa_len, FC_LEN_T transb_len);
#else
  void pdgemm_(char *transa, char *transb, int *m, int *n, int *k, double *alpha, 
    double *a, int *ia, int *ja, int *desca, double *b, int *ib, int *jb, 
    int *descb, double *beta, double *c, int *ic, int *jc, int *descc);
#endif


// SCALAPACK
void pdgesv_(int *n, int *nrhs, double *a, int *ia, int *ja, int *desca, 
  int *ipiv, double *b, int *ib, int *jb, int *descb, int *info);
#ifdef FC_LEN_T
  void pdgesvd_(const char *restrict jobu, const char *restrict jobvt,
    const int *restrict m, const int *restrict n, double *restrict a,
    const int *restrict ia, const int *restrict ja, const int *restrict desca,
    double *restrict s, double *restrict u, const int *restrict iu,
    const int *restrict ju, const int *restrict descu, double *restrict vt,
    const int *restrict ivt, const int *restrict jvt, const int *restrict descvt,
    double *restrict work, const int *restrict lwork, int *restrict info,
    FC_LEN_T jobu_len, FC_LEN_T jobvt_len);
#else
  void pdgesvd_(const char *restrict jobu, const char *restrict jobvt,
    const int *restrict m, const int *restrict n, double *restrict a,
    const int *restrict ia, const int *restrict ja, const int *restrict desca,
    double *restrict s, double *restrict u, const int *restrict iu,
    const int *restrict ju, const int *restrict descu, double *restrict vt,
    const int *restrict ivt, const int *restrict jvt, const int *restrict descvt,
    double *restrict work, const int *restrict lwork, int *restrict info);
#endif
#ifdef FC_LEN_T
  void pdsyev_(char *jobz, char *uplo, int *n, double *a, int *ia, int *ja, 
    int *desca, double *w, double *z, int *iz, int *jz, int *descz, double *work, 
    int *lwork, int *info,
    FC_LEN_T jobz_len, FC_LEN_T uplo_len);
#else
  void pdsyev_(char *jobz, char *uplo, int *n, double *a, int *ia, int *ja, 
    int *desca, double *w, double *z, int *iz, int *jz, int *descz, double *work, 
    int *lwork, int *info);
#endif
#ifdef FC_LEN_T
  void pdsyevr_(char *jobz, char *range, char *uplo, int  *n, double *a, int *ia,
    int *ja, int *desca, double *vl, double *vu, int *il, int *iu, int *m, int *nz,
    double *w, double *z, int *iz, int *jz, int *descz, double *work, int *lwork,
    int *iwork, int  *liwork, int *info,
    FC_LEN_T jobz_len, FC_LEN_T range_len, FC_LEN_T uplo_len);
#else
  void pdsyevr_(char *jobz, char *range, char *uplo, int  *n, double *a, int *ia,
    int *ja, int *desca, double *vl, double *vu, int *il, int *iu, int *m, int *nz,
    double *w, double *z, int *iz, int *jz, int *descz, double *work, int *lwork,
    int *iwork, int  *liwork, int *info);
#endif
void pdgetrf_(const int *const restrict m, const int *const restrict n,
  double *const restrict a, const int *const restrict ia,
  const int *const restrict ja, const int *const restrict desca, 
  int *const restrict ipiv, int *const restrict info);
#ifdef FC_LEN_T
  void pdpotrf_(char *uplo, int *n, double *a, int *ia, int *ja, int *desca, 
    int *info,
    FC_LEN_T uplo_len);
#else
  void pdpotrf_(char *uplo, int *n, double *a, int *ia, int *ja, int *desca, 
    int *info);
#endif
#ifdef FC_LEN_T
  void pdsyevx_(char *jobz, char *range, char *uplo, int *n, double *a, int *ia, 
    int *ja, int *desca, double *vl, double *vu, int *il, int *iu, 
    double *abstol, int *m, int *nz, double *w, double *orfac, double *z, int *iz,
    int *jz, int *descz, double *work, int *lwork, int *iwork, int *liwork, 
    int *ifail, int *iclustr, double *gap, int *info,
    FC_LEN_T jobz_len, FC_LEN_T range_len, FC_LEN_T uplo_len);
#else
  void pdsyevx_(char *jobz, char *range, char *uplo, int *n, double *a, int *ia, 
    int *ja, int *desca, double *vl, double *vu, int *il, int *iu, 
    double *abstol, int *m, int *nz, double *w, double *orfac, double *z, int *iz,
    int *jz, int *descz, double *work, int *lwork, int *iwork, int *liwork, 
    int *ifail, int *iclustr, double *gap, int *info);
#endif
#ifdef FC_LEN_T
  void pdtrcon_(char *norm, char *uplo, char *diag, int *n, double *a, int *ia, 
    int *ja, int *desca, double *rcond, double *work, int *lwork, int *iwork, 
    int *liwork, int *info,
    FC_LEN_T norm_len, FC_LEN_T uplo_len, FC_LEN_T diag_len);
#else
  void pdtrcon_(char *norm, char *uplo, char *diag, int *n, double *a, int *ia, 
    int *ja, int *desca, double *rcond, double *work, int *lwork, int *iwork, 
    int *liwork, int *info);
#endif
#ifdef FC_LEN_T
  void pdormqr_(char *side, char *trans, int *m, int *n, int *k, double *a, 
    int *ia, int *ja, int *desca, double *tau, double *c, int *ic, int *jc, 
    int *descc, double *work, int *lwork, int *info,
    FC_LEN_T side_len, FC_LEN_T trans_len);
#else
  void pdormqr_(char *side, char *trans, int *m, int *n, int *k, double *a, 
    int *ia, int *ja, int *desca, double *tau, double *c, int *ic, int *jc, 
    int *descc, double *work, int *lwork, int *info);
#endif
void pdorgqr_(int *m, int *n, int *k, double *a, int *ia, int *ja, int *desca, 
  double *tau, double *work, int *lwork, int *info);
void pdgelqf_(int *m, int *n, double *a, int *ia, int *ja, int *desca,
  double *tau, double *work, int *lwork, int *info);
void pdorglq_(int *m, int *n, int *k, double *a, int *ia, int *ja, int *desca,
  double *tau, double *work, int *lwork, int *info);


// TOOLS
#ifdef FC_LEN_T
  void bprnt_(int *m, int *n, double *a, int *ia, int *ja, int *desca,
    int *irprnt, int *icprnt, char *cmatnm, int *nout, double *work,
    FC_LEN_T cmatnm_len);
#else
  void bprnt_(int *m, int *n, double *a, int *ia, int *ja, int *desca,
    int *irprnt, int *icprnt, char *cmatnm, int *nout, double *work);
#endif
void descinit_(int *desc, int *m, int *n, int *mb, int *nb, int *isrc,
  int *icsrc, int *ictxt, int *lld, int *info);


#endif
