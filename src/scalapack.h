// BLACS



// PBLAS
void pdtran_(int *m, int *n, double *alpha, double *a, int *ia, int *ja, int *desca, double *beta, double *c, int *ic, int *jc, int *descc);
void pdgemm_(char *transa, char *transb, int *m, int *n, int *k, double *alpha, double *a, int *ia, int *ja, int *desca, double *b, int *ib, int *jb, int *descb, double *beta, double *c, int *ic, int *jc, int *descc);


// SCALAPACK
void pdgesv_(int *n, int *nrhs, double *a, int *ia, int *ja, double *desca, int *ipiv, double *b, int *ib, int *jb, int *descb, int *info);
void pdgesvd_(char *jobu, char *jobvt, int *m, int *n, double *a, int *ia, int *ja, int *desca, double *s, double *u, int *iu, int *ju, int *descu, double *vt, int *ivt, int *jvt, int *descvt, double *work, int *lwork, int *info);
void pdsyev_(char *jobz, char *uplo, int *n, double *a, int *ia, int *ja, int *desca, double *w, double *z, int *iz, int *jz, int *descz, double *work, int *lwork, int *info);
void pdgetrf_(int *m, int *n, double *a, int *ia, int *ja, int *desca, int *ipiv, int *info);
void pdpotrf_(char *uplo, int *n, double *a, int *ia, int *ja, int *desca, int *info);
void pdsyevx_(char *jobz, char *range, char *uplo, int *n, double *a, int *ia, int *ja, int *desca, double *vl, double *vu, int *il, int *iu, double *abstol, int *m, int *nz, double *w, double *orfac, double *z, int *iz, int *jz, int *descz, double *work, int *lwork, int *iwork, int *liwork, int *ifail, int *iclustr, double *gap, int *info);
void pdtrcon_(char *norm, char *uplo, char *diag, int *n, double *a, int *ia, int *ja, int *desca, double *rcond, double *work, int *lwork, int *iwork, int *liwork, int *info);
void pdormqr_(char *side, char *trans, int *m, int *n, int *k, double *a, int *ia, int *ja, int *desca, double *tau, double *c, int *ic, int *jc, int *descc, double *work, int *lwork, int *info);
void pdorgqr_(int *m, int *n, int *k, double *a, int *ia, int *ja, int *desca, double *tau, double *work, int *lwork, int *info);


// TOOLS

