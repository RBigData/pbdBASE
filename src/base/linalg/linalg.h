// auxil.f90
void matnorm_(double *value, char *norm, int *m, int *n, double *a,
  int *ia, int *ja, int *desca);
void condnum_(char *norm, int *m, int *n, double *a, int *ia, int *ja,
  int *desca, double *rcond, int *info);

// pcrossprod.f90
void pdcrossprod_(char *uplo, char *trans, double *alpha, double *x,
  int *ix, int *jx, int *descx, double *c, int *ic, int *jc, int *descc);
void pdchtri_(char *uplo, double *x, int *ix, int *jx, int *descx,
  double *c, int *ic, int *jc, int *descc, int *info);
void pdinvip_(double *x, int *ix, int *jx, int *descx, int *info);
void pdinv_(double *x, int *ix, int *jx, int *descx, double *inv, int *info);

// plm.f
void rpdormqr_(char *side, char *trans, int *m, int *n, int *k, 
  double *a, int *ia, int *ja, int *desca, double *tau, double *c,
  int *ic, int *jc, int *descc, double *work, int *lwork, int *info);
void rpdgeqpf_(double *tol, int *m, int *n, double *a, int *ia, 
  int *ja, int *desca, int *ipiv, double *tau, double *work, 
  int *lwork, int *rank, int *info);
void rpdgels_(double *tol, char *trans, int *m, int *n, int *nrhs,
  double *a, int *ia, int *ja, int *desca, double *b, int *ib, 
  int *jb, int *descb, double *eff, double *ft, double *rsd, 
  double *tau, double *work, int *lwork, int *ipiv, int *rank,
  int *info);

// prblas.f90
void rl2blas_(double *x, int *ix, int *jx, int *descx, double *vec,
  int *lvec, int *fun);
void rl2insert_(double *x, int *ix, int *jx, int *descx, double *vec,
  int *lvec, int *indi, int *lindi, int *indj, int *lindj);
void rcolcpy_(double *x, int *descx, int *xcols, double *y, int *descy,
  int *ycols, int *lcols);
void rcolcpy2_(double *x, int *descx, int *xcols, int *lxcols, double *y,
  int *descy, int *ycols, int *lycols);
void rrowcpy_(double *x, int *descx, int *xrows, double *y, int *descy,
  int *yrows, int *lrows);
void rrowcpy2_(double *x, int *descx, int *xrows, int *lxrows, double *y,
  int *descy, int *yrows, int *lyrows);
void pdmvsum_(double *x, int *descx, double *y, int *descy);

