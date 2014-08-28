#ifndef __DMAT_MATEXP_H__
#define __DMAT_MATEXP_H__


// *LAPACK functions
void dgemm_(char *transa, char *transb, int *m, int *n, int *k, double *alpha, double *a, int *lda, double *b, int *ldb, double *beta, double *c, int *ldc);
void dlacpy_(char *uplo, int *m, int *n, double *a, int *lda, double *b, int *ldb);
void dscal_(int *n, double *a, double *x, int *incx);
void dgesv_(int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info);


void pdgemm_(char *transa, char *transb, int *m, int *n, int *k, double *alpha, double *a, int *ia, int *ja, int *desca, double *b, int *ib, int *jb, int *descb, double *beta, double *c, int *ic, int *jc, int *descc);
void pdlacpy_(char *uplo, int *m, int *n, double *a, int *ia, int *ja, int *desca, double *b, int *ib, int *jb, int *descb);

// pbdBASE functions
void pdims_(int *desc, int *ldm, int *blacs);
void l2gpair_(int *i, int *j, int *gi, int *gj, int *desc, int *blacs);


// Defines
#define SGNEXP(x,pow) (x==0?(pow==0?1:0):(x>0?1:(pow%2==0?1:(-1))))
#define MIN(a,b) (a<b?a:b)


// Functions
void matexp(int n, const int p, double *x, double *ret);

void p_matexp_pade(double *A, int *desca, int p, double *N, double *D);
void p_matpow_by_squaring(double *A, int *desca, int b, double *P);


// Constants
extern const double matexp_pade_coefs[14];


#endif
