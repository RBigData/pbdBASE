//WCC: Header file for all C code called by .Call to wrap .Fortran.

#ifndef __BASE_GLOBAL__
#define __BASE_GLOBAL__

#include <R.h>
#include <Rinternals.h>

/* Obtain character pointers. */
#define CHARPT(x,i)	((char*)CHAR(STRING_ELT(x,i)))

#ifdef  __cplusplus
extern "C" {
#endif

/* Functions in "mpi_blacs.f". */
extern void F77_NAME(mpi_blacs_initialize)(int * nprow, int *npcol, int *ictxt,
				int *myrow, int *mycol);
extern void F77_NAME(row_col_sums)(int *ictxt, const char *scope, int *m,
				int *n, int *lda, double *a);



/* Functions in "pblas_level3.f". */
//extern void F77_NAME(rpdtran)(double *a, double *c, int *ictxt, int *myrow,
//				int *mycol, int *desca, int *descc,
//				int *m, int *n);
extern void F77_NAME(rpdgemm)(double *a, double *b, double *c,
				int *ictxt, int *myrow, int *mycol,
				int *desca, int *descb, int *descc,
				int *m, int *n, int *k);

extern void F77_NAME(pdtran)(int m, int n, double one, 
        double *a, int ia, int ja, int *desca, double zero,
        double *c, int ic, int jc, int *descc);



/* Functions in "scalapack.f". */
extern void F77_NAME(rpdgesv)(double *a, double *b, int *ictxt, int *myrow,
				int *mycol, int *desca, int *descb, int *n,
				int *nrhs, int *mxldims, int *info);
extern void F77_NAME(rpdgesvdsz)(int *m, int *n, int *asize, int *ictxt,
				int *myrow, int *mycol, int *desca, int *descu,
				int *descvt, double *temp, int *info,
				const char *jobu, const char *jobvt);
extern void F77_NAME(rpdgesvd)(int *m, int *n, int *asize, int *ictxt,
				int *myrow, int *mycol, double *a, int *desca,
				double *d, double *u, int *descu, double *vt,
				int *descvt, int *info, int *lwork,
				const char *jobu, const char *jobvt);
extern void F77_NAME(rpdgetrisz)(int *ictxt, int *myrow, int *mycol, int *desca,
				int *n, int *info, double *temp, int *itemp);
extern void F77_NAME(rpdgetri)(double *a, int *ictxt, int *myrow, int *mycol,
				int *desca, int *n, int *info, int *lwork,
				int *liwork);
extern void F77_NAME(rpdgetrf)(double *a, int *ictxt, int *myrow, int *mycol,
				int *desca, int *m, int *n, int *lipiv,
				int *info);
extern void F77_NAME(rpdpotrf)(double *a, int *ictxt, int *myrow, int *mycol,
				int *desca, int *n, const char *uplo,
				int *info);



/* Functions in "scalapack_utility_fort.f". */
extern void F77_NAME(rpdlaprnt)(int *m, int *n, double *a, int *ia, int *ja,
				int *desca, int *irprnt, int *icprnt,
				const char *cmatnm, int *nout, int *ictxt,
				int *myrow, int *mycol);
//extern void F77_NAME(rpdgemr2d)(int *m, int *n, double *x, int *descx,
//				double *b, int *descb, int *ctxt, int *locrx,
//				int *loccx, int *locrb, int *loccb);
extern void Cpdgemr2d(int m, int n, double *x, int ia, int ja, 
        int * descx, double *b, int ib, int jb, int *descb, int ctxt);


#ifdef  __cplusplus
}
#endif

#endif
