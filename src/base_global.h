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


/* Functions in "scalapack.f". */
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
