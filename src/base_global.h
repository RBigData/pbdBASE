#ifndef __BASE_GLOBAL__
#define __BASE_GLOBAL__

#include <R.h>
#include <Rinternals.h>

/* Obtain character pointers. */
#define CHARPT(x,i)	((char*)CHAR(STRING_ELT(x,i)))

#ifdef  __cplusplus
extern "C" {
#endif

/* return 1 if x is 0, x otherwise */
/* need this when alloc'ing work vectors for fortran */
#define nonzero(x) (x?x:1)

/* Functions in "mpi_blacs.f". */
extern void F77_NAME(mpi_blacs_initialize)(int * nprow, int *npcol, int *ictxt,
				int *myrow, int *mycol);

#ifdef  __cplusplus
}
#endif

#endif
