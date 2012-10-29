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

#ifdef  __cplusplus
}
#endif

#endif
