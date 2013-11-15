#ifndef __BASE_GLOBAL__
#define __BASE_GLOBAL__


#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>


// Obtain int/double pointers
#define INT(x,i) (INTEGER(x)[i])
#define RL(x,i) (REAL(x)[i])

// Obtain character pointers
#define CHARPT(x,i)	((char*)CHAR(STRING_ELT(x,i)))

#ifdef  __cplusplus
extern "C" {
#endif

#define nonzero(x) (x?x:1)

/* Functions aux.f. */
extern void F77_NAME(matnorm)(double* val, char* norm, int* m, int* n,
  double* a, int* ia, int* ja, int* desca);

extern void F77_NAME(condnum)(char* norm, int* m, int* n, double* a, 
  int* ia, int* ja, int* desca, double* ret, int* info);

/* Functions in "mpi_blacs.f". */
extern void F77_NAME(mpi_blacs_initialize)(int * nprow, int *npcol, int *ictxt,
  int *myrow, int *mycol);

#ifdef  __cplusplus
}
#endif


#endif
