#ifndef __BASE_GLOBAL__
#define __BASE_GLOBAL__


#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

#include "scalapack.h"

#include "base/linalg/linalg.h"
#include "base/stats/stats.h"
#include "base/utils/utils.h"


// Obtain character pointers
#define CHARPT(x,i)	((char*)CHAR(STRING_ELT(x,i)))

#define nonzero(x) (x?x:1)

#define MIN(a,b) (a<b?a:b)

/* Functions in "mpi_blacs.f". */
extern void F77_NAME(mpi_blacs_initialize)(int * nprow, int *npcol, int *ictxt,
  int *myrow, int *mycol);


#endif
