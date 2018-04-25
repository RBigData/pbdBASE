#ifndef __PBDBASE_PACKAGE__
#define __PBDBASE_PACKAGE__


#include <RNACI.h>
#include <mpi.h>

#include "base/linalg/linalg.h"
#include "base/stats/stats.h"
#include "base/utils/utils.h"
#include "scalapack.h"

#define CHARPT(x,i)	((char*)CHAR(STRING_ELT(x,i)))
#define nonzero(x) (x?x:1)
#define MIN(a,b) (a<b?a:b)
#define UNUSED(x) (void)(x)

extern void F77_NAME(mpi_blacs_initialize)(int * nprow, int *npcol, int *ictxt, int *myrow, int *mycol);

void pdsweep(double *restrict x, const int ix, const int jx, int *restrict descx, double *restrict vec, const int lvec, const int margin, const char fun);


#endif
