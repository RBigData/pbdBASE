#include <R.h>
#include <Rinternals.h>
#include "base_global.h"

// PDTRAN
SEXP R_PDTRAN(SEXP M, SEXP N, SEXP A, SEXP DESCA, SEXP CLDIM, SEXP DESCC)
{
  SEXP C;
  PROTECT(C = allocMatrix(REALSXP, INTEGER(CLDIM)[0], INTEGER(CLDIM)[1]));
  
  const int IJ = 1;
  const double one = 1.0;
  const double zero = 0.0;
  
  F77_CALL(pdtran)(INTEGER(M), INTEGER(N), &one, 
      REAL(A), &IJ, &IJ, INTEGER(DESCA), &zero, 
      REAL(C), &IJ, &IJ, INTEGER(DESCC));
  
  UNPROTECT(1);
  return(C);
}

// PDGEMM
SEXP R_PDGEMM(SEXP TRANSA, SEXP TRANSB, SEXP M, SEXP N, SEXP K,
  SEXP A, SEXP DESCA, SEXP B, SEXP DESCB, SEXP CLDIM, SEXP DESCC)
{
  SEXP C;
  PROTECT(C = allocMatrix(REALSXP, INTEGER(CLDIM)[0], INTEGER(CLDIM)[1]));
  
  const double alpha = 1.0;
  const double beta = 0.0;
  const int IJ = 1;
  
  F77_CALL(pdgemm)(CHARPT(TRANSA, 0), CHARPT(TRANSB, 0),
    INTEGER(M), INTEGER(N), INTEGER(K), &alpha, 
    REAL(A), &IJ, &IJ, INTEGER(DESCA),
    REAL(B), &IJ, &IJ, INTEGER(DESCB), &beta,
    REAL(C), &IJ, &IJ, INTEGER(DESCC));
  
  UNPROTECT(1);
  return(C);
}

// PDSYRK
static void c_l2g_coord(int* ret, int i, int j, int* dim, int* bldim, int* procs, int myproc)
{
  const int nprocs = procs[0] * procs[1];
  ret[0] = nprocs*bldim[0] * (i-1)/bldim[0] + (i-1)%bldim[0] + ((nprocs+myproc)%nprocs)*bldim[0] + 1;
  ret[1] = nprocs*bldim[1] * (j-1)/bldim[1] + (j-1)%bldim[1] + ((nprocs+myproc)%nprocs)*bldim[1] + 1;
}


SEXP R_PDSYRK(SEXP UPLO, SEXP TRANS, SEXP N, SEXP K, 
  SEXP A, SEXP DESCA, SEXP CLDIM, SEXP DESCC, 
  SEXP dim, SEXP bldim, SEXP myproc, SEXP procs, SEXP src)
{
  SEXP C;
  PROTECT(C = allocMatrix(REALSXP, INTEGER(CLDIM)[0], INTEGER(CLDIM)[1]));
  
  int i, j;
  const double alpha = 1.0;
  const double beta = 0.0;
  const int IJ = 1;
  const int size = INTEGER(N)[0];
  const char trans = 'T';
  int ret[6];
  
  pdsyrk_(CHARPT(UPLO, 0), CHARPT(TRANS, 0),
    &size, INTEGER(K), &alpha, 
    REAL(A), &IJ, &IJ, INTEGER(DESCA), &beta,
    REAL(C), &IJ, &IJ, INTEGER(DESCC));
  
  // Copy over the upper triangle to the lower
  pdgeadd_(&trans, &size, &size, &alpha, REAL(C), &IJ, &IJ, INTEGER(DESCC), 
    &alpha, REAL(C), &IJ, &IJ, INTEGER(DESCC));
  
  for (j=0; j<size; j++){
    for (i=0; i<size; i++){
      c_l2g_coord(ret, i, j, INTEGER(dim), INTEGER(bldim), INTEGER(procs), INTEGER(myproc));
      if (INTEGER(myproc)[0]==ret[2] && INTEGER(myproc)[1]==ret[3]){
/*        Rprintf("%d  ", ret[4]);*/
/*        Rprintf("%d\n", ret[5]);*/
        REAL(C)[ret[4] + size*ret[5]] *= 0.5;
      }
    }
  }
  
  UNPROTECT(1);
  return(C);
}


