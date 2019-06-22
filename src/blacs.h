#ifndef __PBDBASE_BLACS__
#define __PBDBASE_BLACS__


// BLACS
void sl_init_(int *ictxt, int *nprow, int *npcol);
void blacs_gridinfo_(int *ictxt, int *nprow, int *npcol, int *myrow,
  int *mycol);
void Cigsum2d(int ConTxt, char *scope, char *top, int m, int n, int *A,
  int lda, int rdest, int cdest);
void Cdgsum2d(int ConTxt, char *scope, char *top, int m, int n, double *A,
  int lda, int rdest, int cdest);
void Cigamx2d(int ConTxt, char *scope, char *top, int m, int n, int *A,
  int lda, int *rA, int *cA, int ldia, int rdest, int cdest);
void Cdgamx2d(int ConTxt, char *scope, char *top, int m, int n, double *A,
  int lda, int *rA, int *cA, int ldia, int rdest, int cdest);
void Cigamn2d(int ConTxt, char *scope, char *top, int m, int n, int *A,
  int lda, int *rA, int *cA, int ldia, int rdest, int cdest);
void Cdgamn2d(int ConTxt, char *scope, char *top, int m, int n, double *A,
  int lda, int *rA, int *cA, int ldia, int rdest, int cdest);
void Cdgesd2d(int ConTxt, int m, int n, double *A, int lda,
  int rdest, int cdest);
void Cdgerv2d(int ConTxt, int m, int n, double *A, int lda, int rsrc, int csrc);

// REDIST
void pigemr2d_(int *m, int *n, int *a, int *ia, int *ja, int *desca,
  int *b, int *ib, int *jb, int *descb, int *ictxt);
void Cpigemr2d(int m, int n, int *a, int ia, int ja, int *desca,
  int *b, int ib, int jb, int *descb, int ictxt);
void pdgemr2d_(int *m, int *n, double *a, int *ia, int *ja, int *desca,
  double *b, int *ib, int *jb, int *descb, int *ictxt);
void Cpdgemr2d(int m, int n, double *a, int ia, int ja, int *desca,
  double *b, int ib, int jb, int *descb, int ictxt);

// ???
int numroc_(const int* n, const int* nb, const int* iproc, const int* isrcproc, const int* nprocs);
void Cblacs_gridinfo(int ConTxt, int *nprow, int *npcol, int *myrow, int *mycol);


#endif
