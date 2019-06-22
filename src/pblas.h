#ifndef __PBDBASE_PBLAS__
#define __PBDBASE_PBLAS__


// PBLAS
void pdtran_(int *m, int *n, double *alpha, double *a, int *ia, int *ja, 
  int *desca, double *beta, double *c, int *ic, int *jc, int *descc);
void pdgemm_(char *transa, char *transb, int *m, int *n, int *k, double *alpha, 
  double *a, int *ia, int *ja, int *desca, double *b, int *ib, int *jb, 
  int *descb, double *beta, double *c, int *ic, int *jc, int *descc);


#endif
