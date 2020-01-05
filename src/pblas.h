#ifndef __PBDBASE_PBLAS__
#define __PBDBASE_PBLAS__


// For C/Fortran char* string lengths using size_t
#ifdef USE_FC_LEN_T
  #include <stddef.h>
  #include <Rconfig.h>    // this defines FC_LEN_T
  #include <string.h>
#endif


// PBLAS
void pdtran_(int *m, int *n, double *alpha, double *a, int *ia, int *ja, 
  int *desca, double *beta, double *c, int *ic, int *jc, int *descc);
#ifdef FC_LEN_T
  void pdgemm_(char *transa, char *transb, int *m, int *n, int *k, double *alpha, 
    double *a, int *ia, int *ja, int *desca, double *b, int *ib, int *jb, 
    int *descb, double *beta, double *c, int *ic, int *jc, int *descc,
    FC_LEN_T transa_len, FC_LEN_T transb_len);
#else
  void pdgemm_(char *transa, char *transb, int *m, int *n, int *k, double *alpha, 
    double *a, int *ia, int *ja, int *desca, double *b, int *ib, int *jb, 
    int *descb, double *beta, double *c, int *ic, int *jc, int *descc);
#endif


#endif
