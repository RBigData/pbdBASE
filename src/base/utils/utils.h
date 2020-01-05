#ifndef __PBDBASE_BASE_UTILS__
#define __PBDBASE_BASE_UTILS__


// For C/Fortran char* string lengths using size_t
#ifdef USE_FC_LEN_T
  #include <stddef.h>
  #include <Rconfig.h>    // this defines FC_LEN_T
  #include <string.h>
#endif


// blacs_util.f90
void optimalgrid_(int *nprocs, int *nrows, int *ncols);
#ifdef FC_LEN_T
  void dallreduce_(double *x, int *descx, char *op, char *scope, FC_LEN_T op_len, FC_LEN_T scope_len);
#else
  void dallreduce_(double *x, int *descx, char *op, char *scope);
#endif
#ifdef FC_LEN_T
  void ireduce_(int *x, int *descx, char *op, int *rdest, int *cdest, char *scope, FC_LEN_T op_len, FC_LEN_T scope_len);
#else
  void ireduce_(int *x, int *descx, char *op, int *rdest, int *cdest, char *scope);
#endif


// dmat_redist.f
void mksubmat_(double *gblx, double *subx, int *descx);
void mkgblmat_(double *gbls, double *subx, int *descx, int *rdest, int *cdest);


// indices.f90
void numrocwrap_(int *n, int *nb, int *iproc, int *nprocs, int *num);
void pdims_(const int *const restrict desc, int *const restrict ldm, int *const restrict blacs);
void l2gpair_(const int *const restrict i, const int *const restrict j, int *const restrict gi, int *const restrict gj, const int *const restrict desc, const int *const restrict blacs);
void g2lpair_(int *i, int *j, int *gi, int *gj, int *desc, int *blacs);


// putil.f
#ifdef FC_LEN_T
  void ptri2zero_(char *uplo, char *diag, double *x, int *descx, FC_LEN_T uplo_len, FC_LEN_T diag_len);
#else
  void ptri2zero_(char *uplo, char *diag, double *x, int *descx);
#endif
#ifdef FC_LEN_T
  void pdmksym_(char *uplo, double *x, int *ix, int *jx, int *descx, FC_LEN_T uplo_len);
#else
  void pdmksym_(char *uplo, double *x, int *ix, int *jx, int *descx);
#endif
void pdgdgtk_(double *x, int *ix, int *jx, int *descx, double *diag, int *rdest, int *cdest);
void pddiagmk_(double *x, int *ix, int *jx, int *descx, double *diag, int *ldiag);


// scale.f90
#ifdef FC_LEN_T
  void pdsweep_(double *x, int *ix, int *jx, int *descx, double *vec, int *lvec, int *margin, char *fun, FC_LEN_T fun_len);
#else
  void pdsweep_(double *x, int *ix, int *jx, int *descx, double *vec, int *lvec, int *margin, char *fun);
#endif


// special.f90
void dhilbmk_(int *n, double *x);
void pdhilbmk_(double *x, int *descx);
void pdmkcpn1_(double *x, int *descx, double *coef);


// util.f90
#ifdef FC_LEN_T
  void dmksym_(char *triang, int *m, int *n, double *x, FC_LEN_T triang_len);
#else
  void dmksym_(char *triang, int *m, int *n, double *x);
#endif


#endif
