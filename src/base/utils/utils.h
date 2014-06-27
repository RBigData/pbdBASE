// blacs_util.f90
void optimalgrid_(int *nprocs, int *nrows, int *ncols);
void dallreduce_(double *x, int *descx, char *op, char *scope);
void dreduce_(double *x, int *descx, char *op, int *rdest, int *cdest, char *scope);
void iallreduce_(int *x, int *descx, char *op, char *scope);
void ireduce_(int *x, int *descx, char *op, int *rdest, int *cdest, char *scope);


// dmat_redist.f
void mksubmat_(double *gblx, double *subx, int *descx);
void mkgblmat_(double *gbls, double *subx, int *descx, int *rdest, int *cdest);

// indices.f90
void numrocwrap_(int *n, int *nb, int *iproc, int *nprocs, int *num);
void pdims_(int *desc, int *ldm, int *blacs);
void l2gpair_(int *i, int *j, int *gi, int *gj, int *desc, int *blacs);
void g2lpair_(int *i, int *j, int *gi, int *gj, int *desc, int *blacs);

// putil.f
void ptri2zero_(char *uplo, char *diag, double *x, int *descx);
void pdmksym_(char *uplo, double *x, int *ix, int *jx, int *descx);
void pdgdgtk_(double *x, int *ix, int *jx, int *descx, double *diag, int *rdest, int *cdest);
void pddiagmk_(double *x, int *ix, int *jx, int *descx, double *diag, int *ldiag);

// scale.f90
void pdsweep_(double *x, int *ix, int *jx, int *descx, double *vec, int *lvec, int *margin, char *fun);

// special.f90
void dhilbmk_(int *n, double *x);
void pdhilbmk_(double *x, int *descx);
void pdmkcpn1_(double *x, int *descx, double *coef);

// util.f90
void dmksym_(char *triang, int *m, int *n, double *x);
