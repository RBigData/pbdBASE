//WCC: wrapers for "mpi_blacs.f".

#include <R.h>
#include <Rinternals.h>
#include "base_global.h"

SEXP R_mpi_blacs_initialize(SEXP NPROW, SEXP NPCOL, SEXP ICTXT, SEXP MYROW,
		SEXP MYCOL){
	F77_CALL(mpi_blacs_initialize)(INTEGER(NPROW), INTEGER(NPCOL),
			INTEGER(ICTXT), INTEGER(MYROW), INTEGER(MYCOL));
	return(R_NilValue);
} /* End of R_mpi_blacs_initialize(). */

SEXP R_row_col_sums(SEXP ICTXT, SEXP SCOPE, SEXP M, SEXP N, SEXP A, SEXP LDA){
	int i;
	double *pt_A = REAL(A), *pt_OUT;
	SEXP OUT;

	PROTECT(OUT = allocMatrix(REALSXP, INTEGER(M)[0], INTEGER(N)[0]));
	pt_OUT = REAL(OUT);
	for(i = 0; i < INTEGER(M)[0] * INTEGER(N)[0]; i++){
		*pt_OUT = *pt_A;
		pt_OUT++;
		pt_A++;
	}
	F77_CALL(row_col_sums)(INTEGER(ICTXT), CHARPT(SCOPE, 0),
			INTEGER(M), INTEGER(N), INTEGER(LDA), REAL(OUT));
	UNPROTECT(1);
	return(OUT);
} /* End of R_pdgemm(). */

