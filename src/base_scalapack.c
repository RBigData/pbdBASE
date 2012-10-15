//WCC: wrapers for "scalapack.f".

#include <R.h>
#include <Rinternals.h>
#include "base_global.h"

SEXP R_PDGESV(SEXP A, SEXP ALDIM, SEXP B, SEXP BLDIM,
		SEXP ICTXT, SEXP MYROW, SEXP MYCOL,
		SEXP DESCA, SEXP DESCB, SEXP N, SEXP NRHS,
		SEXP MXLDIMS){
	int i, *pt_ALDIM = INTEGER(ALDIM), *pt_BLDIM = INTEGER(BLDIM);
	double *pt_ORG, *pt_COPY, *A_OUT;
	SEXP RET, RET_NAMES, INFO, B_OUT;

	/* Protect R objects. */
	PROTECT(RET = allocVector(VECSXP, 2));
	PROTECT(RET_NAMES = allocVector(STRSXP, 2));
	PROTECT(INFO = allocVector(INTSXP, 1));
	PROTECT(B_OUT = allocMatrix(REALSXP, pt_BLDIM[0], pt_BLDIM[1]));

	SET_VECTOR_ELT(RET, 0, INFO);
	SET_VECTOR_ELT(RET, 1, B_OUT);
	SET_STRING_ELT(RET_NAMES, 0, mkChar("info")); 
	SET_STRING_ELT(RET_NAMES, 1, mkChar("B")); 
	setAttrib(RET, R_NamesSymbol, RET_NAMES);

	/* Set INFO, allocate A_OUT, copy A->A_OUT, and copy B -> B_OUT. */
	INTEGER(INFO)[0] = 0;
	A_OUT = (double *) R_alloc(pt_ALDIM[0] * pt_ALDIM[1], sizeof(double));
	pt_ORG = REAL(A);
	pt_COPY = A_OUT;
	for(i = 0; i < pt_ALDIM[0] * pt_ALDIM[1]; i++){
		*pt_COPY = *pt_ORG;
		pt_ORG++;
		pt_COPY++;
	}
	pt_ORG = REAL(B);
	pt_COPY = REAL(B_OUT);
	for(i = 0; i < pt_BLDIM[0] * pt_BLDIM[1]; i++){
		*pt_COPY = *pt_ORG;
		pt_ORG++;
		pt_COPY++;
	}

	/* Call Fortran. */
	F77_CALL(rpdgesv)(A_OUT, REAL(B_OUT),
			INTEGER(ICTXT), INTEGER(MYROW), INTEGER(MYCOL),
			INTEGER(DESCA), INTEGER(DESCB),
			INTEGER(N), INTEGER(NRHS), INTEGER(MXLDIMS),
			INTEGER(INFO));

	/* Return. */
	UNPROTECT(4);
	// Free(A_OUT);		// .Call takes care of this.
	return(RET);
} /* End of R_PDGESV(). */


SEXP R_PDGESVDSZ(SEXP M, SEXP N, SEXP ASIZE, SEXP ICTXT, SEXP MYROW,
                SEXP MYCOL, SEXP DESCA, SEXP DESCU, SEXP DESCVT,
                SEXP JOBU, SEXP JOBVT){
	SEXP RET, RET_NAMES, INFO, TEMP;

	/* Protect R objects. */
	PROTECT(RET = allocVector(VECSXP, 2));
	PROTECT(RET_NAMES = allocVector(STRSXP, 2));
	PROTECT(INFO = allocVector(INTSXP, 1));
	PROTECT(TEMP = allocVector(REALSXP, 1));

	SET_VECTOR_ELT(RET, 0, INFO);
	SET_VECTOR_ELT(RET, 1, TEMP);
	SET_STRING_ELT(RET_NAMES, 0, mkChar("info")); 
	SET_STRING_ELT(RET_NAMES, 1, mkChar("TEMP")); 
	setAttrib(RET, R_NamesSymbol, RET_NAMES);

	/* Set INFO. */
	INTEGER(INFO)[0] = 0;
	REAL(TEMP)[0] = 0.0;

	/* Call Fortran. */
        F77_CALL(rpdgesvdsz)(INTEGER(M), INTEGER(N), INTEGER(ASIZE),
                        INTEGER(ICTXT), INTEGER(MYROW), INTEGER(MYCOL),
                        INTEGER(DESCA), INTEGER(DESCU), INTEGER(DESCVT),
                        REAL(TEMP), INTEGER(INFO), CHARPT(JOBU, 0),
                        CHARPT(JOBVT, 0));

	/* Return. */
	UNPROTECT(4);
        return(RET);
} /* End of R_PDGESVDSZ(). */


SEXP R_PDGESVD(SEXP M, SEXP N, SEXP ASIZE, SEXP ICTXT, SEXP MYROW,
                SEXP MYCOL, SEXP A, SEXP DESCA, SEXP ALDIM, SEXP ULDIM,
                SEXP DESCU, SEXP VTLDIM, SEXP DESCVT, SEXP LWORK,
		SEXP JOBU, SEXP JOBVT){
	int i, *pt_ALDIM = INTEGER(ALDIM);
	double *A_OUT, *pt_ORG, *pt_COPY;
	SEXP RET, RET_NAMES, INFO, D, U, VT;

	/* Protect R objects. */
	PROTECT(RET = allocVector(VECSXP, 4));
	PROTECT(RET_NAMES = allocVector(STRSXP, 4));
	PROTECT(INFO = allocVector(INTSXP, 1));
	PROTECT(D = allocVector(REALSXP, INTEGER(ASIZE)[0]));
	PROTECT(U = allocMatrix(REALSXP, INTEGER(ULDIM)[0], INTEGER(ULDIM)[1]));
	PROTECT(VT = allocMatrix(REALSXP,
			INTEGER(VTLDIM)[0], INTEGER(VTLDIM)[1]));

	SET_VECTOR_ELT(RET, 0, INFO);
	SET_VECTOR_ELT(RET, 1, D);
	SET_VECTOR_ELT(RET, 2, U);
	SET_VECTOR_ELT(RET, 3, VT);
	SET_STRING_ELT(RET_NAMES, 0, mkChar("info")); 
	SET_STRING_ELT(RET_NAMES, 1, mkChar("d")); 
	SET_STRING_ELT(RET_NAMES, 2, mkChar("u")); 
	SET_STRING_ELT(RET_NAMES, 3, mkChar("vt")); 
	setAttrib(RET, R_NamesSymbol, RET_NAMES);

	/* Set INFO. */
	INTEGER(INFO)[0] = 0;
	A_OUT = (double *) R_alloc(pt_ALDIM[0] * pt_ALDIM[1], sizeof(double));
	pt_ORG = REAL(A);
	pt_COPY = A_OUT;
	for(i = 0; i < pt_ALDIM[0] * pt_ALDIM[1]; i++){
		*pt_COPY = *pt_ORG;
		pt_ORG++;
		pt_COPY++;
	}

	/* Call Fortran. */
        F77_CALL(rpdgesvd)(INTEGER(M), INTEGER(N), INTEGER(ASIZE),
                        INTEGER(ICTXT), INTEGER(MYROW), INTEGER(MYCOL),
                        A_OUT, INTEGER(DESCA), REAL(D), REAL(U),
			INTEGER(DESCU), REAL(VT), INTEGER(DESCVT),
			INTEGER(INFO), INTEGER(LWORK), CHARPT(JOBU, 0),
			CHARPT(JOBVT, 0));

	/* Return. */
	UNPROTECT(6);
        return(RET);
} /* End of R_PDGESVDSZ(). */


SEXP R_PDGETRISZ(SEXP ICTXT, SEXP MYROW, SEXP MYCOL, SEXP DESCA, SEXP N){
	SEXP RET, RET_NAMES, INFO, TEMP, ITEMP;

	/* Protect R objects. */
	PROTECT(RET = allocVector(VECSXP, 3));
	PROTECT(RET_NAMES = allocVector(STRSXP, 3));
	PROTECT(INFO = allocVector(INTSXP, 1));
	PROTECT(TEMP = allocVector(REALSXP, 1));
	PROTECT(ITEMP = allocVector(INTSXP, 1));

	SET_VECTOR_ELT(RET, 0, INFO);
	SET_VECTOR_ELT(RET, 1, TEMP);
	SET_VECTOR_ELT(RET, 2, ITEMP);
	SET_STRING_ELT(RET_NAMES, 0, mkChar("info")); 
	SET_STRING_ELT(RET_NAMES, 1, mkChar("TEMP")); 
	SET_STRING_ELT(RET_NAMES, 2, mkChar("ITEMP")); 
	setAttrib(RET, R_NamesSymbol, RET_NAMES);

	/* Set INFO and return R objects. */
	INTEGER(INFO)[0] = 0;
	REAL(TEMP)[0] = 0.0;
	INTEGER(ITEMP)[0] = 0;

	/* Call Fortran. */
        F77_CALL(rpdgetrisz)(INTEGER(ICTXT), INTEGER(MYROW), INTEGER(MYCOL),
                        INTEGER(DESCA), INTEGER(N), INTEGER(INFO),
                        REAL(TEMP), INTEGER(ITEMP));

	/* Return. */
	UNPROTECT(5);
        return(RET);
} /* End of R_PDGETRISZ(). */


SEXP R_PDGETRI(SEXP A, SEXP CLDIM, SEXP ICTXT, SEXP MYROW, SEXP MYCOL,
		SEXP DESCA, SEXP N, SEXP LWORK, SEXP LIWORK){
	int i, *pt_CLDIM = INTEGER(CLDIM);
	double *pt_A, *pt_C;
	SEXP RET, RET_NAMES, INFO, C;

	/* Protect R objects. */
	PROTECT(RET = allocVector(VECSXP, 2));
	PROTECT(RET_NAMES = allocVector(STRSXP, 2));
	PROTECT(INFO = allocVector(INTSXP, 1));
	PROTECT(C = allocMatrix(REALSXP, pt_CLDIM[0], pt_CLDIM[1]));

	SET_VECTOR_ELT(RET, 0, INFO);
	SET_VECTOR_ELT(RET, 1, C);
	SET_STRING_ELT(RET_NAMES, 0, mkChar("info")); 
	SET_STRING_ELT(RET_NAMES, 1, mkChar("A")); 
	setAttrib(RET, R_NamesSymbol, RET_NAMES);

	/* Copy A -> C and set INFO and return R objects. */
	INTEGER(INFO)[0] = 0;
	pt_A = REAL(A);
	pt_C = REAL(C);
	for(i = 0; i < pt_CLDIM[0] * pt_CLDIM[1]; i++){
		*pt_C = *pt_A;
		pt_A++;
		pt_C++;
	}

	/* Call Fortran. */
        F77_CALL(rpdgetri)(REAL(C), INTEGER(ICTXT), INTEGER(MYROW),
			INTEGER(MYCOL), INTEGER(DESCA), INTEGER(N),
			INTEGER(INFO), INTEGER(LWORK), INTEGER(LIWORK));

	/* Return. */
	UNPROTECT(4);
        return(RET);
} /* End of R_PDGETRI(). */


SEXP R_PDGETRF(SEXP A, SEXP CLDIM, SEXP ICTXT, SEXP MYROW, SEXP MYCOL,
		SEXP DESCA, SEXP M, SEXP N, SEXP LIPIV){
	int i, *pt_CLDIM = INTEGER(CLDIM);
	double *pt_A, *pt_C;
	SEXP RET, RET_NAMES, INFO, C;

	/* Protect R objects. */
	PROTECT(RET = allocVector(VECSXP, 2));
	PROTECT(RET_NAMES = allocVector(STRSXP, 2));
	PROTECT(INFO = allocVector(INTSXP, 1));
	PROTECT(C = allocMatrix(REALSXP, pt_CLDIM[0], pt_CLDIM[1]));

	SET_VECTOR_ELT(RET, 0, INFO);
	SET_VECTOR_ELT(RET, 1, C);
	SET_STRING_ELT(RET_NAMES, 0, mkChar("info")); 
	SET_STRING_ELT(RET_NAMES, 1, mkChar("A")); 
	setAttrib(RET, R_NamesSymbol, RET_NAMES);

	/* Set INFO and Copy A -> C. */
	INTEGER(INFO)[0] = 0;
	pt_A = REAL(A);
	pt_C = REAL(C);
	for(i = 0; i < pt_CLDIM[0] * pt_CLDIM[1]; i++){
		*pt_C = *pt_A;
		pt_A++;
		pt_C++;
	}

	/* Call Fortran. */
        F77_CALL(rpdgetrf)(REAL(C), INTEGER(ICTXT), INTEGER(MYROW),
			INTEGER(MYCOL), INTEGER(DESCA), INTEGER(M), INTEGER(N),
			INTEGER(LIPIV), INTEGER(INFO));

	/* Return. */
	UNPROTECT(4);
        return(RET);
} /* End of R_PDGETRF(). */


SEXP R_PDPOTRF(SEXP A, SEXP CLDIM, SEXP ICTXT, SEXP MYROW, SEXP MYCOL,
		SEXP DESCA, SEXP N, SEXP UPLO){
	int i, *pt_CLDIM = INTEGER(CLDIM);
	double *pt_A, *pt_C;
	SEXP RET, RET_NAMES, INFO, C;

	/* Protect R objects. */
	PROTECT(RET = allocVector(VECSXP, 2));
	PROTECT(RET_NAMES = allocVector(STRSXP, 2));
	PROTECT(INFO = allocVector(INTSXP, 1));
	PROTECT(C = allocMatrix(REALSXP, pt_CLDIM[0], pt_CLDIM[1]));

	SET_VECTOR_ELT(RET, 0, INFO);
	SET_VECTOR_ELT(RET, 1, C);
	SET_STRING_ELT(RET_NAMES, 0, mkChar("info")); 
	SET_STRING_ELT(RET_NAMES, 1, mkChar("A")); 
	setAttrib(RET, R_NamesSymbol, RET_NAMES);

	/* Copy A -> C and set INFO and return R objects. */
	INTEGER(INFO)[0] = 0;
	pt_A = REAL(A);
	pt_C = REAL(C);
	for(i = 0; i < pt_CLDIM[0] * pt_CLDIM[1]; i++){
		*pt_C = *pt_A;
		pt_A++;
		pt_C++;
	}

	/* Call Fortran. */
        F77_CALL(rpdpotrf)(REAL(C), INTEGER(ICTXT), INTEGER(MYROW),
			INTEGER(MYCOL), INTEGER(DESCA), INTEGER(N),
			CHARPT(UPLO, 0), INTEGER(INFO));

	/* Return. */
	UNPROTECT(4);
        return(RET);
} /* End of R_PDPOTRF(). */
