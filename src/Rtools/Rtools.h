#ifndef __RTOOLS_H__
#define __RTOOLS_H__


#include <R.h>
#include <Rinternals.h>
#include <stdarg.h>

#include <string.h>


#define INT(x,i) (INTEGER(x)[i])
#define DBL(x,i) (REAL(x)[i])
#define STR(x,i) ((char*)CHAR(STRING_ELT(x,i)))

#define INTP(x) (INTEGER(x))
#define DBLP(x) (REAL(x))


#define PT(x,ct) PROTECT((x)); (ct)++
#define UNPT(ct) (UNPROTECT(ct))


// tools.c
SEXP make_list_names(int n, ...);
SEXP make_list(SEXP R_list_names, ...);
SEXP make_list_nonames(int n, ...);

// alloc.c
SEXP Rvecalloc(int n, char *type);
SEXP Rmatalloc(int m, int n, char *type);

// printing.c
void PRINT(SEXP x);

#endif
