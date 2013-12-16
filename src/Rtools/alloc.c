// Copyright 2013, Schmidt

#include "Rtools.h"


SEXP Rvecalloc(int n, char *type)
{
  SEXP RET;
  
  if (strcmp(type, "vec") == 0)
    RET = allocVector(VECSXP, n);
  else if (strcmp(type, "int") == 0)
    RET = allocVector(INTSXP, n);
  else if (strcmp(type, "double") == 0 || strcmp(type, "dbl") == 0)
    RET = allocVector(REALSXP, n);
  else if (strcmp(type, "str") == 0 || strcmp(type, "char*") == 0)
    RET = allocVector(STRSXP, n);
  
  return RET;
}

SEXP Rmatalloc(int m, int n, char *type)
{
  SEXP RET;
  
  if (strcmp(type, "vec") == 0)
    RET = allocMatrix(VECSXP, m, n);
  else if (strcmp(type, "int") == 0)
    RET = allocMatrix(INTSXP, m, n);
  else if (strcmp(type, "double") == 0 || strcmp(type, "dbl") == 0)
    RET = allocMatrix(REALSXP, m, n);
  else if (strcmp(type, "str") == 0 || strcmp(type, "char*") == 0)
    RET = allocMatrix(STRSXP, m, n);
  
  return RET;
}

