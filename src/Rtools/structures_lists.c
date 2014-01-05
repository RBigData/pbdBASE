// Copyright 2013, Schmidt

#include "Rtools.h"


// Build lists
SEXP make_list_names(int n, ...)
{
  int i;
  char *tmp;
  SEXP R_list_names;
  va_list listPointer;
  
  PROTECT(R_list_names = allocVector(STRSXP, n));
  
  va_start(listPointer, n);
  
  for(i=0; i<n; i++)
  {
    tmp = va_arg(listPointer, char *);
    
    SET_STRING_ELT(R_list_names, i, mkChar(tmp));
  }
  
  va_end(listPointer);
  
  UNPROTECT(1);
  return R_list_names;
}



SEXP make_list(SEXP R_list_names, ...)
{
  int i;
  const int n = LENGTH(R_list_names);
  SEXP tmp, R_list;
  va_list listPointer;
  
  PROTECT(R_list = allocVector(VECSXP, n));
  
  va_start(listPointer, R_list_names);
  
  for(i=0; i<n; i++)
  {
    tmp = va_arg(listPointer, SEXP);
    
    SET_VECTOR_ELT(R_list, i, tmp);
  }
  
  va_end(listPointer);
  
  setAttrib(R_list, R_NamesSymbol, R_list_names);
  
  UNPROTECT(1);
  return R_list;
}



SEXP make_list_nonames(int n, ...)
{
  int i;
  SEXP tmp, R_list;
  va_list listPointer;
  
  PROTECT(R_list = allocVector(VECSXP, n));
  
  va_start(listPointer, n);
  
  for(i=0; i<n; i++)
  {
    tmp = va_arg(listPointer, SEXP);
    
    SET_VECTOR_ELT(R_list, i, tmp);
  }
  
  va_end(listPointer);
  
  UNPROTECT(1);
  return R_list;
}

