#include <R.h>
#include <Rinternals.h>
#include <stdarg.h>


// Build lists
SEXP make_list_names(int n, ...)
{
  int i;
  SEXP tmp, R_list_names;
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
  
  va_start(listPointer, n);
  
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



// Example usage
SEXP listtest()
{
  SEXP a, b;
  SEXP R_list, R_list_names;
  
  PROTECT(a = allocVector(INTSXP, 2));
  PROTECT(b = allocVector(REALSXP, 1));
  
  INTEGER(a)[0] = 1;
  INTEGER(a)[1] = 2;
  
  REAL(b)[0] = -10.10214;
  
  R_list_names = make_list_names(2, "a", "b");
  R_list = make_list(R_list_names, a, b);
  
  UNPROTECT(2);
  return R_list;
}
