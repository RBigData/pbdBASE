// Copyright 2013, Schmidt

#include "Rtools.h"


void PRINT(SEXP x)
{
  int ptct = 0;
  SEXP basePackage;
  
  PT(basePackage, ptct);
  basePackage = eval( lang2( install("getNamespace"), ScalarString(mkChar("base")) ), R_GlobalEnv );
  
  eval( lang2( install("print"), x), basePackage);
  
  UNPT(ptct);
}


