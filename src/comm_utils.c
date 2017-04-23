/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

// Copyright 2014, 2016 Schmidt

#include "pbdBASE.h"


SEXP COMM_STOP(char *msg)
{
  SEXP mpiPackage;
  SEXP Rmsg;
  SEXP ret;
  
  PROTECT(mpiPackage = 
    eval( lang2( install("getNamespace"), ScalarString(mkChar("pbdMPI")) ), R_GlobalEnv )
  );
  
  PROTECT(Rmsg = allocVector(STRSXP, 1));
  SET_STRING_ELT(Rmsg, 0, mkChar(msg));
  
  ret = eval( lang2( install("comm.stop"), Rmsg ), mpiPackage );
  
  UNPROTECT(2);
  return ret;
}



SEXP COMM_WARNING(char *msg)
{
  SEXP mpiPackage;
  SEXP Rmsg;
  SEXP ret;
  
  PROTECT(mpiPackage = 
    eval( lang2( install("getNamespace"), ScalarString(mkChar("pbdMPI")) ), R_GlobalEnv )
  );
  
  PROTECT(Rmsg = allocVector(STRSXP, 1));
  SET_STRING_ELT(Rmsg, 0, mkChar(msg));
  
  ret = eval( lang2( install("comm.warning"), Rmsg ), mpiPackage );
  
  UNPROTECT(2);
  return ret;
}



SEXP COMM_PRINT(SEXP x)
{
  SEXP mpiPackage;
  SEXP ret;
  
  PROTECT(mpiPackage = 
    eval( lang2( install("getNamespace"), ScalarString(mkChar("pbdMPI")) ), R_GlobalEnv )
  );
  
  ret = eval( lang2( install("comm.print"), x ), mpiPackage );
  
  UNPROTECT(2);
  return ret;
}
