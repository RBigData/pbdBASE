/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

// Copyright 2014, 2016 Schmidt

#include "pbdBASE.h"


SEXP COMM_STOP(char *msg)
{
  SEXP mpiPackage, fun_install, expr;
  SEXP Rmsg;
  SEXP ret;
  
  PROTECT(Rmsg = allocVector(STRSXP, 1));
  SET_STRING_ELT(Rmsg, 0, mkChar(msg));
  
  PROTECT(mpiPackage = evalfun_stringarg("getNamespace", "pbdMPI"));
  PROTECT(fun_install = install("comm.stop"));
  PROTECT(expr = lang2(fun_install, Rmsg));
  ret = eval(expr, mpiPackage);
  
  UNPROTECT(4);
  return ret;
}



SEXP COMM_WARNING(char *msg)
{
  SEXP mpiPackage, fun_install, expr;
  SEXP Rmsg;
  SEXP ret;
  
  PROTECT(Rmsg = allocVector(STRSXP, 1));
  SET_STRING_ELT(Rmsg, 0, mkChar(msg));
  
  PROTECT(mpiPackage = evalfun_stringarg("getNamespace", "pbdMPI"));
  PROTECT(fun_install = install("comm.warning"));
  PROTECT(expr = lang2(fun_install, Rmsg));
  ret = eval(expr, mpiPackage);
  
  UNPROTECT(4);
  return ret;
}



SEXP COMM_PRINT(SEXP x)
{
  SEXP mpiPackage, fun_install, expr;
  SEXP ret;
  
  PROTECT(mpiPackage = evalfun_stringarg("getNamespace", "pbdMPI"));
  PROTECT(fun_install = install("comm.print"));
  PROTECT(expr = lang2(fun_install, x));
  ret = eval(expr, mpiPackage);
  
  UNPROTECT(3);
  return ret;
}
