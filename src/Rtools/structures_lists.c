/*  Copyright (c) 2014, Schmidt
    All rights reserved.
    
    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met: 
    
    1. Redistributions of source code must retain the above copyright notice, this
       list of conditions and the following disclaimer. 
    2. Redistributions in binary form must reproduce the above copyright notice,
       this list of conditions and the following disclaimer in the documentation
       and/or other materials provided with the distribution. 
    
    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
    ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
    WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
    ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
    (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
    LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
    ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
    
    The views and conclusions contained in the software and documentation are those
    of the authors and should not be interpreted as representing official policies, 
    either expressed or implied, of the FreeBSD Project.  */

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

