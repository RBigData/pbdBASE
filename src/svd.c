/*  Copyright (c) 2012-2013, 2017 by Drew Schmidt <wrathematics@gmail.com>
    All rights reserved.
    
    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:
    
    1. Redistributions of source code must retain the above copyright notice,
    this list of conditions and the following disclaimer.
    
    2. Redistributions in binary form must reproduce the above copyright
    notice, this list of conditions and the following disclaimer in the
    documentation and/or other materials provided with the distribution.
    
    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
    PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
    CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
    PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
    PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
    LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "pbdBASE.h"


/* SVD */
SEXP R_PDGESVD(SEXP M, SEXP N, SEXP ASIZE, SEXP A, SEXP DESCA, 
    SEXP ULDIM, SEXP DESCU, SEXP VTLDIM, SEXP DESCVT, SEXP JOBU, SEXP JOBVT, 
    SEXP INPLACE)
{
  R_INIT;
  double *A_OUT;
  int IJ = 1, temp_lwork = -1;
  double temp_A = 0, temp_work = 0, *WORK;
  SEXP RET, RET_NAMES, INFO, D, U, VT;
  
  newRvec(INFO, 1, "int");
  newRvec(D, INT(ASIZE, 0), "dbl");
  newRmat(U, INT(ULDIM, 0), INT(ULDIM, 1), "dbl");
  newRmat(VT, INT(VTLDIM, 0), INT(VTLDIM, 1), "dbl");
  
  
  // Query size of workspace
  INT(INFO, 0) = 0;
  
  pdgesvd_(STR(JOBU, 0), STR(JOBVT, 0),
    INTP(M), INTP(N),
    &temp_A, &IJ, &IJ, INTP(DESCA),
    &temp_A, &temp_A, &IJ, &IJ, INTP(DESCU),
    &temp_A, &IJ, &IJ, INTP(DESCVT),
    &temp_work, &temp_lwork, INTP(INFO));
      
  // Allocate work vector and calculate svd
  temp_lwork = (int) temp_work;
  temp_lwork = nonzero(temp_lwork);
  
  WORK = (double *) R_alloc(temp_lwork, sizeof(double));
  
  A_OUT = (double *) R_alloc(nrows(A)*ncols(A), sizeof(double));
  memcpy(A_OUT, REAL(A), nrows(A)*ncols(A)*sizeof(double));
  
  pdgesvd_(STR(JOBU, 0), STR(JOBVT, 0),
    INTP(M), INTP(N),
    A_OUT, &IJ, &IJ, INTP(DESCA),
    REAL(D), REAL(U), &IJ, &IJ, INTP(DESCU),
    REAL(VT), &IJ, &IJ, INTP(DESCVT),
    WORK, &temp_lwork, INTP(INFO));
  
  // Manage return
  RET_NAMES = make_list_names(4, "info", "d", "u", "vt");
  RET = make_list(RET_NAMES, 4, INFO, D, U, VT);
  
  R_END;
  return RET;
}
