/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

// Copyright 2013, Chen


// For C/Fortran char* string lengths using size_t
#ifdef USE_FC_LEN_T
  #include <stddef.h>
  #include <Rconfig.h>    // this defines FC_LEN_T
  #include <string.h>
#endif

// R.h needs to be included after Rconfig.h
#include <R.h>

#ifdef FC_LEN_T
  void bprntc_(const char *data_c, int *data_i, int *data_j, double *data_d,
    FC_LEN_T data_c_len)
#else
  void bprntc_(const char *data_c, int *data_i, int *data_j, double *data_d)
#endif
{
  Rprintf("%s[%6d,%6d]=%30.18f\n", data_c, *data_i, *data_j, *data_d);
}

