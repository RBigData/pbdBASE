/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

// Copyright 2013, Chen

#include <R.h>

void bprntc_(const char *data_c, int *data_i, int *data_j, double *data_d)
{
  Rprintf("%s[%6d,%6d]=%30.18f\n", data_c, *data_i, *data_j, *data_d);
}
