/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

// Copyright 2013, Schmidt

// Pade' Coefficients for p==q==13

// Generated from bc via:
// scale=50
// define f(x) {
// if (x <= 1) return (1);
// return (x*f(x-1));
// }
// m=13
// for(j=0; j<=m; j++) {f(2*m-j)*f(m)/f(2*m)/f(m-j)/f(j);}
const double matexp_pade_coefs[14] =
{
  1.0,
  0.5,
  0.12,
  1.833333333333333333333e-2,
  1.992753623188405797101e-3,
  1.630434782608695652174e-4,
  1.035196687370600414079e-5,
  5.175983436853002070393e-7,
  2.043151356652500817261e-8,
  6.306022705717595115002e-10,
  1.483770048404140027059e-11,
  2.529153491597965955215e-13,
  2.810170546219962172461e-15,
  1.544049750670308885967e-17
};
