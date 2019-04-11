/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

// Copyright 2016, Schmidt

#include "../../blacs.h"

void dmat_ldimget(int *desc, int* nrows, int* ncols)
{
  int M = desc[2];
  int N = desc[3];
  int Mb = desc[4];
  int Nb = desc[5];
  
  int ictxt = desc[1];
  int rsrc = desc[6];
  int csrc = desc[7];
  
  int nprow, npcol, myprow, mypcol;
  Cblacs_gridinfo(ictxt, &nprow, &npcol, &myprow, &mypcol);
  
  *nrows = numroc_(&M, &Mb, &myprow, &rsrc, &nprow);
  *ncols = numroc_(&N, &Nb, &mypcol, &csrc, &npcol);
  
  if (*nrows < 1 || *ncols < 1)
  {
    *nrows = 0;
    *ncols = 0;
  }
}



// Inspired by https://andyspiros.wordpress.com/2011/07/08/an-example-of-blacs-with-c/
/**
* Distribute a global matrix stored on rank 0 to all processors
* in the grid.
* 
* @param desc
* Descriptor array for desired output distributed matrix.
* @param 
*/
int dmat_as_ddmatrix(int *desc, double *A_global, double *A_local)
{
  const int M = desc[2], N = desc[3];
  const int Mb = desc[4], Nb = desc[5];
  int nrows, ncols;         // size of A_local
  
  int row, col;             // index over global matrix 
  int recvr = 0, recvc = 0; // local matrix index
  
  int ictxt = desc[1];      // blacs context
  int nprow, npcol;         // number process rows/cols
  int myprow, mypcol;       // current process row/col
  
  int rdest = 0, cdest = 0; // process grid row x col destination
  int nr, nc;               // message size; max Mb/Nb
  
  dmat_ldimget(desc, &nrows, &ncols);
  if (nrows < 1 || ncols < 1)
    return -1;
  
  Cblacs_gridinfo(ictxt, &nprow, &npcol, &myprow, &mypcol);
  const int owner = (myprow == 0 && mypcol == 0);
  
  
  for (row=0; row<M; row+=Mb, rdest=(rdest+1)%nprow)
  {
    cdest = 0;
    nr = Mb;
    // Is this the last row block?
    if (M - row < Mb)
      nr = M - row;
    
    for (col=0; col<N; col+=Nb, cdest=(cdest+1)%npcol)
    {
      nc = Nb;
      // Is this the last col block?
      if (N - col < Nb)
        nc = N - col;
    
      if (owner)
      {
        // Send a nr x nc submatrix to process (rdest, cdest)
        Cdgesd2d(ictxt, nr, nc, A_global + (row + M*col), M, rdest, cdest);
      }

      if (myprow == rdest && mypcol == cdest)
      {
        // Receive nr x nc submatrix to recvr x recvc in A_local
        Cdgerv2d(ictxt, nr, nc, A_local + (recvr + nrows*recvc), nrows, 0, 0);
        recvc = (recvc+nc)%ncols;
      }
    }
    
    if (myprow == rdest)
      recvr = (recvr+nr)%nrows;
  }
  
  return 0;
}
