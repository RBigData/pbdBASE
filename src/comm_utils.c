/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

// Copyright 2014, 2016 Schmidt

#include <stdint.h>
#if (defined(__MINGW32__) || defined(__MINGW64__))
  #include <_mingw.h>
#endif

#include <mpi.h>
#include <Rinternals.h>


void comm_stop(char *msg)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  if (rank == 0)
    error(msg);
}



void comm_warning(char *msg)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  if (rank == 0)
    warning(msg);
}
