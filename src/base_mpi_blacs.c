/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

// Copyright 2012-2015, Schmidt
#if (defined(__MINGW32__) || defined(__MINGW64__))
  #include <_mingw.h>
#endif

#include <mpi.h>
#include "base/utils/utils.h"

// R.h and Rinternals.h needs to be included after Rconfig.h
#include "pbdBASE.h"
#include <RNACI.h>


// From pbdMPI pkg_global.h and pkg_tools.c
#define MPI_APTS_R_NAME		".__MPI_APTS__"

typedef struct _rmpi_array_pointers	rmpi_array_pointers;
struct _rmpi_array_pointers{
	MPI_Comm *comm;
	MPI_Status *status;
	MPI_Datatype *datatype;
	MPI_Info *info;
	MPI_Request *request;
	int *COMM_MAXSIZE;
	int *STATUS_MAXSIZE;
	int *REQUEST_MAXSIZE;
};

void Cblacs_get(int ConTxt, int what, int *val);
void Cblacs_gridinit(int *ConTxt, char *order, int nprow, int npcol);
void Cblacs_gridinfo(int ConTxt, int *nprow, int *npcol, int *myrow, int *mycol);
void Cblacs_exit(int NotDone);
void Cblacs_gridexit(int Contxt);
int Csys2blacs_handle(MPI_Comm Comm);
void Cfree_blacs_system_handle(int Shandle);

SEXP R_optimal_grid(SEXP NPROCS)
{
  R_INIT;
  SEXP NPROW, NPCOL, RET, RET_NAMES;
  
  newRvec(NPROW, 1, "int", TRUE);
  newRvec(NPCOL, 1, "int", TRUE);
  
  optimalgrid_(INTP(NPROCS), INTP(NPROW), INTP(NPCOL));
  
  make_list_names(RET_NAMES, 2, "nprow", "npcol");
  make_list(RET, RET_NAMES, 2, NPROW, NPCOL);
  
  R_END;
  return RET;
}

SEXP R_blacs_gridinit(SEXP NPROW_in, SEXP NPCOL_in, SEXP SHANDLE)
{
  R_INIT;
  SEXP NPROW, NPCOL, MYROW, MYCOL, RET, RET_NAMES, ICTXT;
  newRvec(NPROW, 1, "int");
  newRvec(NPCOL, 1, "int");
  newRvec(MYROW, 1, "int");
  newRvec(MYCOL, 1, "int");
  newRvec(ICTXT, 1, "int");
  
  INT(NPROW) = INT(NPROW_in);
  INT(NPCOL) = INT(NPCOL_in);
  INT(ICTXT) = INT(SHANDLE);
  
  char order = 'R';
  
  Cblacs_gridinit(INTP(ICTXT), &order, INT(NPROW), INT(NPCOL));
  
  Cblacs_gridinfo(INT(ICTXT), INTP(NPROW), INTP(NPCOL), INTP(MYROW), INTP(MYCOL));
  
  make_list_names(RET_NAMES, 5, "NPROW", "NPCOL", "ICTXT", "MYROW", "MYCOL");
  make_list(RET, RET_NAMES, 5, NPROW, NPCOL, ICTXT, MYROW, MYCOL);
  R_END;
  return(RET);
}

SEXP R_blacs_init(SEXP NPROW_in, SEXP NPCOL_in, SEXP ICTXT_in)
{
  R_INIT;
  SEXP SHANDLE;
  newRvec(SHANDLE, 1, "int");
  Cblacs_get(INT(ICTXT_in), 0, INTP(SHANDLE));
  R_END;
  return(R_blacs_gridinit(NPROW_in, NPCOL_in, SHANDLE));
}

SEXP R_blacs_exit(SEXP CONT)
{
  Cblacs_exit(INT(CONT));
  
  return R_NilValue;
}

SEXP R_blacs_gridexit(SEXP CONT)
{
  Cblacs_gridexit(INT(CONT));
  
  return R_NilValue;
}

SEXP R_sys2blacs_handle(SEXP COMM)
{
  SEXP new_ctxt = PROTECT(Rf_allocVector(INTSXP, 1));
  SEXP R_APTS = findVar(install(MPI_APTS_R_NAME), R_GlobalEnv);
  rmpi_array_pointers *MPI_APTS_ptr = (rmpi_array_pointers*) R_ExternalPtrAddr(R_APTS);
  MPI_Comm comm = MPI_APTS_ptr->comm[INT(COMM)];
  INT(new_ctxt) = Csys2blacs_handle(comm);
  UNPROTECT(1);
  return new_ctxt;
}

SEXP R_free_blacs_system_handle(SEXP SHANDLE)
{
  Cfree_blacs_system_handle(INT(SHANDLE));
  return R_NilValue;
}
