/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/* WCC: These functions are to export and access pointers in R.
 *
 * Wei-Chen Chen, Mar 2013.
 * Modified June 2015, Higgs and Schmidt.
 */

#include "pkg_global.h"


void set_BLACS_APTS_in_R(){
	/* Define R objects. */
	SEXP R_apts;

	/* Protect R objects. */
	PROTECT(R_apts = R_MakeExternalPtr(&BLACS_APTS, R_NilValue, R_NilValue));

	/* Assign an R object in ".GlobalEnv". */
	defineVar(install(BLACS_APTS_R_NAME), R_apts, R_GlobalEnv);

	/* These are only saw by new pakcages. */
	BLACS_APTS.BI_MaxNCtxt = &BI_MaxNCtxt;
	BLACS_APTS.BI_MaxNSysCtxt = &BI_MaxNSysCtxt;
	BLACS_APTS.BI_Iam = &BI_Iam;
	BLACS_APTS.BI_Np = &BI_Np;
	BLACS_APTS.BI_ReadyB = BI_ReadyB;
	BLACS_APTS.BI_ActiveQ = BI_ActiveQ;
	BLACS_APTS.BI_AuxBuff = &BI_AuxBuff;
	BLACS_APTS.BI_MyContxts = BI_MyContxts;
	BLACS_APTS.BI_SysContxts = BI_SysContxts;
	BLACS_APTS.BI_COMM_WORLD = BI_COMM_WORLD;
	BLACS_APTS.BI_Stats = BI_Stats;

	#if (BLACS_APTS_DEBUG & 1) == 1
	int myrank;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	if(myrank == 0){
		REprintf("  %s (v): %d %d %d %d %d.\n", __FILE__, BI_MaxNCtxt,
			BI_MaxNSysCtxt, BI_Iam, BI_Np, BI_AuxBuff);
/* Not a good idea to print NULL pointers.
		REprintf("  %s (v): %d %d %d %d.\n", __FILE__, *BI_ReadyB,
			*BI_ActiveQ, **BI_MyContxts, *BI_COMM_WORLD);
		REprintf("  %s (v): %d %d.\n", __FILE__, *BI_SysContxts,
			*BI_Stats);
*/
		REprintf("  %s (v): %d %d %d.\n", __FILE__, BI_AuxBuff.Len,
			BI_AuxBuff.nAops, BI_AuxBuff.N);
		REprintf("  %s (a): %x %x %x %x %x.\n", __FILE__, &BI_MaxNCtxt,
			&BI_MaxNSysCtxt, &BI_Iam, &BI_Np, &BI_AuxBuff);
		REprintf("  %s (a): %x %x %x %x.\n", __FILE__, BI_ReadyB,
			BI_ActiveQ, *BI_MyContxts, BI_COMM_WORLD);
		REprintf("  %s (a): %x %x.\n", __FILE__, BI_SysContxts,
			BI_Stats);
	}
	#endif

	/* Unprotect R objects. */
	UNPROTECT(1);
} /* End of set_BLACS_APTS_in_R(). */

void get_BLACS_APTS_from_R(){
        /* Define an R object. */
        SEXP R_apts;

        /* Get an R object from ".GlobalEnv". */
        R_apts = findVar(install(BLACS_APTS_R_NAME), R_GlobalEnv);
        if(R_apts == R_NilValue){
                error(".__BLACS_APTS__ is not found in .GlobalEnv");
        }
        /* Get pointers. */
        BLACS_APTS_ptr = R_ExternalPtrAddr(R_apts);

        /* These are only saw by "pbdMPI" not "Rmpi". */
	BI_MaxNCtxt = (int) *BLACS_APTS_ptr->BI_MaxNCtxt;
	BI_MaxNSysCtxt = (int) *BLACS_APTS_ptr->BI_MaxNSysCtxt;
	BI_Iam = (int) *BLACS_APTS_ptr->BI_Iam;
	BI_Np = (int) *BLACS_APTS_ptr->BI_Np;
	BI_ReadyB = (BLACBUFF*) BLACS_APTS_ptr->BI_ReadyB;
	BI_ActiveQ = (BLACBUFF*) BLACS_APTS_ptr->BI_ActiveQ;

	/* This is wrong, and SUNCC doesn't like this. */
	// BI_AuxBuff = (BLACBUFF) *BLACS_APTS_ptr->BI_AuxBuff;
	BI_AuxBuff.Buff = (char*) BLACS_APTS_ptr->BI_AuxBuff->Buff;
	BI_AuxBuff.Len = (int) BLACS_APTS_ptr->BI_AuxBuff->Len;
	BI_AuxBuff.nAops = (int) BLACS_APTS_ptr->BI_AuxBuff->nAops;
	BI_AuxBuff.Aops = (MPI_Request*) BLACS_APTS_ptr->BI_AuxBuff->Aops;
	BI_AuxBuff.dtype = (MPI_Datatype) BLACS_APTS_ptr->BI_AuxBuff->dtype;
	BI_AuxBuff.N = (int) BLACS_APTS_ptr->BI_AuxBuff->N;
	BI_AuxBuff.prev = (BLACBUFF*) BLACS_APTS_ptr->BI_AuxBuff->prev;
	BI_AuxBuff.next = (BLACBUFF*) BLACS_APTS_ptr->BI_AuxBuff->next;

	BI_MyContxts = (BLACSCONTEXT**) BLACS_APTS_ptr->BI_MyContxts;
	BI_SysContxts = (MPI_Comm*) BLACS_APTS_ptr->BI_SysContxts;
	BI_COMM_WORLD = (int*) BLACS_APTS_ptr->BI_COMM_WORLD;
	BI_Stats = (MPI_Status*) BLACS_APTS_ptr->BI_Stats;

	#if (BLACS_APTS_DEBUG & 1) == 1
	int myrank;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	if(myrank == 0){
		REprintf("  %s (v): %d %d %d %d %d.\n", __FILE__, BI_MaxNCtxt,
			BI_MaxNSysCtxt, BI_Iam, BI_Np, BI_AuxBuff);
		REprintf("  %s (v): %d %d %d.\n", __FILE__, BI_AuxBuff.Len,
			BI_AuxBuff.nAops, BI_AuxBuff.N);
		REprintf("  %s (a): %x %x %x %x %x.\n", __FILE__, &BI_MaxNCtxt,
			&BI_MaxNSysCtxt, &BI_Iam, &BI_Np, &BI_AuxBuff);
		REprintf("  %s (a): %x %x %x %x.\n", __FILE__, BI_ReadyB,
			BI_ActiveQ, *BI_MyContxts, BI_COMM_WORLD);
		REprintf("  %s (a): %x %x.\n", __FILE__, BI_SysContxts,
			BI_Stats);
	}
	#endif
} /* End of get_BLACS_APTS_from_R(). */
