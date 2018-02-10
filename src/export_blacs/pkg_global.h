/* WCC: point structure to array pointers.
 *
 * Wei-Chen Chen, May 2012.
 */

#ifndef __SPMD_GLOBAL__
#define __SPMD_GLOBAL__

#include <R.h>
#include <Rinternals.h>
#include <mpi.h>

// pkg_Bdef.h & pkg_Bconf.h are copied from ScaLAPACK, but they should be
// exported by pbdSLAP eventually.
#include "pkg_Bdef.h"

#ifdef __GNUC__
#define VARIABLE_IS_NOT_USED __attribute__ ((unused))
#else
#define VARIABLE_IS_NOT_USED
#endif

/* For loading. */
#define BLACS_APTS_R_NAME	".__BLACS_APTS__"

/* For debugging. */
#define BLACS_APTS_DEBUG	0

/* Declared in "BI_GlobalVars.c", or similar main c functions. */
extern int BI_MaxNCtxt, BI_MaxNSysCtxt, BI_Iam, BI_Np;
extern BLACBUFF *BI_ReadyB, *BI_ActiveQ, BI_AuxBuff;
extern BLACSCONTEXT **BI_MyContxts;
extern MPI_Comm *BI_SysContxts;
extern int *BI_COMM_WORLD;
extern MPI_Status *BI_Stats;

/* Collections. */
typedef struct _blacs_array_pointers	blacs_array_pointers;
struct _blacs_array_pointers{
	int *BI_MaxNCtxt, *BI_MaxNSysCtxt, *BI_Iam, *BI_Np;
	BLACBUFF *BI_ReadyB, *BI_ActiveQ, *BI_AuxBuff;
	BLACSCONTEXT **BI_MyContxts;
	MPI_Comm *BI_SysContxts;
	int *BI_COMM_WORLD;
	MPI_Status *BI_Stats;
};
blacs_array_pointers BLACS_APTS, *BLACS_APTS_ptr;

/* In "pkg_tools.c". */
void set_BLACS_APTS_in_R();
void get_BLACS_APTS_from_R();

#endif

