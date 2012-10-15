#include <R.h>

void F77_SUB(bprnt_c)(const char *data_c, int *data_i, int *data_j,
		double *data_d){
	if(strlen(data_c) > 255){
		Rprintf("invalid name length in bprnt_c");
	} else{
		Rprintf("%s[%6d,%6d]=%30.18f\n",
			data_c, *data_i, *data_j, *data_d);
	}
} /* End of bprnt_c(). */

