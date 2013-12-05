#ifndef __BASE_SCALAPACK_H__
#define __BASE_SCALAPACK_H__


#define TRUE 1
#define FALSE 0


#define NONZERO(x) (x?x:1)

#ifndef USING_R
  #define PRINT printf
#else
  #define PRINT Rprintf
#endif


typedef struct
{
  double *Data;
  int *dim;
  int lda;
} matrix;

// Return pointers
#define DATA(x) ((x)->Data)
#define DIM(x) ((x)->dim)
// Returns value
#define NROWS(x) ((x)->dim[0])
#define NCOLS(x) ((x)->dim[1])
#define LDA(x) ((x)->lda)



typedef struct
{
  double *Data;
  int *desc;
  int *ldim;
} ddmatrix;


// Return pointers
#define P_DATA(x) ((x)->Data)
#define P_DESC(x) ((x)->desc)
#define P_DIM(x) ((x)->(*desc+2)) // 2:3
#define P_LDIM(x) ((x)->ldim)
#define P_BLDIM(x) ((x)->desc+4) // 4:5
// Returns value
#define P_NROWS(x) ((x)->desc[2])
#define P_NCOLS(x) ((x)->desc[3])
#define P_CTXT(x) ((x)->desc[1])


#endif
