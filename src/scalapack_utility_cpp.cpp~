
#include <Rcpp.h>

  // ----------------------------------------------- //
 /* Global-to-local and local-to-global coordinates */
// ----------------------------------------------- //

void
  g2l_coord( std::vector<int> &ret,
            int i, int j, 
            Rcpp::IntegerVector &dim, Rcpp::IntegerVector &bldim,
            Rcpp::IntegerVector &procs, Rcpp::IntegerVector &src
           )
{
  // matrix block position
  ret[0] = i / (procs[0] * bldim[0]);
  ret[1] = j / (procs[1] * bldim[1]);
  
  // process grid block
  ret[2] = (src[0] + i / bldim[0]) % procs[0];
  ret[3] = (src[1] + j / bldim[1]) % procs[1];
  
  // local coordinates
  ret[4] = i % bldim[0] + bldim[0] * ret[0];
  ret[5] = j % bldim[1] + bldim[1] * ret[1];
}

void
  l2g_coord(std::vector<int> &ret,
            int i, int j, 
            Rcpp::IntegerVector &dim, Rcpp::IntegerVector &bldim,
            Rcpp::IntegerVector &procs, int myproc
           )
{
  const int nprocs = procs[0] * procs[1];
  ret[0] = nprocs*bldim[0] * (i-1)/bldim[0] + (i-1)%bldim[0] + ((nprocs+myproc)%nprocs)*bldim[0] + 1;
  ret[1] = nprocs*bldim[1] * (j-1)/bldim[1] + (j-1)%bldim[1] + ((nprocs+myproc)%nprocs)*bldim[1] + 1;
}

//////////////////////////////////////////////////////////////
//#include <Rcpp.h>

  // ----------------------------------------------- //
 /* Global-to-local and local-to-global coordinates */
// ----------------------------------------------- //

RcppExport SEXP rcpp_l2g_coord( SEXP ind_, 
                             SEXP dim_, SEXP bldim_,
                             SEXP procs_, SEXP myproc_
                             )
{
  Rcpp::IntegerVector ind(ind_);
  Rcpp::IntegerVector dim(dim_);
  Rcpp::IntegerVector bldim(bldim_);
  Rcpp::IntegerVector procs(procs_);
  int myproc = *(INTEGER(myproc_));

  std::vector<int> ret(2);

  l2g_coord(ret, ind[0], ind[1], dim, bldim, procs, myproc);
  
  return Rcpp::wrap( ret );
}

RcppExport SEXP rcpp_g2l_coord( SEXP ind_, 
                             SEXP dim_, SEXP bldim_,
                             SEXP procs_, SEXP src_
                             )
{
  Rcpp::IntegerVector ind(ind_);
  Rcpp::IntegerVector dim(dim_);
  Rcpp::IntegerVector bldim(bldim_);
  Rcpp::IntegerVector procs(procs_);
  Rcpp::IntegerVector src(src_);

  std::vector<int> ret(6);
 
  g2l_coord(ret, ind[0], ind[1], dim, bldim, procs, src);
  
  return Rcpp::wrap( ret );
}

  // --------------------------------------------- //
 /* create local subA from global A or vice versa */
// --------------------------------------------- //

// A --> subA
RcppExport SEXP block_submat(SEXP A_, 
                             SEXP dim_, SEXP ldim_, SEXP bldim_,
                             SEXP procs_, SEXP myproc_, SEXP src_
                             )
{
  Rcpp::NumericMatrix A(A_);
  Rcpp::IntegerVector dim(dim_);
  Rcpp::IntegerVector ldim(ldim_);
  Rcpp::IntegerVector bldim(bldim_);
  Rcpp::IntegerVector procs(procs_);
  Rcpp::IntegerVector myproc(myproc_);
  Rcpp::IntegerVector src(src_);

  Rcpp::NumericMatrix subA(ldim[0], ldim[1]);
  
  int i, j;
  std::vector<int> ret(6);

  // see http://www.netlib.org/scalapack/slug/node76.html#SECTION04432000000000000000
  for (j=0; j<dim[1]; j++){
    for (i=0; i<dim[0]; i++){
      g2l_coord(ret, i, j, dim, bldim, procs, src);
      if (myproc[0]==ret[2] && myproc[1]==ret[3]){ // fill subA
        subA(ret[4], ret[5]) = A(i,j);
      }
    }
  }
  
  return subA;
}

// subA --> A
RcppExport SEXP submat_to_gmat(SEXP subA_, 
                             SEXP dim_, SEXP bldim_,
                             SEXP procs_, SEXP myproc_, SEXP src_
                             )
{
  Rcpp::NumericMatrix subA(subA_);
  Rcpp::IntegerVector dim(dim_);
  Rcpp::IntegerVector bldim(bldim_);
  Rcpp::IntegerVector procs(procs_);
  Rcpp::IntegerVector myproc(myproc_);
  Rcpp::IntegerVector src(src_);
  
  Rcpp::NumericMatrix A(dim[0], dim[1]);

  int i, j;
  std::vector<int> ret(6);
  
  for (j=0; j<dim[1]; j++){
    for (i=0; i<dim[0]; i++){
      g2l_coord(ret, i, j, dim, bldim, procs, src);
      if (myproc[0]==ret[2] && myproc[1]==ret[3]){ // fill subA
        A(i,j) = subA(ret[4], ret[5]);
      }
    }
  }
  
  return A;
}


  // ----------------------------------------------------- //
 /* Fill lower triangle of A, distributed as subA, with 0 */
// ----------------------------------------------------- //

RcppExport SEXP fill_lower_zero(SEXP subA_, 
                             SEXP dim_, SEXP bldim_,
                             SEXP procs_, SEXP myproc_, SEXP src_
                             )
{
  Rcpp::NumericMatrix subA(subA_);
  Rcpp::IntegerVector dim(dim_);
  Rcpp::IntegerVector bldim(bldim_);
  Rcpp::IntegerVector procs(procs_);
  Rcpp::IntegerVector myproc(myproc_);
  Rcpp::IntegerVector src(src_);

  int i, j;
  std::vector<int> ret(6);
  
  for (j=0; j<dim[1]-1; j++){
    for (i=j+1; i<dim[0]; i++){
      g2l_coord(ret, i, j, dim, bldim, procs, src);
      if (myproc[0]==ret[2] && myproc[1]==ret[3]){ // fill subA
        subA(ret[4], ret[5]) = 0;
      }
    }
  }
  
  return subA;
}


RcppExport SEXP optimal_process_grid (SEXP nprocs_)
{
  const int nprocs = *(INTEGER(nprocs_));
  Rcpp::IntegerVector procs(2);

  const int root = (int)std::sqrt(nprocs);
  int test = 0;
  
  for (int i=0; i<root; i++){
    test = root - i;
    
    if (nprocs % test == 0)
      break;
  }
  
  procs[0] = test;
  procs[1] = nprocs/test;
  
  return procs;
}


  // ------------------------ //
 /* Matrix/vector operations */
// ------------------------- //

// Calculate a "mod" b, where b is numeric
double fp_modulo (double a, double b)
{
  double ret=a;

  if (b==0){
    return R_NaN;
  }  
  else if (b > 0){
    if (a >= 0)
      ret = std::fmod(a, b);
    else
      ret = b-(std::fmod(-a, b));
  } else {
    if (a==0)
      ret = 0;
    else if (a > 0)
      ret = -(-b - (std::fmod(a, -b)));
    else
      ret = - (std::fmod(-a, -b));
  }

  return ret;
}

RcppExport SEXP ddmatrix_vecops(SEXP subA_, SEXP dim_, SEXP bldim_,
                             SEXP vec_, SEXP nvec_,
                             SEXP procs_, SEXP myproc_, SEXP src_,
                             SEXP FUN_
                             )
{
  Rcpp::NumericMatrix subA(subA_);
  Rcpp::IntegerVector dim(dim_);
  Rcpp::IntegerVector bldim(bldim_);
  Rcpp::NumericVector vec(vec_);
  const int nvec = *(INTEGER(nvec_));
  Rcpp::IntegerVector procs(procs_);
  Rcpp::IntegerVector myproc(myproc_);
  Rcpp::IntegerVector src(src_);
  const int FUN = *(INTEGER(FUN_));

  Rcpp::NumericMatrix cp(subA.nrow(), subA.ncol());

  int i, j, k;
  
  std::vector<int> ret(6);

  for (j=0; j<dim[1]; j++){
    for (i=0; i<dim[0]; i++){
      g2l_coord(ret, i, j, dim, bldim, procs, src);
      if (myproc[0]==ret[2] && myproc[1]==ret[3]){

        k = (i+j*dim[0]) % nvec;

        // arithmetic
        if (FUN==0)
          cp(ret[4], ret[5]) = subA(ret[4], ret[5]) + vec[k];
        else if (FUN==1)
          cp(ret[4], ret[5]) = subA(ret[4], ret[5]) - vec[k];
        else if (FUN==2)
          cp(ret[4], ret[5]) = subA(ret[4], ret[5]) * vec[k];
        else if (FUN==3)
          cp(ret[4], ret[5]) = subA(ret[4], ret[5]) / vec[k];
        else if (FUN==4)
          cp(ret[4], ret[5]) = std::pow(subA(ret[4], ret[5]), vec[k]);
        else if (FUN==5)
          cp(ret[4], ret[5]) = fp_modulo(subA(ret[4], ret[5]), vec[k]); // fp_modulo(subA(ret[4], ret[5]), vec[k]);
        else if (FUN==6)
          cp(ret[4], ret[5]) = fp_modulo(vec[k], subA(ret[4], ret[5]));
        // logical
        else if (FUN==7)
          cp(ret[4], ret[5]) = (double)(subA(ret[4], ret[5]) < vec[k]);
        else if (FUN==8)
          cp(ret[4], ret[5]) = (double)(subA(ret[4], ret[5]) > vec[k]);
        else if (FUN==9)
          cp(ret[4], ret[5]) = (double)(subA(ret[4], ret[5]) <= vec[k]);
        else if (FUN==10)
          cp(ret[4], ret[5]) = (double)(subA(ret[4], ret[5]) >= vec[k]);
        else if (FUN==11)
          cp(ret[4], ret[5]) = (double)(subA(ret[4], ret[5]) == vec[k]);
      }
    }
  }
  
  return cp;
}


// insertion
RcppExport SEXP ddmatrix_insert(SEXP subA_, SEXP dim_, SEXP bldim_,
                             SEXP vec_, SEXP nvec_,
                             SEXP I_, SEXP ni_, SEXP J_, SEXP nj_,
                             SEXP procs_, SEXP myproc_, SEXP src_
                             )
{
  Rcpp::NumericMatrix subA(subA_);
  Rcpp::IntegerVector dim(dim_);
  Rcpp::IntegerVector bldim(bldim_);
  Rcpp::NumericVector vec(vec_);
  const int nvec = *(INTEGER(nvec_));
  Rcpp::IntegerVector I(I_);
  const int ni = *(INTEGER(ni_));
  Rcpp::IntegerVector J(J_);
  const int nj = *(INTEGER(nj_));
  Rcpp::IntegerVector procs(procs_);
  Rcpp::IntegerVector myproc(myproc_);
  Rcpp::IntegerVector src(src_);

  int i, j, k;
  
  std::vector<int> ret(6);

  for (j=0; j<nj; j++){
    for (i=0; i<ni; i++){
      g2l_coord(ret, I[i]-1, J[j]-1, dim, bldim, procs, src);
      if (myproc[0]==ret[2] && myproc[1]==ret[3]){

        k = (i+j*dim[0]) % nvec;

        subA(ret[4], ret[5]) = vec[k];
      }
    }
  }
  
  return subA;
}
