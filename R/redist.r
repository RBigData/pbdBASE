#' base.redist
#' 
#' Redistribute a matrix from rank 0 to all ranks in block cyclic
#' fashion.
#' 
#' @param desc
#' ScaLAPACK descriptor array.
#' @param A
#' Matrix.
#' 
#' @return A block cyclic matrix of the input matrix A from rank 0.
#' 
#' @examples
#' spmd.code <- "
#'   suppressMessages(library(pbdMPI))
#'   suppressMessages(library(pbdBASE))
#'   init.grid()
#'
#'   ### Set data matrix A and desc.
#'   A <- matrix(as.double(1:30), nrow = 6, ncol = 5)
#'   if (comm.rank() != 0)
#'     A <- matrix(as.double(0), nrow = 6, ncol = 5)
#'   dim <- dim(A)
#'   bldim <- c(3L, 3L)
#'   ldim <- base.numroc(dim = dim, bldim = bldim)
#'   desc <- base.descinit(dim = dim, bldim = bldim, ldim = ldim)
#'
#'   ### Redistribute from rank 0.
#'   dA <- base.redist(desc, A)
#'   comm.print(dA, all.rank = TRUE)
#'
#'   finalize()
#' "
#' pbdMPI::execmpi(spmd.code = spmd.code, nranks = 2L)
#' 
#' @useDynLib pbdBASE R_redist
#' @export
base.redist <- function(desc, A)
{
  if (!is.double(A))
    storage.mode(A) <- "double"
  
  .Call("R_redist", desc, A, PACKAGE="pbdBASE")
}
