#' base.redist
#' 
#' Redistribute a matrix from rank 0 to all ranks in block cyclic
#' fashion.
#' 
#' @param A
#' Matrix.
#' @param descx
#' ScaLAPACK descriptor array.
#' 
#' @export
base.redist <- function(desc, A)
{
  if (!is.double(A))
    storage.mode(A) <- "double"
  
  .Call("R_redist", desc, A, PACKAGE="pbdBASE")
}
