#' @export
base.pdclvar <- function(x, descx)
{
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  ret <- .Call(R_PDCLVAR, x, as.integer(descx), as.integer(dim(x)[2L]))
  
  return( ret )
}
