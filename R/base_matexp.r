base.matexp <- function(A, p=6, t=1)
{
  if (!is.double(A))
    storage.mode(A) <- "double"
  
  R <- .Call(R_matexp, A, as.integer(p), as.double(t))
  
  return( R )
}

