base.matpow_by_squaring <- function(A, b=1)
{
  b <- as.integer(b)
  
  if (!is.double(A))
    storage.mode(A) <- "double"
  
  ret <- .Call("R_matpow_by_squaring", A, b, PACKAGE="pbdBASE")
  
  return( ret )
}



base.matexp_pade <- function(A)
{
  if (!is.double(A))
    storage.mode(A) <- "double"
  
  out <- .Call("R_matexp_pade", A, PACKAGE="pbdBASE")
  
  N <- out$N
  D <- out$D
  
  R <- solve(D) %*% N
  
  return( R )
}


