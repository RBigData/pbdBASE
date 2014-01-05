base.p_matpow_by_squaring_wrap <- function(A, desca, b=1)
{
  b <- as.integer(b)
  desca <- as.integer(desca)
  
  if (!is.double(A))
    storage.mode(A) <- "double"
  
  ret <- .Call("R_p_matpow_by_squaring", A, desca, b, PACKAGE="pbdBASE")
  
  return( ret )
}



base.p_matexp_pade_wrap <- function(A, desca)
{
  desca <- as.integer(desca)
  
  if (!is.double(A))
    storage.mode(A) <- "double"
  
  ret <- .Call("R_p_matexp_pade", A, desca, PACKAGE="pbdBASE")
  
  return( ret )
}

