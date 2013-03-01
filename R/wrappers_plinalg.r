# crossprod/tcrossprod
base.crossprod <- function(uplo, trans, x, descx, descc)
{
  trans <- toupper(trans)
  uplo <- toupper(uplo)
  
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  cldim <- base.numroc(descc[3:4], descc[5:6], ICTXT=descc[2])
  
  ret <- .Call("R_PDCROSSPROD",
                  uplo, trans, x, as.integer(descx),
                  as.integer(cldim), as.integer(descc),
                  PACKAGE="pbdBASE")
  
  return( ret )
}


# FIXME move row/col adjustment down to Fortran (currently in calling DMAT chol2inv method)
base.pdchtri <- function(uplo, x, descx, descc)
{
  uplo <- toupper(uplo)
  
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  cldim <- base.numroc(descc[3:4], descc[5:6], ICTXT=descc[2])
  
  ret <- .Call("R_PDCHTRI", 
                uplo, x, as.integer(dim(x)), as.integer(descx), 
                as.integer(cldim), as.integer(descc),
                PACKAGE="pbdBASE")
  
  return( ret )
}



base.pdclvar <- function(x, descx)
{
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  ret <- .Call("R_PDCLVAR", 
              x, as.integer(descx), as.integer(dim(x)[2L]),
              PACKAGE="pbdBASE")
  
  return( ret )
}




