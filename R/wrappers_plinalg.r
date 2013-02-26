# ------------------------------------------------
# PDSVRK:  Symmetric Rank-k Update
# ------------------------------------------------

base.crossprod <- function(trans, x, descx, descc)
{
  trans <- toupper(trans)
  
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  cldim <- base.numroc(descc[3:4], descc[5:6], ICTXT=descc[2])
  
  ret <- .Call("R_PDCROSSPROD",
                  trans, x, as.integer(descx),
                  as.integer(cldim), as.integer(descc),
                  PACKAGE="pbdBASE")
  
  return( ret )
}

base.pdchtri <- function(x, descx, descc)
{
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  cldim <- base.numroc(descc[3:4], descc[5:6], ICTXT=descc[2])
  
  ret <- .Call("R_PDCHTRI", 
                x, as.integer(dim(x)), as.integer(descx), as.integer(cldim), as.integer(descc),
                PACKAGE="pbdBASE")
  
  return( ret )
}



