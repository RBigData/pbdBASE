base.rl2blas <- function(x, descx, vec, FUN)
{
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  if (!is.double(vec))
    storage.mode(vec) <- "double"
  
  ret <- .Call(R_RL2BLAS, 
               x, as.integer(dim(x)), as.integer(descx), vec, as.integer(length(vec)), as.integer(FUN))
  
  return(ret)
}



# matrix-vector insertion
base.rl2insert <- function(x, descx, vec, i, j)
{
  dim <- descx[3L:4L]
  
  if (i[1L] < 0){
    new <- 1L:dim[1L]
    i <- new[-which(new %in% abs(i))] # FIXME make this less stupid
  }
  
  if (j[1L] < 0){
    new <- 1L:dim[2L]
    j <- new[-which(new %in% abs(j))] # FIXME make this less stupid
  }
  
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  if (!is.double(vec))
    storage.mode(vec) <- "double"
  
  ret <- .Call(R_RL2INSERT, 
               x, as.integer(dim(x)), as.integer(descx), vec, as.integer(length(vec)), as.integer(i), as.integer(length(i)), as.integer(j), as.integer(length(j)))
  
  return( ret )
}



base.rcolcpy <- function(x, descx, y, descy, xcol, ycol)
{
  if (!is.double(x))
    storage.mode(x) <- "double"
  if (!is.double(y))
    storage.mode(y) <- "double"
  
  ret <- .Call(R_RCOLCPY, 
               x, as.integer(dim(x)), as.integer(descx), as.integer(xcol), y, as.integer(descy), as.integer(ycol), as.integer(length(ycol)))
  
  return( ret )
}



base.rcolcpy2 <- function(x, descx, y, descy, xcol, ycol)
{
  if (!is.double(x))
    storage.mode(x) <- "double"
  if (!is.double(y))
    storage.mode(y) <- "double"
  
  ret <- .Call(R_RCOLCPY2, 
               x, as.integer(dim(x)), as.integer(descx), as.integer(xcol), as.integer(length(xcol)), y, as.integer(descy), as.integer(ycol), as.integer(length(ycol)))
  
  return( ret )
}



base.rrowcpy <- function(x, descx, y, descy, xrow, yrow)
{
  if (!is.double(x))
    storage.mode(x) <- "double"
  if (!is.double(y))
    storage.mode(y) <- "double"
  
  ret <- .Call(R_RROWCPY, 
               x, as.integer(dim(x)), as.integer(descx), as.integer(xrow), y, as.integer(descy), as.integer(yrow), as.integer(length(yrow)))
  
  return( ret )
}



base.rrowcpy2 <- function(x, descx, y, descy, xrow, yrow)
{
  if (!is.double(x))
    storage.mode(x) <- "double"
  if (!is.double(y))
    storage.mode(y) <- "double"
  
  ret <- .Call(R_RROWCPY2, 
               x, as.integer(dim(x)), as.integer(descx), as.integer(xrow), as.integer(length(xrow)), y, as.integer(descy), as.integer(yrow), as.integer(length(yrow)))
  
  return( ret )
}



base.pdmvsum <- function(x, descx, y, descy)
{
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  if (!is.double(y))
    storage.mode(y) <- "double"
  
  ret <- .Call(R_PDMVSUM, 
               x, as.integer(dim(x)), as.integer(descx), y, as.integer(descy))
  
  return(ret)
}

