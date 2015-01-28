# Sums
base.igsum2d <- function(ICTXT, SCOPE, m, n, x, lda, RDEST, CDEST)
{
  if (!is.matrix(x) && !is.vector(x))
    comm.stop("ERROR : object 'x' must be of class matrix or vector.")
  
  if (!is.integer(x))
    storage.mode(x) <- "integer"
  
  out <- .Call(R_igsum2d1, as.integer(ICTXT), as.character(SCOPE), 
                as.integer(m), as.integer(n), x, as.integer(lda), 
                as.integer(RDEST), as.integer(CDEST))
  
  return( out )
}

base.dgsum2d <- function(ICTXT, SCOPE, m, n, x, lda, RDEST, CDEST)
{
  if (!is.matrix(x) && !is.vector(x))
    comm.stop("ERROR : object 'x' must be of class matrix or vector.")
  
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  out <- .Call(R_dgsum2d1, as.integer(ICTXT), as.character(SCOPE), 
                as.integer(m), as.integer(n), x, as.integer(lda), 
                as.integer(RDEST), as.integer(CDEST))
  
  return( out )
}


# Max value
base.igamx2d <- function(ICTXT, SCOPE, m, n, x, lda, RDEST, CDEST)
{
  if (!is.matrix(x) && !is.vector(x))
    comm.stop("ERROR : object 'x' must be of class matrix or vector.")
  
  if (!is.integer(x))
    storage.mode(x) <- "integer"
  
  out <- .Call(R_igamx2d1, as.integer(ICTXT), as.character(SCOPE), 
                as.integer(m), as.integer(n), x, as.integer(lda), 
                as.integer(RDEST), as.integer(CDEST))
  
  return( out )
}

base.dgamx2d <- function(ICTXT, SCOPE, m, n, x, lda, RDEST, CDEST)
{
  if (!is.matrix(x) && !is.vector(x))
    comm.stop("ERROR : object 'x' must be of class matrix or vector.")
  
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  out <- .Call(R_dgamx2d1, as.integer(ICTXT), as.character(SCOPE), 
                as.integer(m), as.integer(n), x, as.integer(lda), 
                as.integer(RDEST), as.integer(CDEST))
  
  return( out )
}


# Min value
base.igamn2d <- function(ICTXT, SCOPE, m, n, x, lda, RDEST, CDEST)
{
  if (!is.matrix(x) && !is.vector(x))
    comm.stop("ERROR : object 'x' must be of class matrix or vector.")
  
  if (!is.integer(x))
    storage.mode(x) <- "integer"
  
  out <- .Call(R_igamn2d1, as.integer(ICTXT), as.character(SCOPE), 
                as.integer(m), as.integer(n), x, as.integer(lda), 
                as.integer(RDEST), as.integer(CDEST))
  
  return( out )
}

base.dgamn2d <- function(ICTXT, SCOPE, m, n, x, lda, RDEST, CDEST)
{
  if (!is.matrix(x) && !is.vector(x))
    comm.stop("ERROR : object 'x' must be of class matrix or vector.")
  
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  out <- .Call(R_dgamn2d1, as.integer(ICTXT), as.character(SCOPE), 
                as.integer(m), as.integer(n), x, as.integer(lda), 
                as.integer(RDEST), as.integer(CDEST))
  
  return( out )
}


# point to point communication
base.dgesd2d <- function(ICTXT, SCOPE, m, n, x, lda, RDEST, CDEST)
{
  if (!is.matrix(x) && !is.vector(x))
    comm.stop("ERROR : object 'x' must be of class matrix or vector.")
  
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  out <- .Call(R_dgesd2d1, as.integer(ICTXT), as.integer(m), as.integer(n), 
                x, as.integer(lda), as.integer(RDEST), as.integer(CDEST))
  
  return( out )
}

base.dgerv2d <- function(ICTXT, SCOPE, m, n, x, lda, RDEST, CDEST)
{
  if (!is.matrix(x) && !is.vector(x))
    comm.stop("ERROR : object 'x' must be of class matrix or vector.")
  
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  out <- .Call(R_dgerv2d1, as.integer(ICTXT), as.integer(m), as.integer(n), 
                x, as.integer(lda), as.integer(RDEST), as.integer(CDEST))
  
  return( out )
}


